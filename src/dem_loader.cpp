#include "dem_loader.hpp"
#include <gdal_priv.h>
#include <stdexcept>
#include <cmath>
#include <iostream>

static void gdal_init() {
    static bool done = false;
    if (!done) { GDALAllRegister(); done = true; }
}

DemData DemLoader::load(const std::string& path) {
    gdal_init();

    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(path.c_str(), GA_ReadOnly));
    if (!ds)
        throw std::runtime_error("Cannot open DEM: " + path);

    double gt_arr[6];
    ds->GetGeoTransform(gt_arr);

    GeoTransform gt;
    gt.origin_lon = gt_arr[0];
    gt.origin_lat = gt_arr[3];
    gt.pixel_w    = gt_arr[1];
    gt.pixel_h    = gt_arr[5];

    int nx = ds->GetRasterXSize();
    int ny = ds->GetRasterYSize();

    GDALRasterBand* band = ds->GetRasterBand(1);
    double nodata_val = band->GetNoDataValue();

    cv::Mat elev(ny, nx, CV_32F);
    CPLErr err = band->RasterIO(GF_Read, 0, 0, nx, ny,
                                elev.ptr<float>(), nx, ny,
                                GDT_Float32, 0, 0);
    GDALClose(ds);

    if (err != CE_None)
        throw std::runtime_error("RasterIO failed for: " + path);

    // Заменить nodata на 0
    float nd = static_cast<float>(nodata_val);
    for (auto it = elev.begin<float>(); it != elev.end<float>(); ++it)
        if (*it == nd || std::isnan(*it)) *it = 0.0f;

    std::cerr << "[DEM] Загружено " << nx << "×" << ny
              << " пикселей из " << path << "\n";
    std::cerr << "[DEM] Геотрансформ: origin=("
              << gt.origin_lon << ", " << gt.origin_lat << ") "
              << "pixel=(" << gt.pixel_w << ", " << gt.pixel_h << ")\n";

    return DemData{elev, gt};
}

std::pair<double,double> DemLoader::geo_to_pixel(const GeoTransform& gt, double lat, double lon) {
    double px = (lon - gt.origin_lon) / gt.pixel_w;
    double py = (lat - gt.origin_lat) / gt.pixel_h;
    return {px, py};
}

float DemLoader::sample(const DemData& dem, double lat, double lon, float nodata) {
    auto [px, py] = geo_to_pixel(dem.gt, lat, lon);
    int x0 = static_cast<int>(px);
    int y0 = static_cast<int>(py);
    int nx = dem.elev.cols;
    int ny = dem.elev.rows;
    if (x0 < 0 || y0 < 0 || x0 >= nx - 1 || y0 >= ny - 1)
        return nodata;
    float dx = static_cast<float>(px - x0);
    float dy = static_cast<float>(py - y0);
    float v00 = dem.elev.at<float>(y0,     x0);
    float v10 = dem.elev.at<float>(y0,     x0 + 1);
    float v01 = dem.elev.at<float>(y0 + 1, x0);
    float v11 = dem.elev.at<float>(y0 + 1, x0 + 1);
    return v00*(1-dx)*(1-dy) + v10*dx*(1-dy) + v01*(1-dx)*dy + v11*dx*dy;
}

double DemLoader::meters_per_pixel_x(const GeoTransform& gt, double lat) {
    constexpr double R = 6378137.0;
    return std::abs(gt.pixel_w) * (M_PI / 180.0) * R * std::cos(lat * M_PI / 180.0);
}

double DemLoader::meters_per_pixel_y(const GeoTransform& gt) {
    constexpr double R = 6378137.0;
    return std::abs(gt.pixel_h) * (M_PI / 180.0) * R;
}
