#pragma once
#include <string>
#include <utility>
#include <opencv2/core.hpp>

struct GeoTransform {
    double origin_lon;   // долгота левого верхнего пикселя
    double origin_lat;   // широта левого верхнего пикселя
    double pixel_w;      // размер пикселя по долготе (градусы, > 0)
    double pixel_h;      // размер пикселя по широте (градусы, < 0 — растёт вниз)
};

struct DemData {
    cv::Mat elev;         // float32, высоты в метрах
    GeoTransform gt;
};

class DemLoader {
public:
    // Загрузить GeoTIFF или HGT, вернуть float32 cv::Mat + геотрансформацию
    static DemData load(const std::string& path);

    // Пересчёт (lat, lon) → пиксельные координаты (может быть дробным)
    static std::pair<double,double> geo_to_pixel(const GeoTransform& gt, double lat, double lon);

    // Билинейная интерполяция высоты в (lat, lon)
    static float sample(const DemData& dem, double lat, double lon, float nodata = 0.0f);

    // Размеры в метрах на пиксель (приближение для средней широты)
    static double meters_per_pixel_x(const GeoTransform& gt, double lat);
    static double meters_per_pixel_y(const GeoTransform& gt);
};
