#include "keypoint_db.hpp"
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <cmath>
#include <algorithm>

KeypointDB::KeypointDB(const DemData& dem, int window, float min_prominence) {
    int W = 2 * window + 1;
    cv::Mat se = cv::getStructuringElement(cv::MORPH_RECT, {W, W});

    cv::Mat dil, ero;
    cv::dilate(dem.elev, dil, se);
    cv::erode(dem.elev,  ero,  se);

    int rows = dem.elev.rows, cols = dem.elev.cols;
    kps_.reserve(16384);

    for (int r = window; r < rows - window; ++r) {
        const float* h = dem.elev.ptr<float>(r);
        const float* d = dil.ptr<float>(r);
        const float* e = ero.ptr<float>(r);
        for (int c = window; c < cols - window; ++c) {
            float hv = h[c];
            if (hv < 1.0f) continue;

            float prominence = d[c] - e[c];
            if (prominence < min_prominence) continue;

            bool is_peak = (d[c] == hv);
            bool is_pit  = (e[c] == hv);
            if (!is_peak && !is_pit) continue;

            Keypoint kp;
            kp.lat  = static_cast<float>(dem.gt.origin_lat + r * dem.gt.pixel_h);
            kp.lon  = static_cast<float>(dem.gt.origin_lon + c * dem.gt.pixel_w);
            kp.elev = hv;
            kp.type = is_peak ? KpType::PEAK : KpType::PIT;
            kps_.push_back(kp);
        }
    }

    std::cerr << "[KPD] Ключевые точки: " << kps_.size()
              << " (peaks+pits, W=" << W << " min_prom=" << min_prominence << " м)\n";

    // DEM extent (lat increases upward, pixel_h < 0)
    lat0_ = static_cast<float>(dem.gt.origin_lat + rows * dem.gt.pixel_h);
    lon0_ = static_cast<float>(dem.gt.origin_lon);

    build_grid(0.05f);  // ~5.5 km cells
}

void KeypointDB::build_grid(float cell_deg) {
    cell_deg_ = cell_deg;

    float lat_max = lat0_;
    float lon_max = lon0_;
    for (const auto& kp : kps_) {
        lat_max = std::max(lat_max, kp.lat);
        lon_max = std::max(lon_max, kp.lon);
    }

    grows_ = static_cast<int>((lat_max - lat0_) / cell_deg_) + 2;
    gcols_ = static_cast<int>((lon_max - lon0_) / cell_deg_) + 2;
    grid_.assign(grows_ * gcols_, {});

    for (uint32_t i = 0; i < static_cast<uint32_t>(kps_.size()); ++i) {
        int gr = static_cast<int>((kps_[i].lat - lat0_) / cell_deg_);
        int gc = static_cast<int>((kps_[i].lon - lon0_) / cell_deg_);
        gr = std::clamp(gr, 0, grows_ - 1);
        gc = std::clamp(gc, 0, gcols_ - 1);
        grid_[gr * gcols_ + gc].push_back(i);
    }
}

std::vector<const Keypoint*> KeypointDB::query(float lat, float lon,
                                                 float target_elev, float elev_tol,
                                                 float radius_m, KpType type) const {
    float r_deg = radius_m / 111320.0f;
    int r0 = std::max(0,          static_cast<int>((lat - r_deg - lat0_) / cell_deg_));
    int r1 = std::min(grows_ - 1, static_cast<int>((lat + r_deg - lat0_) / cell_deg_));
    int c0 = std::max(0,          static_cast<int>((lon - r_deg - lon0_) / cell_deg_));
    int c1 = std::min(gcols_ - 1, static_cast<int>((lon + r_deg - lon0_) / cell_deg_));

    float r2 = r_deg * r_deg;
    std::vector<const Keypoint*> result;

    for (int gr = r0; gr <= r1; ++gr) {
        for (int gc = c0; gc <= c1; ++gc) {
            for (uint32_t idx : grid_[gr * gcols_ + gc]) {
                const Keypoint& kp = kps_[idx];
                if (kp.type != type) continue;
                if (std::abs(kp.elev - target_elev) > elev_tol) continue;
                float dlat = kp.lat - lat, dlon = kp.lon - lon;
                if (dlat * dlat + dlon * dlon > r2) continue;
                result.push_back(&kp);
            }
        }
    }
    return result;
}
