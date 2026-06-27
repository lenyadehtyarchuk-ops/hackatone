#include "ridge_map.hpp"
#include <opencv2/imgproc.hpp>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric>

static constexpr float PI = static_cast<float>(M_PI);

// ──────────────────────────────────────────────────────────────────────────────
// Hessian eigenvalue computation
// H = [[dxx, dxy], [dxy, dyy]]  (x = col = East, y = row = South)
// λ1 ≤ λ2 (smaller / larger eigenvalue)
// ──────────────────────────────────────────────────────────────────────────────
static inline void hessian_eigen(float dxx, float dxy, float dyy,
                                  float& lam1, float& lam2,
                                  float& dir_deg_out) {
    float T    = dxx + dyy;
    float disc = std::sqrt(std::max(0.0f, T * T * 0.25f - (dxx * dyy - dxy * dxy)));
    lam1 = T * 0.5f - disc;  // smaller eigenvalue
    lam2 = T * 0.5f + disc;  // larger eigenvalue

    // Ridge direction = eigenvector of λ2 (the "flat" direction along the ridge).
    // In image coords: eigvec(λ2) = [dxy, λ2 - dxx].
    // x = col direction = East; y = row direction = South → North = -y.
    // Geographic bearing of ridge:
    //   East component  = dxy
    //   North component = -(λ2 - dxx) = dxx - λ2
    float ex = dxy, ey_north = dxx - lam2;
    float angle_deg = std::atan2(ex, ey_north) * 180.0f / PI;  // atan2(E, N)
    // Ridge is bidirectional → map to [0, 180)
    dir_deg_out = std::fmod(angle_deg + 360.0f, 180.0f);
}

// ──────────────────────────────────────────────────────────────────────────────
// Bilinear interpolation of a float Mat
// ──────────────────────────────────────────────────────────────────────────────
static float interp(const cv::Mat& m, float r, float c) {
    int r0 = static_cast<int>(r), c0 = static_cast<int>(c);
    if (r0 < 0 || c0 < 0 || r0 >= m.rows - 1 || c0 >= m.cols - 1)
        return 0.0f;
    float dr = r - r0, dc = c - c0;
    return m.at<float>(r0, c0)   * (1-dr)*(1-dc)
         + m.at<float>(r0, c0+1) * (1-dr)* dc
         + m.at<float>(r0+1, c0) *  dr   *(1-dc)
         + m.at<float>(r0+1,c0+1)*  dr   * dc;
}

// ──────────────────────────────────────────────────────────────────────────────
// RidgeMap constructor
// ──────────────────────────────────────────────────────────────────────────────
RidgeMap::RidgeMap(const DemData& dem,
                   float gauss_sigma_px,
                   float curv_thresh_rel,
                   float ratio_thresh) {
    // 1. Smooth DEM
    cv::Mat sm;
    cv::GaussianBlur(dem.elev, sm, cv::Size(0, 0), gauss_sigma_px, gauss_sigma_px,
                     cv::BORDER_REFLECT);

    // 2. Second derivatives via Sobel (ksize=5 for better accuracy)
    cv::Mat dx, dy;
    cv::Sobel(sm, dx, CV_32F, 1, 0, 5, 1.0, 0.0, cv::BORDER_REFLECT);
    cv::Sobel(sm, dy, CV_32F, 0, 1, 5, 1.0, 0.0, cv::BORDER_REFLECT);

    cv::Mat dxx, dyy, dxy;
    cv::Sobel(dx, dxx, CV_32F, 1, 0, 5, 1.0, 0.0, cv::BORDER_REFLECT);
    cv::Sobel(dy, dyy, CV_32F, 0, 1, 5, 1.0, 0.0, cv::BORDER_REFLECT);
    cv::Sobel(dx, dxy, CV_32F, 0, 1, 5, 1.0, 0.0, cv::BORDER_REFLECT);

    int R = sm.rows, C = sm.cols;

    // 3. Compute eigenvalues + direction at every pixel
    cv::Mat lam1_map(R, C, CV_32F);
    cv::Mat lam2_map(R, C, CV_32F);
    cv::Mat dir_map(R, C, CV_32F);

    float max_ridge_curv = 0.0f, max_valley_curv = 0.0f;

    for (int r = 0; r < R; ++r) {
        const float* dxxr = dxx.ptr<float>(r);
        const float* dyyr = dyy.ptr<float>(r);
        const float* dxyr = dxy.ptr<float>(r);
        float* l1 = lam1_map.ptr<float>(r);
        float* l2 = lam2_map.ptr<float>(r);
        float* dr = dir_map.ptr<float>(r);
        for (int c = 0; c < C; ++c) {
            hessian_eigen(dxxr[c], dxyr[c], dyyr[c], l1[c], l2[c], dr[c]);
            max_ridge_curv  = std::max(max_ridge_curv,  -l1[c]);
            max_valley_curv = std::max(max_valley_curv,  l2[c]);
        }
    }

    float ridge_thresh  = curv_thresh_rel * max_ridge_curv;
    float valley_thresh = curv_thresh_rel * max_valley_curv;

    std::cerr << "[RDG] max ridge curv=" << max_ridge_curv
              << " thresh=" << ridge_thresh
              << "  max valley curv=" << max_valley_curv
              << " thresh=" << valley_thresh << "\n";

    // 4. NMS: keep only local maxima along the cross-ridge (perpendicular) direction.
    //    Perpendicular to ridge = direction of eigvec(λ1) = [-dxy, dxx-λ1] (image coords).
    pts_.reserve(65536);
    int border = 3;

    for (int r = border; r < R - border; ++r) {
        const float* l1 = lam1_map.ptr<float>(r);
        const float* l2 = lam2_map.ptr<float>(r);
        const float* dr = dir_map.ptr<float>(r);
        const float* dxxr = dxx.ptr<float>(r);
        const float* dxyr = dxy.ptr<float>(r);
        const float* hv  = dem.elev.ptr<float>(r);
        for (int c = border; c < C - border; ++c) {
            if (hv[c] < 1.0f) continue;  // nodata

            float l1v = l1[c], l2v = l2[c];
            bool is_ridge  = (-l1v > ridge_thresh)  && (-l1v > -l2v * (1.0f / ratio_thresh));
            bool is_valley = ( l2v > valley_thresh) && ( l2v >  l1v * (1.0f / ratio_thresh));
            if (!is_ridge && !is_valley) continue;

            // Direction across the feature (perpendicular):
            // eigvec(λ1) = [-dxy, dxx-λ1] in image(col, row) coords
            float ex_perp = -dxyr[c];
            float ey_perp = dxxr[c] - l1v;
            float plen = std::sqrt(ex_perp*ex_perp + ey_perp*ey_perp);
            if (plen < 1e-6f) continue;
            ex_perp /= plen; ey_perp /= plen;

            // Curvature map for NMS
            const cv::Mat& curv_map = is_ridge ? lam1_map : lam2_map;
            float cv = is_ridge ? -l1v : l2v;  // curvature magnitude at (r,c)

            // Sample neighbors at ±1 step along perpendicular
            float r1 = r + ey_perp, c1 = static_cast<float>(c) + ex_perp;
            float r2 = r - ey_perp, c2 = static_cast<float>(c) - ex_perp;
            float cv1 = is_ridge ? -interp(curv_map, r1, c1) : interp(curv_map, r1, c1);
            float cv2 = is_ridge ? -interp(curv_map, r2, c2) : interp(curv_map, r2, c2);

            // NMS: suppress if not the local maximum
            if (cv < cv1 || cv < cv2) continue;

            RidgePoint pt;
            pt.lat       = static_cast<float>(dem.gt.origin_lat + r * dem.gt.pixel_h);
            pt.lon       = static_cast<float>(dem.gt.origin_lon + c * dem.gt.pixel_w);
            pt.elev      = hv[c];
            pt.type      = is_ridge ? FeatureType::RIDGE : FeatureType::VALLEY;
            pt.dir_deg   = dr[c];
            pt.curvature = cv;
            pts_.push_back(pt);
        }
    }

    std::cerr << "[RDG] Линии рельефа: " << pts_.size()
              << " точек (хребты + впадины)\n";

    // Compute DEM extent for grid
    lat0_ = static_cast<float>(dem.gt.origin_lat + R * dem.gt.pixel_h);
    lon0_ = static_cast<float>(dem.gt.origin_lon);
    build_grid(0.02f);  // ~2.2 km cells
}

void RidgeMap::build_grid(float cell_deg) {
    cell_deg_ = cell_deg;
    float lat_max = lat0_, lon_max = lon0_;
    for (const auto& p : pts_) {
        lat_max = std::max(lat_max, p.lat);
        lon_max = std::max(lon_max, p.lon);
    }
    grows_ = static_cast<int>((lat_max - lat0_) / cell_deg_) + 2;
    gcols_ = static_cast<int>((lon_max - lon0_) / cell_deg_) + 2;
    grid_.assign(grows_ * gcols_, {});

    for (uint32_t i = 0; i < static_cast<uint32_t>(pts_.size()); ++i) {
        int gr = std::clamp(static_cast<int>((pts_[i].lat - lat0_) / cell_deg_), 0, grows_-1);
        int gc = std::clamp(static_cast<int>((pts_[i].lon - lon0_) / cell_deg_), 0, gcols_-1);
        grid_[gr * gcols_ + gc].push_back(i);
    }
}

std::vector<RidgeMap::Match> RidgeMap::query(
        float lat, float lon, float radius_m,
        float target_elev, float elev_tol,
        float uav_hdg_deg, float hdg_tol_deg,
        FeatureType type) const {

    float r_deg = radius_m / 111320.0f;
    int r0 = std::max(0, static_cast<int>((lat - r_deg - lat0_) / cell_deg_));
    int r1 = std::min(grows_-1, static_cast<int>((lat + r_deg - lat0_) / cell_deg_));
    int c0 = std::max(0, static_cast<int>((lon - r_deg - lon0_) / cell_deg_));
    int c1 = std::min(gcols_-1, static_cast<int>((lon + r_deg - lon0_) / cell_deg_));
    float r2_deg = r_deg * r_deg;

    float hdg_rad = uav_hdg_deg * PI / 180.0f;
    float cos_lat = std::cos(lat * PI / 180.0f);

    std::vector<Match> result;
    for (int gr = r0; gr <= r1; ++gr) {
        for (int gc = c0; gc <= c1; ++gc) {
            for (uint32_t idx : grid_[gr * gcols_ + gc]) {
                const RidgePoint& pt = pts_[idx];
                if (pt.type != type) continue;

                float elev_err = std::abs(pt.elev - target_elev);
                if (elev_err > elev_tol) continue;

                float dlat = pt.lat - lat;
                float dlon = (pt.lon - lon) * cos_lat;
                if (dlat*dlat + dlon*dlon > r2_deg) continue;
                float dist_m = std::sqrt(dlat*dlat + dlon*dlon) * 111320.0f;

                // Heading compatibility: UAV heading should be ≈ perpendicular to ridge dir.
                // Cross-product |sin(hdg - ridge_dir)| should be large (close to ±1).
                float ridge_rad = pt.dir_deg * PI / 180.0f;
                float cross = std::sin(hdg_rad - ridge_rad);  // ≈ ±1 for perpendicular
                float hdg_err = (std::acos(std::min(1.0f, std::abs(cross))) * 180.0f / PI);
                // hdg_err = 0 means perpendicular, 90 means parallel
                if (hdg_err > hdg_tol_deg) continue;

                float score = std::exp(-elev_err * elev_err / (2.0f * 30.0f * 30.0f))
                            * std::exp(-dist_m   * dist_m   / (2.0f * 2000.0f * 2000.0f))
                            * std::abs(cross);   // bonus for perpendicular crossings
                result.push_back({&pt, dist_m, elev_err, hdg_err, score});
            }
        }
    }
    std::sort(result.begin(), result.end(),
              [](const Match& a, const Match& b){ return a.score > b.score; });
    return result;
}

cv::Mat RidgeMap::make_overlay(const DemData& dem, int /*scale*/) const {
    int R = dem.elev.rows, C = dem.elev.cols;
    cv::Mat out(R, C, CV_8UC3, cv::Scalar(0, 0, 0));

    auto geo_to_px = [&](float lat, float lon) -> cv::Point {
        int col = static_cast<int>((lon - dem.gt.origin_lon) / dem.gt.pixel_w);
        int row = static_cast<int>((lat - dem.gt.origin_lat) / dem.gt.pixel_h);
        return {col, row};
    };

    for (const auto& pt : pts_) {
        auto px = geo_to_px(pt.lat, pt.lon);
        if (px.x < 0 || px.y < 0 || px.x >= C || px.y >= R) continue;
        cv::Vec3b& pix = out.at<cv::Vec3b>(px.y, px.x);
        if (pt.type == FeatureType::RIDGE)
            pix = {0, 80, 255};   // orange-red for ridges
        else
            pix = {255, 80, 0};   // blue for valleys
    }
    return out;
}
