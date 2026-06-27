#include "contour_pf.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

static constexpr float DEG2RAD = static_cast<float>(M_PI / 180.0);

ContourPF::ContourPF(const DemData& dem, const KeypointDB& kpdb,
                     double start_lat, double start_lon,
                     double heading_deg, double speed_mps,
                     int n, double pos_sigma_m, double hdg_sigma_deg,
                     const RidgeMap* rm,
                     double yaw_rate_dps, double yaw_sigma_dps)
    : dem_(dem), kpdb_(kpdb), rm_(rm), rng_(42), n_(n)
{
    p_.resize(n);
    std::normal_distribution<float> pd(0.0f, static_cast<float>(pos_sigma_m));
    std::normal_distribution<float> hd(0.0f, static_cast<float>(hdg_sigma_deg));
    std::normal_distribution<float> sd(0.0f, 3.0f);
    std::normal_distribution<float> yd(0.0f, static_cast<float>(yaw_sigma_dps));
    float cos_lat = std::cos(static_cast<float>(start_lat) * DEG2RAD);
    float uw = 1.0f / static_cast<float>(n);
    for (auto& pt : p_) {
        pt.lat          = static_cast<float>(start_lat) + pd(rng_) / 111320.0f;
        pt.lon          = static_cast<float>(start_lon) + pd(rng_) / (111320.0f * cos_lat);
        pt.heading_deg  = std::fmod(static_cast<float>(heading_deg) + hd(rng_) + 360.0f, 360.0f);
        pt.speed_mps    = std::max(5.0f, static_cast<float>(speed_mps) + sd(rng_));
        pt.yaw_rate_dps = static_cast<float>(yaw_rate_dps) + yd(rng_);
        pt.weight       = uw;
    }
}

void ContourPF::predict(double dt_s, double hdg_noise, double spd_noise) {
    std::normal_distribution<float> hn(0.0f, static_cast<float>(hdg_noise));
    std::normal_distribution<float> sn(0.0f, static_cast<float>(spd_noise));
    // yaw_rate slowly drifts so it can adapt if the path changes character
    std::normal_distribution<float> yn(0.0f, 0.03f);
    float dt = static_cast<float>(dt_s);
    for (auto& pt : p_) {
        pt.yaw_rate_dps += yn(rng_);
        pt.heading_deg = std::fmod(
            pt.heading_deg + pt.yaw_rate_dps * dt + hn(rng_) + 360.0f, 360.0f);
        pt.speed_mps   = std::max(5.0f, pt.speed_mps + sn(rng_));
        float az      = pt.heading_deg * DEG2RAD;
        float dist    = pt.speed_mps * dt;
        float cos_lat = std::cos(pt.lat * DEG2RAD);
        pt.lat += dist * std::cos(az) / 111320.0f;
        pt.lon += dist * std::sin(az) / (111320.0f * cos_lat);
    }
}

void ContourPF::contour_update(float expected_terrain, float sigma_m) {
    float inv2s2 = 1.0f / (2.0f * sigma_m * sigma_m);
    float total_w = 0.0f;

    for (auto& pt : p_) {
        float h = DemLoader::sample(dem_, pt.lat, pt.lon, -9999.0f);
        if (h < -9000.0f) { pt.weight = 0.0f; continue; }
        float r = expected_terrain - h;
        pt.weight *= std::exp(-r * r * inv2s2);
        total_w += pt.weight;
    }

    if (total_w < 1e-30f) {
        float uw = 1.0f / static_cast<float>(n_);
        for (auto& pt : p_) pt.weight = uw;
        return;
    }
    float inv_w = 1.0f / total_w;
    float neff_inv = 0.0f;
    for (auto& pt : p_) {
        pt.weight *= inv_w;
        neff_inv  += pt.weight * pt.weight;
    }
    if (1.0f / neff_inv < 0.5f * static_cast<float>(n_))
        resample();
}

std::optional<AglEvent> ContourPF::detect_event() const {
    int K = static_cast<int>(agl_buf_.size());
    if (K < AGL_WIN) return std::nullopt;

    int mid = AGL_WIN / 2;
    float v_mid = agl_buf_[K - 1 - mid];
    float v_max = v_mid, v_min = v_mid;
    bool  is_min = true, is_max = true;

    for (int i = 0; i < AGL_WIN; ++i) {
        if (i == mid) continue;
        float v = agl_buf_[K - 1 - i];
        v_max = std::max(v_max, v);
        v_min = std::min(v_min, v);
        if (v <= v_mid) is_min = false;
        if (v >= v_mid) is_max = false;
    }

    constexpr float MIN_DEPTH = 8.0f;
    if (is_min && (v_max - v_mid) > MIN_DEPTH) return AglEvent{AglEvType::PEAK, v_mid};
    if (is_max && (v_mid - v_min) > MIN_DEPTH) return AglEvent{AglEvType::PIT,  v_mid};
    return std::nullopt;
}

// Ridge/valley crossing fix: uses RidgeMap lines (much sharper than peak points).
// Elevation E = baro - AGL pins a specific point along the line.
// Heading perpendicularity confirms the crossing direction.
// Boost radius: 150 m (vs 400 m for peaks) because ridge lines are thin.
void ContourPF::apply_ridge_fix(const AglEvent& ev, double baro_alt_m,
                                  float est_lat, float est_lon, float sigma_m) {
    float kp_elev = static_cast<float>(baro_alt_m) - ev.agl_m;
    FeatureType ft = (ev.type == AglEvType::PEAK) ? FeatureType::RIDGE : FeatureType::VALLEY;

    auto est = estimate();
    auto matches = rm_->query(est_lat, est_lon, 5000.0f,
                               kp_elev, sigma_m * 2.5f,
                               static_cast<float>(est.heading_deg), 65.0f,
                               ft);
    if (matches.empty()) return;

    // Keep only matches above score threshold
    constexpr float MIN_SCORE = 0.55f;
    while (!matches.empty() && matches.back().score < MIN_SCORE)
        matches.pop_back();
    if (matches.empty()) return;

    std::cerr << "[RIDGE-FIX] " << (ft == FeatureType::RIDGE ? "RIDGE" : "VALLEY")
              << " AGL=" << ev.agl_m << " → terrain≈" << kp_elev
              << " m, кандидатов=" << matches.size()
              << " best_score=" << matches[0].score
              << " elev_err=" << matches[0].elev_err << " м\n";

    // Boost particles near any top-N matching ridge points (cover full line segment)
    constexpr int   TOP_N       = 12;
    constexpr float BOOST_R_M   = 150.0f;
    float boost_r_deg = BOOST_R_M / 111320.0f;
    float inv2r2 = 1.0f / (2.0f * boost_r_deg * boost_r_deg);
    int   n_use  = std::min(static_cast<int>(matches.size()), TOP_N);

    float total_w = 0.0f;
    for (auto& pt : p_) {
        float best_d2 = 1e9f;
        for (int mi = 0; mi < n_use; ++mi) {
            float dlat = pt.lat - matches[mi].pt->lat;
            float dlon = pt.lon - matches[mi].pt->lon;
            float d2 = dlat * dlat + dlon * dlon;
            if (d2 < best_d2) best_d2 = d2;
        }
        float boost = 0.02f + 0.98f * std::exp(-best_d2 * inv2r2);
        pt.weight *= boost;
        total_w += pt.weight;
    }

    if (total_w < 1e-30f) return;
    float inv_w = 1.0f / total_w;
    float neff_inv = 0.0f;
    for (auto& pt : p_) {
        pt.weight *= inv_w;
        neff_inv  += pt.weight * pt.weight;
    }
    if (1.0f / neff_inv < 0.25f * static_cast<float>(n_))
        resample();
}

// Fallback: peak/pit point boost via KeypointDB (used when no RidgeMap).
void ContourPF::apply_peak_fix(const AglEvent& ev, double baro_alt_m,
                                 float est_lat, float est_lon, float sigma_m) {
    float kp_elev = static_cast<float>(baro_alt_m) - ev.agl_m;
    KpType kpt = (ev.type == AglEvType::PEAK) ? KpType::PEAK : KpType::PIT;

    auto kps = kpdb_.query(est_lat, est_lon, kp_elev, sigma_m * 3.0f, 6000.0f, kpt);
    if (kps.empty()) return;

    std::cerr << "[KP-FIX] " << (ev.type == AglEvType::PEAK ? "PEAK" : "PIT")
              << " AGL=" << ev.agl_m << " → terrain≈" << kp_elev
              << " m, кандидатов=" << kps.size() << "\n";

    float boost_r_deg = 400.0f / 111320.0f;
    float inv2r2 = 1.0f / (2.0f * boost_r_deg * boost_r_deg);
    float total_w = 0.0f;

    for (auto& pt : p_) {
        float best_d2 = 1e9f;
        for (const auto* kp : kps) {
            float dlat = pt.lat - kp->lat, dlon = pt.lon - kp->lon;
            float d2 = dlat * dlat + dlon * dlon;
            if (d2 < best_d2) best_d2 = d2;
        }
        float boost = 0.05f + 0.95f * std::exp(-best_d2 * inv2r2);
        pt.weight *= boost;
        total_w += pt.weight;
    }

    if (total_w < 1e-30f) return;
    float inv_w = 1.0f / total_w;
    float neff_inv = 0.0f;
    for (auto& pt : p_) {
        pt.weight *= inv_w;
        neff_inv  += pt.weight * pt.weight;
    }
    if (1.0f / neff_inv < 0.3f * static_cast<float>(n_))
        resample();
}

void ContourPF::resample() {
    std::vector<CPFParticle> next;
    next.reserve(n_);
    std::uniform_real_distribution<float> ud(0.0f, 1.0f / static_cast<float>(n_));
    float u  = ud(rng_);
    float c  = p_[0].weight;
    int   i  = 0;
    float uw = 1.0f / static_cast<float>(n_);
    for (int j = 0; j < n_; ++j) {
        while (i < n_ - 1 && u > c) c += p_[++i].weight;
        next.push_back(p_[i]);
        next.back().weight = uw;
        u += uw;
    }
    p_ = std::move(next);
}

// Gradient-based slope constraint.
// For each particle: project DEM gradient at particle position onto its heading
// direction → predicted terrain slope for this step.  Weight particles by
// exp(-(meas_slope - pred_slope)² / 2σ²).
// When terrain is flat, pred_slope ≈ 0 ≈ meas_slope → weights nearly uniform
// → no effect.  On hilly terrain, wrong-heading particles get low weight →
// heading converges to the direction that matches the measured slope.
void ContourPF::slope_update(float meas_slope_m, float dt_s, float sigma_m) {
    float inv2s2 = 1.0f / (2.0f * sigma_m * sigma_m);
    constexpr float EPS = 0.0005f;  // ~55 m, ~2 DEM pixels for gradient
    float total_w = 0.0f;

    for (auto& pt : p_) {
        float cos_lat = std::cos(pt.lat * DEG2RAD);
        float eps_lon = EPS / cos_lat;

        float h_n = DemLoader::sample(dem_, pt.lat + EPS, pt.lon,       -9999.0f);
        float h_s = DemLoader::sample(dem_, pt.lat - EPS, pt.lon,       -9999.0f);
        float h_e = DemLoader::sample(dem_, pt.lat,       pt.lon + eps_lon, -9999.0f);
        float h_w = DemLoader::sample(dem_, pt.lat,       pt.lon - eps_lon, -9999.0f);

        if (h_n < -9000.f || h_s < -9000.f || h_e < -9000.f || h_w < -9000.f) continue;

        // DEM gradient in m/m (north and east components)
        float g_north = (h_n - h_s) / (2.0f * EPS * 111320.0f);
        float g_east  = (h_e - h_w) / (2.0f * eps_lon * 111320.0f * cos_lat);

        // Project onto particle heading: pred slope per step
        float az = pt.heading_deg * DEG2RAD;
        float pred_slope = (g_north * std::cos(az) + g_east * std::sin(az))
                           * pt.speed_mps * static_cast<float>(dt_s);

        float r = meas_slope_m - pred_slope;
        pt.weight *= std::exp(-r * r * inv2s2);
        total_w += pt.weight;
    }

    if (total_w < 1e-30f) return;
    float inv_w = 1.0f / total_w;
    for (auto& pt : p_) pt.weight *= inv_w;
}

ContourPF::Estimate ContourPF::step(double agl_m, double baro_alt_m, double dt_s,
                                     double hdg_noise_deg, double spd_noise_mps,
                                     double meas_sigma_m) {
    float terrain_h = static_cast<float>(baro_alt_m - agl_m);
    float meas_slope = (prev_terrain_h_ > -9000.0f)
                       ? (terrain_h - prev_terrain_h_) : 0.0f;
    prev_terrain_h_ = terrain_h;

    predict(dt_s, hdg_noise_deg, spd_noise_mps);
    contour_update(terrain_h, static_cast<float>(meas_sigma_m));
    // slope_sigma is based on AGL sensor noise (~2m), NOT on position tolerance
    // (meas_sigma_m is position tolerance, unrelated to slope measurement noise).
    // sigma_slope = sqrt(2) * agl_noise ≈ 4m: slope is diff of two AGL readings.
    slope_update(meas_slope, dt_s, 4.0f);

    agl_buf_.push_back(static_cast<float>(agl_m));
    if (static_cast<int>(agl_buf_.size()) > AGL_WIN)
        agl_buf_.pop_front();

    auto ev = detect_event();
    if (ev) {
        auto est = estimate();
        float elat = static_cast<float>(est.lat);
        float elon = static_cast<float>(est.lon);
        float sig  = static_cast<float>(meas_sigma_m);
        if (rm_)
            apply_ridge_fix(*ev, baro_alt_m, elat, elon, sig);
        else
            apply_peak_fix(*ev, baro_alt_m, elat, elon, sig);
    }

    return estimate();
}

ContourPF::Estimate ContourPF::estimate() const {
    double lat = 0, lon = 0, spd = 0, sin_h = 0, cos_h = 0, neff_inv = 0;
    for (const auto& pt : p_) {
        double w = pt.weight;
        lat      += pt.lat       * w;
        lon      += pt.lon       * w;
        spd      += pt.speed_mps * w;
        sin_h    += std::sin(pt.heading_deg * DEG2RAD) * w;
        cos_h    += std::cos(pt.heading_deg * DEG2RAD) * w;
        neff_inv += w * w;
    }
    double hdg = std::atan2(sin_h, cos_h) / DEG2RAD;
    if (hdg < 0.0) hdg += 360.0;
    return {lat, lon, hdg, spd, neff_inv > 1e-30 ? 1.0 / neff_inv : 0.0};
}
