#include "terrain_correlator.hpp"
#include <cmath>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <opencv2/imgproc.hpp>

static constexpr double NCC_CONFIDENCE_THRESHOLD = 0.85;
static constexpr double FLAT_TERRAIN_STD_M       = 5.0;  // профиль с σ < 5 м = слишком плоский

TerrainCorrelator::TerrainCorrelator(const DemData& dem, const CorrelatorConfig& cfg)
    : m_dem(&dem), m_cfg(cfg)
{}

std::optional<CorrelationResult> TerrainCorrelator::add_measurement(double radio_alt_m) {
    double terrain_h = m_cfg.baro_alt_m - radio_alt_m;
    m_measured.push_back(terrain_h);

    if (static_cast<int>(m_measured.size()) >= m_cfg.min_profile_len)
        return run_search();
    return std::nullopt;
}

std::vector<float> TerrainCorrelator::extract_profile(
    double start_lat, double start_lon,
    double azimuth_deg, int n_points,
    double step_px_x, double step_px_y) const
{
    double az_rad = azimuth_deg * M_PI / 180.0;
    double dpx =  std::sin(az_rad) * step_px_x;
    double dpy = -std::cos(az_rad) * step_px_y;

    auto [px0, py0] = DemLoader::geo_to_pixel(m_dem->gt, start_lat, start_lon);

    int nx = m_dem->elev.cols;
    int ny = m_dem->elev.rows;

    std::vector<float> profile;
    profile.reserve(n_points);

    double px = px0, py = py0;
    for (int i = 0; i < n_points; ++i) {
        int xi = static_cast<int>(px);
        int yi = static_cast<int>(py);
        float val = 0.0f;
        if (xi >= 0 && yi >= 0 && xi < nx - 1 && yi < ny - 1) {
            float dx = static_cast<float>(px - xi);
            float dy = static_cast<float>(py - yi);
            val = m_dem->elev.at<float>(yi,   xi  ) * (1-dx)*(1-dy)
                + m_dem->elev.at<float>(yi,   xi+1) * dx*(1-dy)
                + m_dem->elev.at<float>(yi+1, xi  ) * (1-dx)*dy
                + m_dem->elev.at<float>(yi+1, xi+1) * dx*dy;
        }
        profile.push_back(val);
        px += dpx;
        py += dpy;
    }
    return profile;
}

std::vector<float> TerrainCorrelator::ncc_sliding_weighted(
    const std::vector<double>& measured,
    const std::vector<float>&  reference,
    const std::vector<double>& weights)
{
    int n = static_cast<int>(measured.size());
    int m = static_cast<int>(reference.size());
    if (m < n) return {};

    // Взвешенные статистики измеренного профиля (не зависят от смещения d)
    double W = 0;
    for (double w : weights) W += w;

    double mean_m = 0;
    for (int k = 0; k < n; ++k) mean_m += weights[k] * measured[k];
    mean_m /= W;

    double var_m = 0;
    for (int k = 0; k < n; ++k) {
        double d = measured[k] - mean_m;
        var_m += weights[k] * d * d;
    }
    double std_m = std::sqrt(var_m / W);
    if (std_m < 1e-6) return std::vector<float>(m - n + 1, 0.0f);

    int n_offsets = m - n + 1;
    std::vector<float> result(n_offsets);

    for (int d = 0; d < n_offsets; ++d) {
        double mean_r = 0;
        for (int k = 0; k < n; ++k) mean_r += weights[k] * reference[d+k];
        mean_r /= W;

        double var_r = 0;
        for (int k = 0; k < n; ++k) {
            double diff = reference[d+k] - mean_r;
            var_r += weights[k] * diff * diff;
        }
        double std_r = std::sqrt(var_r / W);
        if (std_r < 1e-6) { result[d] = 0.0f; continue; }

        double ncc = 0;
        for (int k = 0; k < n; ++k)
            ncc += weights[k] * (measured[k]-mean_m) * (reference[d+k]-mean_r);
        result[d] = static_cast<float>(ncc / (W * std_m * std_r));
    }
    return result;
}

std::optional<CorrelationResult> TerrainCorrelator::run_search() {
    int n = static_cast<int>(m_measured.size());
    if (n < m_cfg.min_profile_len) return std::nullopt;

    double mpp_x  = DemLoader::meters_per_pixel_x(m_dem->gt, m_cfg.start_lat);
    double mpp_y  = DemLoader::meters_per_pixel_y(m_dem->gt);
    double step_m = m_cfg.speed_mps * m_cfg.sample_dt_s;
    double step_px_x = step_m / mpp_x;
    double step_px_y = step_m / mpp_y;

    int search_px = static_cast<int>(m_cfg.search_radius_m / step_m);
    int ref_len   = n + search_px;
    int n_azimuths = 360 / m_cfg.azimuth_step_deg;
    int n_offsets  = search_px + 1;

    // Вычислить σ измеренного профиля (информативность рельефа)
    double mean_m = 0;
    for (double v : m_measured) mean_m += v;
    mean_m /= n;
    double std_m = 0;
    for (double v : m_measured) std_m += (v - mean_m) * (v - mean_m);
    std_m = std::sqrt(std_m / n);

    // Веса: w[k] = |M[k] - mean_M| + epsilon
    // Точки далеко от среднего (пики, долины) получают больший вес
    std::vector<double> weights(n);
    for (int k = 0; k < n; ++k)
        weights[k] = std::abs(m_measured[k] - mean_m) + 1.0;

    cv::Mat corr_map(n_azimuths, n_offsets, CV_32F, cv::Scalar(0));
    double best_ncc = -2.0;
    int best_az = 0, best_off = 0;

    for (int ai = 0; ai < n_azimuths; ++ai) {
        double az = ai * m_cfg.azimuth_step_deg;
        auto ref = extract_profile(m_cfg.start_lat, m_cfg.start_lon,
                                   az, ref_len, step_px_x, step_px_y);
        if (static_cast<int>(ref.size()) < n) continue;

        auto ncc_vec = ncc_sliding_weighted(m_measured, ref, weights);

        for (int d = 0; d < static_cast<int>(ncc_vec.size()) && d < n_offsets; ++d) {
            corr_map.at<float>(ai, d) = ncc_vec[d];
            if (ncc_vec[d] > best_ncc) {
                best_ncc = ncc_vec[d];
                best_az  = ai * m_cfg.azimuth_step_deg;
                best_off = d;
            }
        }
    }

    double best_offset_m = best_off * step_m;
    double az_rad  = best_az * M_PI / 180.0;
    double cos_lat = std::cos(m_cfg.start_lat * M_PI / 180.0);

    double dlat_start = (best_offset_m * std::cos(az_rad)) / 111320.0;
    double dlon_start = (best_offset_m * std::sin(az_rad)) / (111320.0 * cos_lat);

    double travel_m  = (n - 1) * step_m;
    double dlat_cur  = dlat_start + (travel_m * std::cos(az_rad)) / 111320.0;
    double dlon_cur  = dlon_start + (travel_m * std::sin(az_rad)) / (111320.0 * cos_lat);

    bool low_conf = (best_ncc < NCC_CONFIDENCE_THRESHOLD) || (std_m < FLAT_TERRAIN_STD_M);

    CorrelationResult res;
    res.best_azimuth_deg = best_az;
    res.best_offset_m    = best_offset_m;
    res.best_ncc         = best_ncc;
    res.start_lat        = m_cfg.start_lat;
    res.start_lon        = m_cfg.start_lon;
    res.found_lat        = m_cfg.start_lat + dlat_start;
    res.found_lon        = m_cfg.start_lon + dlon_start;
    res.current_lat      = m_cfg.start_lat + dlat_cur;
    res.current_lon      = m_cfg.start_lon + dlon_cur;
    res.ground_speed_mps = m_cfg.speed_mps;
    res.profile_len      = n;
    res.corr_map         = corr_map;
    res.low_confidence   = low_conf;
    res.profile_std      = std_m;

    std::cerr << "[CORR] az=" << best_az << "° NCC=" << best_ncc
              << " σ_profile=" << std_m << " м"
              << (low_conf ? " [LOW CONFIDENCE]" : "")
              << " cur=(" << res.current_lat << ", " << res.current_lon << ")\n";

    return res;
}
