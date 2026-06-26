#include "terrain_correlator.hpp"
#include <cmath>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <opencv2/imgproc.hpp>

static constexpr double NCC_CONFIDENCE_THRESHOLD = 0.80;  // для σ_датчика=2м
static constexpr double FLAT_TERRAIN_STD_M       = 5.0;

TerrainCorrelator::TerrainCorrelator(const DemData& dem, const CorrelatorConfig& cfg)
    : m_dem(&dem), m_cfg(cfg)
{}

std::optional<CorrelationResult> TerrainCorrelator::add_measurement(double radio_alt_m) {
    double terrain_h = m_cfg.baro_alt_m - radio_alt_m;
    m_measured.push_back(terrain_h);

    // Накапливаем градиентный профиль: разность с предыдущей точкой
    if (m_measured.size() >= 2) {
        m_grad.push_back(m_measured.back() - m_measured[m_measured.size() - 2]);
    }

    if (static_cast<int>(m_measured.size()) >= m_cfg.min_profile_len)
        return run_search();
    return std::nullopt;
}

// --- Извлечение профиля из ЦМР ---
struct ExtractedProfile {
    std::vector<float> elev;  // сырые высоты
    std::vector<float> grad;  // градиент: elev[k] - elev[k-1]
};

static ExtractedProfile extract_profile_full(
    const DemData& dem,
    double start_lat, double start_lon,
    double azimuth_deg, int n_points,
    double step_px_x, double step_px_y)
{
    double az_rad = azimuth_deg * M_PI / 180.0;
    double dpx =  std::sin(az_rad) * step_px_x;
    double dpy = -std::cos(az_rad) * step_px_y;

    auto [px0, py0] = DemLoader::geo_to_pixel(dem.gt, start_lat, start_lon);

    int nx = dem.elev.cols;
    int ny = dem.elev.rows;

    ExtractedProfile result;
    result.elev.reserve(n_points);

    double px = px0, py = py0;
    float prev = 0.0f;
    for (int i = 0; i < n_points; ++i) {
        int xi = static_cast<int>(px);
        int yi = static_cast<int>(py);
        float val = 0.0f;
        if (xi >= 0 && yi >= 0 && xi < nx - 1 && yi < ny - 1) {
            float dx = static_cast<float>(px - xi);
            float dy = static_cast<float>(py - yi);
            val = dem.elev.at<float>(yi,   xi  ) * (1-dx)*(1-dy)
                + dem.elev.at<float>(yi,   xi+1) * dx*(1-dy)
                + dem.elev.at<float>(yi+1, xi  ) * (1-dx)*dy
                + dem.elev.at<float>(yi+1, xi+1) * dx*dy;
        }
        result.elev.push_back(val);
        if (i > 0) result.grad.push_back(val - prev);
        prev = val;
        px += dpx;
        py += dpy;
    }
    return result;
}

// --- Взвешенная нормированная кросс-корреляция ---
// measured[0..n-1], reference[0..m-1], m >= n
// Возвращает вектор NCC для каждого смещения d=0..m-n
static std::vector<float> ncc_sliding_weighted(
    const std::vector<double>& measured,
    const std::vector<float>&  reference,
    const std::vector<double>& weights)
{
    int n = static_cast<int>(measured.size());
    int m = static_cast<int>(reference.size());
    if (m < n) return {};

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

// Для градиентного канала measured — m_grad (double), reference — ref.grad (float)
// Градиент на 1 точку короче профиля, поэтому передаём n-1 элементов
static std::vector<float> ncc_grad_channel(
    const std::vector<double>& meas_grad,  // длина n-1
    const std::vector<float>&  ref_grad,   // длина m-1
    const std::vector<double>& weights_g)  // длина n-1
{
    return ncc_sliding_weighted(meas_grad, ref_grad, weights_g);
}

std::optional<CorrelationResult> TerrainCorrelator::run_search() {
    int n = static_cast<int>(m_measured.size());
    if (n < m_cfg.min_profile_len) return std::nullopt;

    double mpp_x  = DemLoader::meters_per_pixel_x(m_dem->gt, m_cfg.start_lat);
    double mpp_y  = DemLoader::meters_per_pixel_y(m_dem->gt);
    double step_m = m_cfg.speed_mps * m_cfg.sample_dt_s;
    double step_px_x = step_m / mpp_x;
    double step_px_y = step_m / mpp_y;

    int search_px  = static_cast<int>(m_cfg.search_radius_m / step_m);
    int ref_len    = n + search_px;
    int n_offsets  = search_px + 1;
    int n_azimuths = 360 / m_cfg.azimuth_step_deg;

    // σ измеренного профиля (для low_confidence)
    double mean_m = 0;
    for (double v : m_measured) mean_m += v;
    mean_m /= n;
    double std_m = 0;
    for (double v : m_measured) std_m += (v - mean_m) * (v - mean_m);
    std_m = std::sqrt(std_m / n);

    // Веса для высотного канала
    std::vector<double> weights_h(n);
    for (int k = 0; k < n; ++k)
        weights_h[k] = std::abs(m_measured[k] - mean_m) + 1.0;

    // Веса для градиентного канала (по самому градиенту)
    int ng = static_cast<int>(m_grad.size());   // ng = n-1
    std::vector<double> weights_g;
    if (ng > 0) {
        double mean_g = 0;
        for (double v : m_grad) mean_g += v;
        mean_g /= ng;
        weights_g.resize(ng);
        for (int k = 0; k < ng; ++k)
            weights_g[k] = std::abs(m_grad[k] - mean_g) + 1.0;
    }

    double alpha_h = 1.0 - m_cfg.gradient_weight;  // вес высотного канала
    double alpha_g = m_cfg.gradient_weight;          // вес градиентного канала

    // Секторный поиск: если prior_heading задан, ограничиваем азимут
    bool use_sector = (m_cfg.prior_heading_deg >= 0.0);

    cv::Mat corr_map(n_azimuths, n_offsets, CV_32F, cv::Scalar(0));
    double best_ncc = -2.0;
    int best_az = 0, best_off = 0;

    for (int ai = 0; ai < n_azimuths; ++ai) {
        double az = ai * m_cfg.azimuth_step_deg;

        // Пропустить, если за пределами сектора (с учётом цикличности 360°)
        if (use_sector) {
            double az_norm = az;
            // Привести az к диапазону [prior_heading - 180, prior_heading + 180]
            double center = m_cfg.prior_heading_deg;
            while (az_norm - center > 180.0)  az_norm -= 360.0;
            while (az_norm - center < -180.0) az_norm += 360.0;
            if (std::abs(az_norm - center) > m_cfg.heading_search_deg) continue;
        }

        auto ref = extract_profile_full(*m_dem,
                                        m_cfg.start_lat, m_cfg.start_lon,
                                        az, ref_len, step_px_x, step_px_y);
        if (static_cast<int>(ref.elev.size()) < n) continue;

        // Высотный канал
        auto ncc_h = ncc_sliding_weighted(m_measured, ref.elev, weights_h);

        // Градиентный канал (если включён и есть данные)
        std::vector<float> ncc_g;
        if (alpha_g > 1e-6 && ng > 0 &&
            static_cast<int>(ref.grad.size()) >= ng) {
            ncc_g = ncc_grad_channel(m_grad, ref.grad, weights_g);
        }

        for (int d = 0; d < static_cast<int>(ncc_h.size()) && d < n_offsets; ++d) {
            float combined = static_cast<float>(alpha_h) * ncc_h[d];
            if (!ncc_g.empty() && d < static_cast<int>(ncc_g.size()))
                combined += static_cast<float>(alpha_g) * ncc_g[d];

            corr_map.at<float>(ai, d) = combined;
            if (combined > best_ncc) {
                best_ncc = combined;
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
