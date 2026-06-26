#include "particle_filter.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

static constexpr double DEG2RAD = M_PI / 180.0;

static void move_latlon(double lat, double lon,
                        double heading_deg, double dist_m,
                        double& out_lat, double& out_lon) {
    double az  = heading_deg * DEG2RAD;
    out_lat = lat + dist_m * std::cos(az) / 111320.0;
    out_lon = lon + dist_m * std::sin(az) /
                    (111320.0 * std::cos(lat * DEG2RAD));
}

ParticleFilter::ParticleFilter(const DemData& dem,
                               double start_lat,
                               double start_lon,
                               double heading_deg,
                               double speed_mps,
                               int    n_particles,
                               int    profile_win,
                               double pos_sigma_m,
                               double hdg_sigma_deg)
    : dem_(dem), rng_(42), n_(n_particles), win_(profile_win), speed_nom_(speed_mps)
{
    p_.resize(n_particles);

    std::normal_distribution<double> pos_d(0.0, pos_sigma_m);
    std::normal_distribution<double> hdg_d(heading_deg, hdg_sigma_deg);
    std::normal_distribution<double> spd_d(speed_mps, 3.0);

    double cos_lat = std::cos(start_lat * DEG2RAD);
    double uw      = 1.0 / n_particles;

    for (auto& pt : p_) {
        pt.lat         = start_lat + pos_d(rng_) / 111320.0;
        pt.lon         = start_lon + pos_d(rng_) / (111320.0 * cos_lat);
        pt.heading_deg = std::fmod(hdg_d(rng_) + 360.0, 360.0);
        pt.speed_mps   = std::max(5.0, spd_d(rng_));
        pt.weight      = uw;
    }
}

void ParticleFilter::predict(double dt_s,
                             double hdg_noise,
                             double spd_noise) {
    std::normal_distribution<double> hn(0.0, hdg_noise);
    std::normal_distribution<double> sn(0.0, spd_noise);
    for (auto& pt : p_) {
        pt.heading_deg = std::fmod(pt.heading_deg + hn(rng_) + 360.0, 360.0);
        pt.speed_mps   = std::max(5.0, pt.speed_mps + sn(rng_));
        move_latlon(pt.lat, pt.lon, pt.heading_deg, pt.speed_mps * dt_s,
                    pt.lat, pt.lon);
        pt.hdg_history.push_front(pt.heading_deg);
        if ((int)pt.hdg_history.size() > win_)
            pt.hdg_history.pop_back();
    }
}

// NCC между двумя профилями длины K (инвариантна к постоянному смещению)
static double ncc_profiles(const double* meas, const double* ref, int K) {
    double mm = 0.0, mr = 0.0;
    for (int k = 0; k < K; ++k) { mm += meas[k]; mr += ref[k]; }
    mm /= K; mr /= K;
    double num = 0.0, sm = 0.0, sr = 0.0;
    for (int k = 0; k < K; ++k) {
        double dm = meas[k] - mm, dr = ref[k] - mr;
        num += dm * dr;
        sm  += dm * dm;
        sr  += dr * dr;
    }
    double denom = std::sqrt(sm * sr);
    return denom < 1e-6 ? 0.0 : num / denom;
}

void ParticleFilter::profile_update(double baro_alt_m,
                                    double dt_s,
                                    double /*meas_sigma_m*/) {
    int K = static_cast<int>(agl_buf_.size());

    // σ измеренного профиля — плоский рельеф → слабая дискриминация
    double mm = 0.0;
    for (int k = 0; k < K; ++k) mm += agl_buf_[k];
    mm /= K;
    double sm = 0.0;
    for (int k = 0; k < K; ++k) sm += (agl_buf_[k]-mm)*(agl_buf_[k]-mm);
    double profile_std = std::sqrt(sm / K);
    if (profile_std < 5.0) return;  // рельеф слишком плоский

    // Измеренный профиль: agl_buf_[K-1] = текущая позиция (k=0)
    std::vector<double> meas(K);
    for (int k = 0; k < K; ++k) meas[k] = agl_buf_[K - 1 - k];

    std::vector<double> ref(K);
    double total_w = 0.0;
    constexpr double LAMBDA = 5.0;

    for (auto& pt : p_) {
        double step_m = pt.speed_mps * dt_s;
        bool valid = true;

        // k=0: текущая позиция
        double plat = pt.lat, plon = pt.lon;
        float h = DemLoader::sample(dem_, plat, plon, -9999.0f);
        if (h < -9000.0f) { pt.weight = 0.0; continue; }
        ref[0] = baro_alt_m - static_cast<double>(h);

        // k=1..K-1: идём назад по реальным курсам частицы (не прямая!)
        for (int k = 1; k < K; ++k) {
            // курс которым летели ДО текущей позиции на k шагов назад
            double hdg_k = (k - 1 < (int)pt.hdg_history.size())
                           ? pt.hdg_history[k - 1]
                           : pt.heading_deg;  // fallback если истории мало
            double az_b = std::fmod(hdg_k + 180.0, 360.0) * DEG2RAD;
            double cos_lat = std::cos(plat * DEG2RAD);
            plat += step_m * std::cos(az_b) / 111320.0;
            plon += step_m * std::sin(az_b) / (111320.0 * cos_lat);
            h = DemLoader::sample(dem_, plat, plon, -9999.0f);
            if (h < -9000.0f) { valid = false; break; }
            ref[k] = baro_alt_m - static_cast<double>(h);
        }

        if (!valid) { pt.weight = 0.0; continue; }

        double ncc = ncc_profiles(meas.data(), ref.data(), K);
        // NCC ∈ [-1,1]: вес = exp(λ·ncc); плохие профили (ncc<0) → малый вес
        pt.weight *= std::exp(LAMBDA * ncc);
        total_w   += pt.weight;
    }

    if (total_w < 1e-300) {
        double uw = 1.0 / n_;
        for (auto& pt : p_) pt.weight = uw;
        return;  // мёртвое счисление без лога (это норма при flat terrain)
    }

    double inv_w = 1.0 / total_w;
    double neff  = 0.0;
    for (auto& pt : p_) {
        pt.weight *= inv_w;
        neff += pt.weight * pt.weight;
    }
    if (1.0 / neff < 0.5 * n_)
        resample();
}

void ParticleFilter::resample() {
    std::vector<Particle> next;
    next.reserve(n_);

    std::uniform_real_distribution<double> ud(0.0, 1.0 / n_);
    double u  = ud(rng_);
    double c  = p_[0].weight;
    int    i  = 0;
    double uw = 1.0 / n_;

    for (int j = 0; j < n_; ++j) {
        while (i < n_ - 1 && u > c)
            c += p_[++i].weight;
        next.push_back(p_[i]);
        next.back().weight = uw;
        u += uw;
    }
    p_ = std::move(next);
}

PFEstimate ParticleFilter::step(double agl_m,
                                double baro_alt_m,
                                double dt_s,
                                double hdg_noise_deg,
                                double spd_noise_mps,
                                double meas_sigma_m) {
    predict(dt_s, hdg_noise_deg, spd_noise_mps);

    agl_buf_.push_back(agl_m);
    if (static_cast<int>(agl_buf_.size()) > win_)
        agl_buf_.pop_front();

    // Обновлять только когда накоплен полный профиль
    if (static_cast<int>(agl_buf_.size()) == win_)
        profile_update(baro_alt_m, dt_s, meas_sigma_m);

    return estimate();
}

PFEstimate ParticleFilter::estimate() const {
    double lat = 0.0, lon = 0.0, spd = 0.0;
    double sin_h = 0.0, cos_h = 0.0, neff_inv = 0.0;

    for (const auto& pt : p_) {
        lat    += pt.lat        * pt.weight;
        lon    += pt.lon        * pt.weight;
        spd    += pt.speed_mps  * pt.weight;
        sin_h  += std::sin(pt.heading_deg * DEG2RAD) * pt.weight;
        cos_h  += std::cos(pt.heading_deg * DEG2RAD) * pt.weight;
        neff_inv += pt.weight * pt.weight;
    }

    double hdg = std::atan2(sin_h, cos_h) / DEG2RAD;
    if (hdg < 0.0) hdg += 360.0;

    return {lat, lon, hdg, spd, 1.0 / neff_inv};
}
