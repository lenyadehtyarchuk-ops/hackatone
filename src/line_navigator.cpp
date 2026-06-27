#include "line_navigator.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

static constexpr double DEG2RAD = M_PI / 180.0;

LineNavigator::LineNavigator(const RidgeMap& rm,
                             double start_lat, double start_lon,
                             double heading_deg, double speed_mps,
                             const LNavConfig& cfg)
    : rm_(rm), cfg_(cfg),
      lat_(start_lat), lon_(start_lon),
      hdg_(heading_deg), spd_(speed_mps),
      hdg_noise_acc_(0.0),
      total_fixes_(0),
      time_s_(0.0)
{}

void LineNavigator::dead_reckon(double dt_s, double hdg_noise, double spd_noise) {
    // Update speed (noisy, bounded)
    spd_ = std::max(5.0, spd_ + spd_noise * 0.0);  // hold speed unless crossing updates
    hdg_ = std::fmod(hdg_ + 360.0, 360.0);

    // Accumulate heading uncertainty sigma (random walk)
    hdg_noise_acc_ = std::sqrt(hdg_noise_acc_*hdg_noise_acc_ + hdg_noise*hdg_noise);

    double az = hdg_ * DEG2RAD;
    double dist = spd_ * dt_s;
    double cos_lat = std::cos(lat_ * DEG2RAD);
    lat_ += dist * std::cos(az) / 111320.0;
    lon_ += dist * std::sin(az) / (111320.0 * cos_lat);
}

std::optional<LineNavigator::AglCrossing> LineNavigator::detect_crossing() const {
    int K = static_cast<int>(agl_buf_.size());
    if (K < cfg_.agl_win) return std::nullopt;

    int mid = cfg_.agl_win / 2;
    float v_mid = agl_buf_[K - 1 - mid];
    float v_max = v_mid, v_min = v_mid;
    bool is_min = true, is_max = true;

    for (int i = 0; i < cfg_.agl_win; ++i) {
        if (i == mid) continue;
        float v = agl_buf_[K - 1 - i];
        v_max = std::max(v_max, v);
        v_min = std::min(v_min, v);
        if (v <= v_mid) is_min = false;
        if (v >= v_mid) is_max = false;
    }

    if (is_min && (v_max - v_mid) > cfg_.agl_event_depth)
        return AglCrossing{FeatureType::RIDGE, v_mid};  // AGL minimum = terrain peak = ridge crossing
    if (is_max && (v_mid - v_min) > cfg_.agl_event_depth)
        return AglCrossing{FeatureType::VALLEY, v_mid};  // AGL maximum = terrain dip = valley crossing
    return std::nullopt;
}

bool LineNavigator::apply_fix(double fix_lat, double fix_lon) {
    // Sanity check: implied speed from previous fix
    if (!fix_history_.empty()) {
        const auto& prev = fix_history_.back();
        double dt = time_s_ - prev.timestamp_s;
        if (dt > 0.1) {
            double dlat = fix_lat - prev.lat;
            double dlon = (fix_lon - prev.lon) * std::cos(prev.lat * DEG2RAD);
            double dist_m = std::sqrt(dlat*dlat + dlon*dlon) * 111320.0;
            double implied_spd = dist_m / dt;
            if (implied_spd > cfg_.max_speed_mps) {
                std::cerr << "[LNAV] Отклонён fix: скорость " << implied_spd << " м/с\n";
                return false;
            }
        }
    }

    // Update heading from consecutive fixes
    if (!fix_history_.empty()) {
        const auto& prev = fix_history_.back();
        double dlat = fix_lat - prev.lat;
        double dlon = (fix_lon - prev.lon) * std::cos(fix_lat * DEG2RAD);
        if (std::abs(dlat) + std::abs(dlon) > 1e-9) {
            double new_hdg = std::atan2(dlon, dlat) / DEG2RAD;
            // Blend: trust new heading estimate, but not 100%
            double diff = std::fmod(new_hdg - hdg_ + 540.0, 360.0) - 180.0;
            hdg_ = std::fmod(hdg_ + 0.7 * diff + 360.0, 360.0);
        }
    }

    lat_ = fix_lat;
    lon_ = fix_lon;
    hdg_noise_acc_ = 0.0;  // reset accumulated heading uncertainty
    ++total_fixes_;

    fix_history_.push_back({fix_lat, fix_lon, time_s_});
    if (static_cast<int>(fix_history_.size()) > 20) fix_history_.pop_front();
    return true;
}

LineNavigator::Estimate LineNavigator::step(double agl_m, double baro_alt_m,
                                             double dt_s,
                                             double hdg_noise_deg,
                                             double spd_noise_mps) {
    time_s_ += dt_s;
    dead_reckon(dt_s, hdg_noise_deg, spd_noise_mps);

    agl_buf_.push_back(static_cast<float>(agl_m));
    if (static_cast<int>(agl_buf_.size()) > cfg_.agl_win)
        agl_buf_.pop_front();

    bool fix_applied = false;
    auto crossing = detect_crossing();
    if (crossing) {
        float target_elev = static_cast<float>(baro_alt_m) - crossing->agl_m;

        // The search is centered on our CURRENT dead-reckoning estimate.
        // The crossing happened ~mid steps ago, but since radius >> position error,
        // the search will still find the correct ridge point.
        auto matches = rm_.query(
                static_cast<float>(lat_),
                static_cast<float>(lon_),
                cfg_.search_radius_m,
                target_elev, cfg_.elev_tol_m,
                static_cast<float>(hdg_),
                cfg_.hdg_tol_deg,
                crossing->type);

        if (!matches.empty()) {
            const auto& best = matches[0];
            std::cerr << "[LNAV-FIX] "
                      << (crossing->type == FeatureType::RIDGE ? "RIDGE" : "VALLEY")
                      << " AGL=" << crossing->agl_m
                      << " terrain≈" << target_elev
                      << " m → fix(" << best.pt->lat << "," << best.pt->lon
                      << ") elev_err=" << best.elev_err
                      << " hdg_err=" << best.hdg_err_deg
                      << "° score=" << best.score << "\n";

            fix_applied = apply_fix(best.pt->lat, best.pt->lon);
        } else {
            std::cerr << "[LNAV] "
                      << (crossing->type == FeatureType::RIDGE ? "RIDGE" : "VALLEY")
                      << " событие: terrain≈" << target_elev
                      << " м — нет совпадения в радиусе " << cfg_.search_radius_m << " м\n";
        }
    }

    return {lat_, lon_, hdg_, spd_, fix_applied, total_fixes_};
}

LineNavigator::Estimate LineNavigator::current() const {
    return {lat_, lon_, hdg_, spd_, false, total_fixes_};
}
