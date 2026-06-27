#pragma once
#include "ridge_map.hpp"
#include <deque>
#include <vector>
#include <optional>

// LineNavigator: GPS-denied navigation using ridge/valley line crossings.
//
// Algorithm:
//  1. Between crossings: dead-reckoning from last known position (heading + speed).
//  2. When a crossing is detected (local AGL extremum):
//     a. Compute expected terrain elevation E = baro_alt - AGL_min/max.
//     b. Query RidgeMap: find ridge/valley points with elevation ≈ E,
//        within search radius, whose direction is ≈ perpendicular to UAV heading.
//     c. Best match → new position fix (equivalent to GPS fix, 30–60 m accuracy).
//     d. Correct dead-reckoning estimate.
//  3. Sequence consistency: maintain last N fixes; if a new fix is spatially
//     inconsistent (implies impossible speed), reject it.
//
// The crossing gives a 2D position because elevation E pins a specific POINT
// on the ridge line (elevation varies along the line), and the perpendicular
// heading further confirms it.
struct LNavConfig {
    float search_radius_m  = 4000.0f;
    float elev_tol_m       = 30.0f;
    float hdg_tol_deg      = 65.0f;
    float agl_event_depth  = 8.0f;
    int   agl_win          = 7;
    float max_speed_mps    = 120.0f;
};

class LineNavigator {
public:
    LineNavigator(const RidgeMap& rm,
                  double start_lat, double start_lon,
                  double heading_deg, double speed_mps,
                  const LNavConfig& cfg = LNavConfig{});

    struct Estimate {
        double lat, lon, heading_deg, speed_mps;
        bool   fix_applied;
        int    total_fixes;
    };

    Estimate step(double agl_m, double baro_alt_m, double dt_s,
                  double hdg_noise_deg = 2.0,
                  double spd_noise_mps = 1.0);

    Estimate current() const;

private:
    const RidgeMap& rm_;
    LNavConfig cfg_;

    // Dead-reckoning state
    double lat_, lon_, hdg_, spd_;
    double hdg_noise_acc_;  // accumulated heading noise sigma
    int    total_fixes_;

    // AGL ring buffer for crossing detection
    std::deque<float> agl_buf_;

    // Recent fixes for consistency check
    struct Fix {
        double lat, lon;
        double timestamp_s;
    };
    std::deque<Fix> fix_history_;
    double time_s_;

    struct AglCrossing {
        FeatureType type;
        float agl_m;
    };
    std::optional<AglCrossing> detect_crossing() const;
    bool apply_fix(double fix_lat, double fix_lon);
    void dead_reckon(double dt_s, double hdg_noise, double spd_noise);
};
