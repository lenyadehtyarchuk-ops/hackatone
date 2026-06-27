#pragma once
#include "dem_loader.hpp"
#include "keypoint_db.hpp"
#include "ridge_map.hpp"
#include <vector>
#include <deque>
#include <random>
#include <optional>

struct CPFParticle {
    float lat, lon, heading_deg, speed_mps, weight;
    float yaw_rate_dps;  // °/s, estimated from GPS history at jammer entry
};

enum class AglEvType { PEAK, PIT };

struct AglEvent {
    AglEvType type;
    float agl_m;
};

// Contour Particle Filter — combined mode.
//
// Every step:
//   1. predict       — move particles (heading + speed + noise)
//   2. contour_update — weight *= exp(-residual²/2σ²), 1 DEM lookup/particle
//   3. ridge_fix     — if AGL extremum detected AND RidgeMap is set:
//                       query ridge lines → boost particles near matching segment
//                       (replaces old peak-only KeypointDB boost)
//
// The RidgeMap fix is much sharper than peak-only because:
//   - elevation E pins a specific POINT on the ridge line
//   - heading perpendicularity confirms the crossing direction
//   → 150 m boost radius vs 400 m for peaks
class ContourPF {
public:
    // rm: optional RidgeMap for ridge/valley crossing fixes (nullptr = peaks only via kpdb)
    // yaw_rate_dps: initial turn rate from GPS history (°/s); spread by yaw_sigma
    ContourPF(const DemData& dem, const KeypointDB& kpdb,
              double start_lat, double start_lon,
              double heading_deg, double speed_mps,
              int    n              = 3000,
              double pos_sigma_m    = 100.0,
              double hdg_sigma_deg  = 30.0,
              const RidgeMap* rm    = nullptr,
              double yaw_rate_dps   = 0.0,
              double yaw_sigma_dps  = 0.2);

    struct Estimate {
        double lat, lon, heading_deg, speed_mps, neff;
    };

    Estimate step(double agl_m, double baro_alt_m, double dt_s,
                  double hdg_noise_deg = 3.0,
                  double spd_noise_mps = 1.0,
                  double meas_sigma_m  = 15.0);

    Estimate estimate() const;
    int n() const { return n_; }

private:
    std::vector<CPFParticle> p_;
    const DemData&    dem_;
    const KeypointDB& kpdb_;
    const RidgeMap*   rm_;
    std::mt19937      rng_;
    int               n_;

    std::deque<float> agl_buf_;
    static constexpr int AGL_WIN = 7;

    float prev_terrain_h_ = -9999.0f;

    void predict(double dt_s, double hdg_noise, double spd_noise);
    void contour_update(float expected_terrain, float sigma_m);
    // Gradient-based slope constraint: weights particles by how well their
    // heading direction matches the measured terrain slope (dh/dt).
    // Analogous to TERCOM's gradient NCC channel but applied per-particle.
    void slope_update(float meas_slope_m, float dt_s, float sigma_m);
    void apply_ridge_fix(const AglEvent& ev, double baro_alt_m,
                         float est_lat, float est_lon, float sigma_m);
    void apply_peak_fix(const AglEvent& ev, double baro_alt_m,
                        float est_lat, float est_lon, float sigma_m);
    std::optional<AglEvent> detect_event() const;
    void resample();
};
