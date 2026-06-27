#pragma once
#include "dem_loader.hpp"
#include <vector>
#include <cstdint>
#include <opencv2/core.hpp>

enum class FeatureType : uint8_t { NONE = 0, RIDGE = 1, VALLEY = 2 };

// One point on a ridge or valley centerline, after Hessian + NMS thinning.
struct RidgePoint {
    float lat, lon, elev;
    FeatureType type;
    float dir_deg;    // geographic direction of the line, [0, 180)
                      // (bidirectional: ridge goes both ways)
    float curvature;  // |λ_min| for ridges, λ_max for valleys
};

class RidgeMap {
public:
    // gauss_sigma_px : Gaussian blur sigma before Hessian (smooths DEM noise)
    // curv_thresh_rel: keep points with curvature > curv_thresh_rel * max_curvature
    // ratio_thresh   : |λ1/λ2| must be < this to qualify (asymmetric curvature)
    RidgeMap(const DemData& dem,
             float gauss_sigma_px  = 3.0f,
             float curv_thresh_rel = 0.03f,
             float ratio_thresh    = 0.4f);

    // Find ridge/valley crossing candidates.
    // lat, lon    : current UAV estimated position
    // radius_m    : spatial search radius
    // target_elev : expected terrain elevation at crossing (baro_alt - AGL)
    // elev_tol    : ±tolerance on elevation match
    // uav_hdg_deg : UAV heading (0-360°)
    // hdg_tol_deg : max deviation from perpendicular crossing (0=exactly perp, 90=any angle)
    // type        : RIDGE or VALLEY
    struct Match {
        const RidgePoint* pt;
        float dist_m;
        float elev_err;
        float hdg_err_deg;  // deviation from ideal perpendicular; 0 = perfectly perp
        float score;        // combined quality (higher = better)
    };

    std::vector<Match> query(float lat, float lon, float radius_m,
                              float target_elev, float elev_tol,
                              float uav_hdg_deg, float hdg_tol_deg,
                              FeatureType type) const;

    std::size_t size() const { return pts_.size(); }
    const std::vector<RidgePoint>& points() const { return pts_; }

    // Overlay ridge/valley pixels on a BGR image (for visualization)
    cv::Mat make_overlay(const DemData& dem, int scale = 1) const;

private:
    std::vector<RidgePoint> pts_;

    // Spatial grid: divide DEM extent into cells, each stores point indices
    float lat0_, lon0_;
    float cell_deg_;
    int   grows_, gcols_;
    std::vector<std::vector<uint32_t>> grid_;

    void build_grid(float cell_deg);
};
