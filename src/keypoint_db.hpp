#pragma once
#include "dem_loader.hpp"
#include <vector>
#include <cstdint>

enum class KpType : uint8_t { PEAK = 0, PIT = 1 };

struct Keypoint {
    float lat, lon, elev;
    KpType type;
};

// Precomputed database of terrain extrema (peaks and pits) extracted from DEM.
// Used by ContourPF for keypoint-event-based weight boosting.
class KeypointDB {
public:
    // window: half-side of the local extremum search box (pixels)
    // min_prominence: minimum elevation range within window to qualify
    KeypointDB(const DemData& dem, int window = 15, float min_prominence = 30.0f);

    // Return pointers to keypoints of given type near (lat,lon) within radius_m,
    // whose elevation is within elev_tol of target_elev.
    std::vector<const Keypoint*> query(float lat, float lon,
                                        float target_elev, float elev_tol,
                                        float radius_m, KpType type) const;

    std::size_t size() const { return kps_.size(); }
    const std::vector<Keypoint>& all() const { return kps_; }

private:
    std::vector<Keypoint> kps_;

    // Spatial grid for O(1) radius pre-filter
    float lat0_, lon0_;
    float cell_deg_;
    int   grows_, gcols_;
    std::vector<std::vector<uint32_t>> grid_;

    void build_grid(float cell_deg);
};
