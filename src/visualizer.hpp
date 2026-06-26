#pragma once
#include "terrain_correlator.hpp"
#include "dem_loader.hpp"
#include <string>
#include <vector>
#include <opencv2/core.hpp>

struct TrajectoryPoint {
    double lat, lon;
    double speed_mps;
    double heading_deg;
    double ncc;
    bool   gps_denied = false;  // true = точка внутри зоны глушения
};

class Visualizer {
public:
    // Сохранить тепловую карту корреляции (360 азимутов × расстояние)
    static void save_correlation_heatmap(const cv::Mat& corr_map,
                                         const std::string& path,
                                         int best_az = -1,
                                         int best_off = -1);

    // Наложить найденную траекторию на карту высот.
    // jammer_zone: {lat1, lon1, lat2, lon2} или пустой вектор
    static void save_trajectory_on_dem(const DemData& dem,
                                        const std::vector<TrajectoryPoint>& traj,
                                        double start_lat, double start_lon,
                                        const std::string& path,
                                        const std::vector<double>& jammer_zone = {});

    // Показать в окне (если есть display)
    static void show(const std::string& title, const cv::Mat& img, int wait_ms = 1);
};
