#pragma once
#include <string>

struct CheckpointManifest {
    std::string dem_path;
    std::string nmea_path;
    std::string out_dir;
    double baro_alt_m   = -1.0;   // -1 = не задан → автооценка
    double start_lat    = 0.0;
    double start_lon    = 0.0;
    double speed_mps    = 40.0;
    double azimuth_deg  = 0.0;
    double search_radius_m = 3000.0;
    int    min_profile  = 35;
};

// Прочитать manifest.ini + heights.txt → flight.nmea, заполнить CheckpointManifest.
bool load_checkpoint_source(const std::string& source_dir, CheckpointManifest& out);
