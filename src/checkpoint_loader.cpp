#include "checkpoint_loader.hpp"
#include <fstream>
#include <iostream>
#include <filesystem>
#include <cmath>
#include <cctype>

namespace fs = std::filesystem;

static std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b-1]))) --b;
    return s.substr(a, b - a);
}

static std::pair<double,double> meters_to_latlon(double x_m, double y_m,
                                                  double ref_lat, double ref_lon) {
    constexpr double R = 6378137.0;
    double lat = ref_lat + y_m / R * (180.0 / M_PI);
    double lon = ref_lon + x_m / (R * std::cos(ref_lat * M_PI / 180.0)) * (180.0 / M_PI);
    return {lat, lon};
}

bool load_checkpoint_source(const std::string& source_dir, CheckpointManifest& out) {
    fs::path dir(source_dir);
    if (!fs::is_directory(dir)) {
        std::cerr << "ERROR: --source не каталог: " << source_dir << "\n";
        return false;
    }

    fs::path ini = dir / "manifest.ini";
    if (!fs::exists(ini)) {
        std::cerr << "ERROR: нет " << ini << "\n";
        return false;
    }

    std::ifstream in(ini);
    std::string line, section;
    auto kv = [&](const std::string& key, std::string& val) {
        auto pos = line.find('=');
        if (pos == std::string::npos) return;
        if (trim(line.substr(0, pos)) == key)
            val = trim(line.substr(pos + 1));
    };

    std::string map_f = "map.tif", heights_f = "heights.txt", out_f = "output";
    std::string ref_lat_s, ref_lon_s, start_x_s = "0", start_y_s = "0";
    std::string hdg_s, baro_s, spd_s, hz_s = "1", minp_s, radius_s;

    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#' || line[0] == ';') continue;
        if (line.front() == '[') { section = line; continue; }
        if (section != "[flight]") continue;
        kv("map", map_f);
        kv("heights", heights_f);
        kv("output", out_f);
        kv("origin_lat", ref_lat_s);
        kv("origin_lon", ref_lon_s);
        kv("start_x_m", start_x_s);
        kv("start_y_m", start_y_s);
        kv("heading_deg", hdg_s);
        kv("baro_alt_m", baro_s);
        kv("speed_mps", spd_s);
        kv("sample_hz", hz_s);
        kv("min_profile", minp_s);
        kv("search_radius_m", radius_s);
    }

    if (ref_lat_s.empty() || ref_lon_s.empty()) {
        std::cerr << "ERROR: manifest.ini: нужны origin_lat, origin_lon\n";
        return false;
    }

    double origin_lat = std::stod(ref_lat_s);
    double origin_lon = std::stod(ref_lon_s);
    double start_x_m  = std::stod(start_x_s);
    double start_y_m  = std::stod(start_y_s);
    double sample_hz  = std::stod(hz_s);

    out.dem_path = (dir / map_f).string();
    out.out_dir  = (dir / out_f).string();
    if (!baro_s.empty())   out.baro_alt_m = std::stod(baro_s);
    if (!spd_s.empty())    out.speed_mps  = std::stod(spd_s);
    if (!hdg_s.empty())    out.azimuth_deg= std::stod(hdg_s);
    if (!minp_s.empty())   out.min_profile = std::stoi(minp_s);
    if (!radius_s.empty()) out.search_radius_m = std::stod(radius_s);

    auto [lat, lon] = meters_to_latlon(start_x_m, start_y_m, origin_lat, origin_lon);
    out.start_lat = lat;
    out.start_lon = lon;

    fs::path heights_path = dir / heights_f;
    if (!fs::exists(heights_path)) {
        std::cerr << "ERROR: нет " << heights_path << "\n";
        return false;
    }

    fs::create_directories(out.out_dir);
    fs::path nmea_out = fs::path(out.out_dir) / "flight.nmea";
    out.nmea_path = nmea_out.string();

    std::ifstream hin(heights_path);
    std::ofstream nout(nmea_out);
    int i = 0, count = 0;
    double t0 = 36000.0, dt = 1.0 / sample_hz;

    while (std::getline(hin, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        double agl = std::stod(line);
        double ts = t0 + i * dt;
        int hh = static_cast<int>(ts / 3600) % 24;
        int mm = static_cast<int>(std::fmod(ts, 3600.0) / 60.0);
        double ss = std::fmod(ts, 60.0);

        char time_buf[16];
        std::snprintf(time_buf, sizeof(time_buf), "%02d%02d%06.3f", hh, mm, ss);

        std::string body = std::string("GPGGA,") + time_buf +
                           ",,,,,0,,," + std::to_string(agl) + ",M," +
                           std::to_string(out.baro_alt_m) + ",M,,";
        int cs = 0;
        for (char c : body) cs ^= static_cast<unsigned char>(c);
        char cs_buf[8];
        std::snprintf(cs_buf, sizeof(cs_buf), "%02X", cs & 0xFF);
        nout << "$" << body << "*" << cs_buf << "\n";
        ++i; ++count;
    }

    if (count == 0) {
        std::cerr << "ERROR: heights.txt пуст\n";
        return false;
    }

    std::cerr << "[SOURCE] Каталог: " << source_dir << "\n";
    std::cerr << "[SOURCE] heights.txt → " << nmea_out
              << " (" << count << " измерений, " << sample_hz << " Гц)\n";
    std::cerr << "[SOURCE] Старт: " << out.start_lat << "°N, "
              << out.start_lon << "°E, курс " << out.azimuth_deg << "°\n";
    return true;
}
