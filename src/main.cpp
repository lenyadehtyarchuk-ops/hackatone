#include <iostream>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stdexcept>

#include "nmea_parser.hpp"
#include "dem_loader.hpp"
#include "terrain_correlator.hpp"
#include "kalman_filter.hpp"
#include "visualizer.hpp"

namespace fs = std::filesystem;

struct Config {
    std::string dem_path;
    std::string nmea_path;
    double      baro_alt_m      = 1500.0;
    double      start_lat       = 0.0;
    double      start_lon       = 0.0;
    double      speed_mps       = 80.0;
    double      search_radius_m = 20000.0;
    int         min_profile     = 40;
    int         az_step         = 1;
    std::string out_dir         = "results";
    bool        show_gui        = false;
    bool        sliding         = false;
    // Зона глушения GPS: lat1,lon1,lat2,lon2 (пустая = нет зоны)
    std::vector<double> jammer_zone;
};

static void print_usage(const char* prog) {
    std::cerr << "Использование:\n"
              << "  " << prog << " --dem <path.tif> --nmea <path.nmea>\n"
              << "               [--baro 1500] [--lat 55.0] [--lon 60.0]\n"
              << "               [--speed 80] [--radius 20000] [--min-profile 40]\n"
              << "               [--az-step 1] [--out results/] [--gui] [--sliding]\n";
}

static Config parse_args(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + a);
            return argv[++i];
        };
        if      (a == "--dem")          cfg.dem_path        = next();
        else if (a == "--nmea")         cfg.nmea_path       = next();
        else if (a == "--baro")         cfg.baro_alt_m      = std::stod(next());
        else if (a == "--lat")          cfg.start_lat       = std::stod(next());
        else if (a == "--lon")          cfg.start_lon       = std::stod(next());
        else if (a == "--speed")        cfg.speed_mps       = std::stod(next());
        else if (a == "--radius")       cfg.search_radius_m = std::stod(next());
        else if (a == "--min-profile")  cfg.min_profile     = std::stoi(next());
        else if (a == "--az-step")      cfg.az_step         = std::stoi(next());
        else if (a == "--out")          cfg.out_dir         = next();
        else if (a == "--gui")          cfg.show_gui        = true;
        else if (a == "--sliding")      cfg.sliding         = true;
        else if (a == "--jammer-zone") {
            // Формат: lat1,lon1,lat2,lon2
            std::string val = next();
            std::replace(val.begin(), val.end(), ',', ' ');
            std::istringstream ss(val);
            double v;
            while (ss >> v) cfg.jammer_zone.push_back(v);
            if (cfg.jammer_zone.size() != 4) {
                std::cerr << "ERROR: --jammer-zone требует 4 числа: lat1,lon1,lat2,lon2\n";
                std::exit(1);
            }
        }
        else if (a == "--help" || a == "-h") { print_usage(argv[0]); std::exit(0); }
        else { std::cerr << "Unknown arg: " << a << "\n"; print_usage(argv[0]); std::exit(1); }
    }
    if (cfg.dem_path.empty() || cfg.nmea_path.empty()) {
        print_usage(argv[0]);
        std::exit(1);
    }
    return cfg;
}

// Перевести (lat, lon) относительно стартовой точки в метры
static std::pair<double,double> latlon_to_meters(double lat, double lon,
                                                  double ref_lat, double ref_lon) {
    constexpr double R = 6378137.0;
    double x = (lon - ref_lon) * (M_PI/180.0) * R * std::cos(ref_lat * M_PI/180.0);
    double y = (lat - ref_lat) * (M_PI/180.0) * R;
    return {x, y};
}

int main(int argc, char** argv) {
    Config cfg = parse_args(argc, argv);

    fs::create_directories(cfg.out_dir);

    // Загрузить ЦМР
    std::cerr << "[MAIN] Загрузка ЦМР: " << cfg.dem_path << "\n";
    DemData dem = DemLoader::load(cfg.dem_path);

    // Если стартовые координаты не заданы — брать центр карты
    if (cfg.start_lat == 0.0 && cfg.start_lon == 0.0) {
        cfg.start_lat = dem.gt.origin_lat + dem.elev.rows * dem.gt.pixel_h / 2.0;
        cfg.start_lon = dem.gt.origin_lon + dem.elev.cols * dem.gt.pixel_w / 2.0;
        std::cerr << "[MAIN] Стартовая точка (центр карты): "
                  << cfg.start_lat << ", " << cfg.start_lon << "\n";
    }

    // Загрузить NMEA
    std::cerr << "[MAIN] Загрузка NMEA: " << cfg.nmea_path << "\n";
    auto fixes = NmeaParser::load_file(cfg.nmea_path);
    if (fixes.empty()) {
        std::cerr << "ERROR: NMEA файл пустой или не содержит GPGGA сообщений\n";
        return 1;
    }

    // Вычислить шаг дискретизации по временным меткам
    double dt_s = 1.0;
    if (fixes.size() > 1) {
        double total_t = fixes.back().timestamp_s - fixes.front().timestamp_s;
        if (total_t > 0)
            dt_s = total_t / (fixes.size() - 1);
    }
    std::cerr << "[MAIN] dt=" << dt_s << " с, " << fixes.size() << " измерений\n";

    // Инициализировать коррелятор
    CorrelatorConfig ccfg;
    ccfg.start_lat        = cfg.start_lat;
    ccfg.start_lon        = cfg.start_lon;
    ccfg.baro_alt_m       = cfg.baro_alt_m;
    ccfg.speed_mps        = cfg.speed_mps;
    ccfg.sample_dt_s      = dt_s;
    ccfg.min_profile_len  = cfg.min_profile;
    ccfg.search_radius_m  = cfg.search_radius_m;
    ccfg.azimuth_step_deg = cfg.az_step;
    ccfg.gradient_weight  = 0.4;   // 40% градиент, 60% высота
    ccfg.prior_heading_deg = -1.0; // -1 = полный 360° поиск до первого решения
    ccfg.heading_search_deg = 45.0;

    TerrainCorrelator correlator(dem, ccfg);

    // Фильтр Калмана в локальных метрах от стартовой точки
    TrnKalmanFilter kf(dt_s);

    std::vector<TrajectoryPoint> trajectory;
    std::optional<CorrelationResult> last_result;

    double ref_lat = cfg.start_lat;
    double ref_lon = cfg.start_lon;

    std::cout << "timestamp_s,found_lat,found_lon,kf_x_m,kf_y_m,"
                 "speed_mps,heading_deg,ncc\n";

    for (size_t i = 0; i < fixes.size(); ++i) {
        const auto& fix = fixes[i];

        kf.predict(dt_s);

        auto res = correlator.add_measurement(fix.radio_alt_m);
        if (!res) continue;

        last_result = res;

        // Текущая позиция ВС (конец профиля)
        auto [x_m, y_m] = latlon_to_meters(res->current_lat, res->current_lon,
                                            ref_lat, ref_lon);

        // Отклонить физически нереалистичные скачки позиции (> 250 м/с между шагами)
        bool accept = true;
        if (kf.is_initialized()) {
            double dx = x_m - kf.x(), dy = y_m - kf.y();
            double implied_speed = std::sqrt(dx*dx + dy*dy) / dt_s;
            if (implied_speed > 250.0) {
                accept = false;
                std::cerr << "[MAIN] Отклонено: implied_speed=" << implied_speed << " м/с\n";
            }
        }

        // Низкая уверенность → не обновлять Калман, держать мёртвое счисление
        if (res->low_confidence) {
            accept = false;
            std::cerr << "[MAIN] Низкая уверенность: NCC=" << res->best_ncc
                      << " σ_рельефа=" << res->profile_std
                      << " м — мёртвое счисление\n";
        }

        // Оценка точности по NCC: чем выше NCC, тем меньше шум измерения
        double noise_m = 200.0 * (1.0 - std::max(0.0, res->best_ncc));
        if (accept) kf.correct(x_m, y_m, noise_m);

        TrajectoryPoint tp;
        tp.lat         = res->current_lat;
        tp.lon         = res->current_lon;
        tp.speed_mps   = kf.is_initialized() ? kf.speed_mps() : res->ground_speed_mps;
        tp.heading_deg = kf.is_initialized() ? kf.heading_deg()
                                              : static_cast<double>(res->best_azimuth_deg);
        tp.ncc         = res->best_ncc;
        tp.gps_denied  = (fix.gps_quality == 0);
        trajectory.push_back(tp);

        std::cout << std::fixed << std::setprecision(6)
                  << fix.timestamp_s << ","
                  << res->current_lat << ","
                  << res->current_lon << ","
                  << kf.x() << ","
                  << kf.y() << ","
                  << std::setprecision(2)
                  << tp.speed_mps << ","
                  << tp.heading_deg << ","
                  << std::setprecision(4)
                  << res->best_ncc << "\n";

        std::cerr << "[MAIN] t=" << fix.timestamp_s
                  << " az=" << res->best_azimuth_deg << "°"
                  << " v=" << tp.speed_mps << " м/с"
                  << " hdg=" << tp.heading_deg << "°"
                  << " NCC=" << res->best_ncc << "\n";

        // В скользящем режиме сдвигаем стартовую точку только при принятом решении.
        // При dead reckoning оставляем последнюю надёжную стартовую точку.
        if (cfg.sliding && accept) {
            ccfg.start_lat = res->current_lat;
            ccfg.start_lon = res->current_lon;
            if (kf.is_initialized())
                ccfg.prior_heading_deg = kf.heading_deg();
            correlator = TerrainCorrelator(dem, ccfg);
        }
    }

    // Визуализация
    if (last_result) {
        std::string heatmap_path = cfg.out_dir + "/correlation_heatmap.png";
        Visualizer::save_correlation_heatmap(
            last_result->corr_map, heatmap_path,
            last_result->best_azimuth_deg / cfg.az_step,
            static_cast<int>(last_result->best_offset_m /
                             (cfg.speed_mps * dt_s)));

        if (!trajectory.empty()) {
            std::string traj_path = cfg.out_dir + "/trajectory.png";
            Visualizer::save_trajectory_on_dem(dem, trajectory,
                                                cfg.start_lat, cfg.start_lon,
                                                traj_path, cfg.jammer_zone);
        }
    }

    if (trajectory.empty()) {
        std::cerr << "[MAIN] Недостаточно данных для корреляции. "
                     "Увеличьте длительность полёта или уменьшите --min-profile.\n";
        return 2;
    }

    // Итоговый вывод
    auto& last = trajectory.back();
    std::cerr << "\n=== РЕЗУЛЬТАТ ===\n"
              << "Найденная позиция: " << last.lat << "° N, " << last.lon << "° E\n"
              << "Вектор скорости:   " << last.speed_mps << " м/с, "
              << last.heading_deg << "° (азимут)\n"
              << "NCC (качество):    " << last.ncc << "\n"
              << "Файлы сохранены в: " << cfg.out_dir << "/\n";

    return 0;
}
