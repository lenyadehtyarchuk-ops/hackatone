#include <iostream>
#include <fstream>
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
#include "particle_filter.hpp"
#include "visualizer.hpp"

namespace fs = std::filesystem;

struct Config {
    std::string dem_path;
    std::string nmea_path;
    double      baro_alt_m      = 1500.0;
    double      start_lat       = 0.0;
    double      start_lon       = 0.0;
    double      speed_mps       = 80.0;
    double      azimuth_deg     = 120.0;  // начальный азимут для PF
    double      search_radius_m = 20000.0;
    int         min_profile     = 40;
    int         az_step         = 1;
    std::string out_dir         = "results";
    bool        show_gui        = false;
    bool        sliding         = false;
    bool        use_pf          = false;  // particle filter вместо TERCOM
    int         pf_n            = 1000;   // число частиц
    double      pf_hdg_noise    = 8.0;   // шум курса (°/шаг)
    double      pf_meas_sigma   = 4.0;   // σ измерения AGL для PF (чуть выше истинного)
    std::vector<double> jammer_zone;
};

static void print_usage(const char* prog) {
    std::cerr << "Использование:\n"
              << "  " << prog << " --dem <path.tif> --nmea <path.nmea>\n"
              << "               [--baro 1500] [--lat 55.0] [--lon 60.0]\n"
              << "               [--speed 80] [--azimuth 120]\n"
              << "               [--radius 20000] [--min-profile 40]\n"
              << "               [--az-step 1] [--out results/] [--gui] [--sliding]\n"
              << "               [--pf] [--pf-n 1000] [--pf-hdg-noise 8] [--pf-meas-sigma 4]\n";
}

static Config parse_args(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + a);
            return argv[++i];
        };
        if      (a == "--dem")            cfg.dem_path        = next();
        else if (a == "--nmea")           cfg.nmea_path       = next();
        else if (a == "--baro")           cfg.baro_alt_m      = std::stod(next());
        else if (a == "--lat")            cfg.start_lat       = std::stod(next());
        else if (a == "--lon")            cfg.start_lon       = std::stod(next());
        else if (a == "--speed")          cfg.speed_mps       = std::stod(next());
        else if (a == "--azimuth")        cfg.azimuth_deg     = std::stod(next());
        else if (a == "--radius")         cfg.search_radius_m = std::stod(next());
        else if (a == "--min-profile")    cfg.min_profile     = std::stoi(next());
        else if (a == "--az-step")        cfg.az_step         = std::stoi(next());
        else if (a == "--out")            cfg.out_dir         = next();
        else if (a == "--gui")            cfg.show_gui        = true;
        else if (a == "--sliding")        cfg.sliding         = true;
        else if (a == "--pf")             cfg.use_pf          = true;
        else if (a == "--pf-n")           cfg.pf_n            = std::stoi(next());
        else if (a == "--pf-hdg-noise")   cfg.pf_hdg_noise    = std::stod(next());
        else if (a == "--pf-meas-sigma")  cfg.pf_meas_sigma   = std::stod(next());
        else if (a == "--jammer-zone") {
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

// ─────────────────── Ветка Particle Filter ───────────────────────────────────

static int run_particle_filter(const Config& cfg, const DemData& dem,
                                const std::vector<NmeaFix>& fixes, double dt_s) {
    std::cerr << "[PF] Запуск Particle Filter: N=" << cfg.pf_n
              << " az0=" << cfg.azimuth_deg << "°"
              << " hdg_σ=" << cfg.pf_hdg_noise << "°/шаг"
              << " meas_σ=" << cfg.pf_meas_sigma << " м\n";

    // PF создаётся по первому переходу GPS→denied (из последней GPS-точки)
    std::unique_ptr<ParticleFilter> pf;

    std::vector<TrajectoryPoint> trajectory;
    trajectory.reserve(fixes.size());

    // CSV для fetch_satellite.py
    std::ofstream csv_out(cfg.out_dir + "/trn_estimate.csv");
    csv_out << "timestamp_s,found_lat,found_lon,heading_deg,speed_mps,neff\n";

    double last_gps_lat = cfg.start_lat, last_gps_lon = cfg.start_lon;
    double last_gps_hdg = cfg.azimuth_deg;
    bool   pf_ready     = false;

    // Кольцевой буфер последних win AGL-измерений, курсов и скоростей из GPS-фазы
    std::deque<double> gps_agl_buf, gps_hdg_buf, gps_spd_buf;
    const int win = cfg.min_profile;

    for (size_t fi = 0; fi < fixes.size(); ++fi) {
        const auto& fix = fixes[fi];
        bool gps_ok = (fix.gps_quality > 0);

        TrajectoryPoint tp;
        tp.gps_denied = !gps_ok;

        if (gps_ok) {
            // GPS доступен — берём координаты напрямую
            tp.lat         = fix.lat;
            tp.lon         = fix.lon;
            tp.speed_mps   = cfg.speed_mps;
            tp.heading_deg = last_gps_hdg;
            tp.ncc         = 1.0;
            // Вычислить heading из двух последних GPS-точек
            if (fi > 0 && fixes[fi-1].gps_quality > 0) {
                double dlat = fix.lat - fixes[fi-1].lat;
                double dlon = (fix.lon - fixes[fi-1].lon)
                              * std::cos(fix.lat * M_PI / 180.0);
                if (std::abs(dlat) + std::abs(dlon) > 1e-9)
                    last_gps_hdg = std::atan2(dlon, dlat) * 180.0 / M_PI;
            }
            last_gps_lat = fix.lat;
            last_gps_lon = fix.lon;
            pf_ready = false;

            // Накапливаем AGL, курс и скорость для предзаполнения PF
            gps_agl_buf.push_back(fix.radio_alt_m);
            gps_hdg_buf.push_front(last_gps_hdg);
            gps_spd_buf.push_front(cfg.speed_mps);  // [0]=последний, как в spd_history
            if ((int)gps_agl_buf.size() > win) gps_agl_buf.pop_front();
            if ((int)gps_hdg_buf.size() > win) gps_hdg_buf.pop_back();
            if ((int)gps_spd_buf.size() > win) gps_spd_buf.pop_back();
        } else {
            // GPS отсутствует — используем PF
            if (!pf_ready) {
                pf = std::make_unique<ParticleFilter>(
                        dem, last_gps_lat, last_gps_lon,
                        last_gps_hdg, cfg.speed_mps,
                        cfg.pf_n,
                        /*profile_win=*/win,
                        /*pos_sigma_m=*/10.0,
                        /*hdg_sigma_deg=*/5.0);
                // Предзаполнить буфер реальными GPS-данными — холодный старт устранён
                pf->prefill(gps_agl_buf, gps_hdg_buf, gps_spd_buf);
                pf_ready = true;
            }
            PFEstimate est = pf->step(fix.radio_alt_m, cfg.baro_alt_m, dt_s,
                                      cfg.pf_hdg_noise, 1.0, cfg.pf_meas_sigma);
            tp.lat        = est.lat;
            tp.lon        = est.lon;
            tp.heading_deg = est.heading_deg;
            tp.speed_mps  = est.speed_mps;
            tp.ncc        = est.neff / cfg.pf_n;
        }
        trajectory.push_back(tp);

        csv_out << std::fixed << std::setprecision(6)
                << fix.timestamp_s << ","
                << tp.lat << ","
                << tp.lon << ","
                << std::setprecision(2)
                << tp.heading_deg << ","
                << tp.speed_mps << ","
                << static_cast<int>(tp.ncc * cfg.pf_n) << "\n";
    }
    csv_out.close();

    // Визуализация
    if (!trajectory.empty()) {
        Visualizer::save_trajectory_on_dem(dem, trajectory,
                                           cfg.start_lat, cfg.start_lon,
                                           cfg.out_dir + "/trajectory.png",
                                           cfg.jammer_zone);
    }

    auto& last = trajectory.back();
    std::cerr << "\n=== РЕЗУЛЬТАТ (PF) ===\n"
              << "Позиция:   " << last.lat << "° N, " << last.lon << "° E\n"
              << "Вектор:    " << last.speed_mps << " м/с, "
              << last.heading_deg << "° (азимут)\n"
              << "Neff/N:    " << last.ncc << "\n"
              << "CSV:       " << cfg.out_dir << "/trn_estimate.csv\n"
              << "Траектория:" << cfg.out_dir << "/trajectory.png\n";
    return 0;
}

// ─────────────────── Ветка TERCOM (sliding NCC) ──────────────────────────────

static int run_tercom(const Config& cfg, const DemData& dem,
                      const std::vector<NmeaFix>& fixes, double dt_s) {
    CorrelatorConfig ccfg;
    ccfg.start_lat        = cfg.start_lat;
    ccfg.start_lon        = cfg.start_lon;
    ccfg.baro_alt_m       = cfg.baro_alt_m;
    ccfg.speed_mps        = cfg.speed_mps;
    ccfg.sample_dt_s      = dt_s;
    ccfg.min_profile_len  = cfg.min_profile;
    ccfg.search_radius_m  = cfg.search_radius_m;
    ccfg.azimuth_step_deg = cfg.az_step;
    ccfg.gradient_weight  = 0.4;
    ccfg.prior_heading_deg  = -1.0;
    ccfg.heading_search_deg = 45.0;

    TerrainCorrelator correlator(dem, ccfg);
    TrnKalmanFilter   kf(dt_s);

    std::vector<TrajectoryPoint>     trajectory;
    std::optional<CorrelationResult> last_result;

    double ref_lat = cfg.start_lat;
    double ref_lon = cfg.start_lon;

    std::ofstream csv_out(cfg.out_dir + "/trn_estimate.csv");
    csv_out << "timestamp_s,found_lat,found_lon,kf_x_m,kf_y_m,"
               "speed_mps,heading_deg,ncc\n";

    for (size_t i = 0; i < fixes.size(); ++i) {
        const auto& fix = fixes[i];
        kf.predict(dt_s);

        auto res = correlator.add_measurement(fix.radio_alt_m);
        if (!res) continue;
        last_result = res;

        auto [x_m, y_m] = latlon_to_meters(res->current_lat, res->current_lon,
                                            ref_lat, ref_lon);
        bool accept = true;
        if (kf.is_initialized()) {
            double dx = x_m - kf.x(), dy = y_m - kf.y();
            if (std::sqrt(dx*dx + dy*dy) / dt_s > 250.0) {
                accept = false;
                std::cerr << "[TERCOM] Отклонено: скорость "
                          << std::sqrt(dx*dx+dy*dy)/dt_s << " м/с\n";
            }
        }
        if (res->low_confidence) accept = false;

        double noise_m = 200.0 * (1.0 - std::max(0.0, res->best_ncc));
        if (accept) {
            if (!kf.is_initialized()) {
                double az_rad = res->best_azimuth_deg * M_PI / 180.0;
                kf.init(x_m, y_m,
                        cfg.speed_mps * std::sin(az_rad),
                        cfg.speed_mps * std::cos(az_rad));
            } else {
                kf.correct(x_m, y_m, noise_m);
            }
        }

        TrajectoryPoint tp;
        tp.lat        = res->current_lat;
        tp.lon        = res->current_lon;
        tp.speed_mps  = kf.is_initialized() ? kf.speed_mps() : res->ground_speed_mps;
        tp.heading_deg = kf.is_initialized() ? kf.heading_deg()
                                             : static_cast<double>(res->best_azimuth_deg);
        tp.ncc        = res->best_ncc;
        tp.gps_denied = (fix.gps_quality == 0);
        trajectory.push_back(tp);

        csv_out << std::fixed << std::setprecision(6)
                << fix.timestamp_s << ","
                << res->current_lat << ","
                << res->current_lon << ","
                << kf.x() << "," << kf.y() << ","
                << std::setprecision(2) << tp.speed_mps << ","
                << tp.heading_deg << ","
                << std::setprecision(4) << res->best_ncc << "\n";

        std::cerr << "[CORR] az=" << res->best_azimuth_deg << "°"
                  << " NCC=" << res->best_ncc
                  << " σ_profile=" << res->profile_std << " м"
                  << " cur=(" << res->current_lat << ", " << res->current_lon << ")\n";

        if (cfg.sliding && accept) {
            ccfg.start_lat = res->current_lat;
            ccfg.start_lon = res->current_lon;
            ccfg.prior_heading_deg = (kf.is_initialized() && kf.speed_mps() > 5.0)
                                   ? kf.heading_deg()
                                   : static_cast<double>(res->best_azimuth_deg);
            ccfg.min_profile_len = std::min(cfg.min_profile, 40);
            correlator = TerrainCorrelator(dem, ccfg);
        }
    }

    if (last_result) {
        Visualizer::save_correlation_heatmap(
            last_result->corr_map,
            cfg.out_dir + "/correlation_heatmap.png",
            last_result->best_azimuth_deg / cfg.az_step,
            static_cast<int>(last_result->best_offset_m / (cfg.speed_mps * dt_s)));
    }

    if (!trajectory.empty()) {
        Visualizer::save_trajectory_on_dem(dem, trajectory,
                                           cfg.start_lat, cfg.start_lon,
                                           cfg.out_dir + "/trajectory.png",
                                           cfg.jammer_zone);
    }

    if (trajectory.empty()) {
        std::cerr << "[TERCOM] Нет данных для корреляции.\n";
        return 2;
    }

    auto& last = trajectory.back();
    std::cerr << "\n=== РЕЗУЛЬТАТ (TERCOM) ===\n"
              << "Найденная позиция: " << last.lat << "° N, " << last.lon << "° E\n"
              << "Вектор скорости:   " << last.speed_mps << " м/с, "
              << last.heading_deg << "° (азимут)\n"
              << "NCC (качество):    " << last.ncc << "\n"
              << "Файлы сохранены в: " << cfg.out_dir << "/\n";
    return 0;
}

// ─────────────────── main ────────────────────────────────────────────────────

int main(int argc, char** argv) {
    Config cfg = parse_args(argc, argv);
    fs::create_directories(cfg.out_dir);

    std::cerr << "[MAIN] Загрузка ЦМР: " << cfg.dem_path << "\n";
    DemData dem = DemLoader::load(cfg.dem_path);

    if (cfg.start_lat == 0.0 && cfg.start_lon == 0.0) {
        cfg.start_lat = dem.gt.origin_lat + dem.elev.rows * dem.gt.pixel_h / 2.0;
        cfg.start_lon = dem.gt.origin_lon + dem.elev.cols * dem.gt.pixel_w / 2.0;
        std::cerr << "[MAIN] Стартовая точка (центр карты): "
                  << cfg.start_lat << ", " << cfg.start_lon << "\n";
    }

    std::cerr << "[MAIN] Загрузка NMEA: " << cfg.nmea_path << "\n";
    auto fixes = NmeaParser::load_file(cfg.nmea_path);
    if (fixes.empty()) {
        std::cerr << "ERROR: NMEA файл пустой\n";
        return 1;
    }

    double dt_s = 1.0;
    if (fixes.size() > 1) {
        double total_t = fixes.back().timestamp_s - fixes.front().timestamp_s;
        if (total_t > 0) dt_s = total_t / (fixes.size() - 1);
    }
    std::cerr << "[MAIN] dt=" << dt_s << " с, " << fixes.size() << " измерений\n";

    if (cfg.use_pf)
        return run_particle_filter(cfg, dem, fixes, dt_s);
    else
        return run_tercom(cfg, dem, fixes, dt_s);
}
