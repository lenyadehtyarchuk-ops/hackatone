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
#include "contour_pf.hpp"
#include "keypoint_db.hpp"
#include "ridge_map.hpp"
#include "line_navigator.hpp"
#include "visualizer.hpp"
#include "checkpoint_loader.hpp"

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
    bool        use_pf          = false;  // profile particle filter
    bool        use_cpf         = false;  // contour particle filter (новый, быстрый)
    int         pf_n            = 1000;   // число частиц для PF
    double      pf_hdg_noise    = 8.0;   // шум курса (°/шаг)
    double      pf_meas_sigma   = 4.0;   // σ для PF (профильный)
    int         cpf_n           = 3000;  // число частиц для CPF
    double      cpf_sigma       = 15.0;  // σ измерения для CPF (м)
    int         cpf_kp_window   = 15;   // окно поиска экстремумов ЦМР (пиксели)
    bool        use_lnav        = false; // line navigator (хребты/впадины)
    bool        use_combined    = false; // CPF + RidgeMap ridge fixes
    float       lnav_gauss      = 3.0f; // Gaussian sigma для Hessian (пикселей)
    float       lnav_thresh     = 0.03f;// порог кривизны (доля от максимума)
    std::vector<double> jammer_zone;
    std::string source_dir;
};

static void apply_checkpoint(const CheckpointManifest& m, Config& cfg) {
    cfg.dem_path        = m.dem_path;
    cfg.nmea_path       = m.nmea_path;
    cfg.out_dir         = m.out_dir;
    cfg.baro_alt_m      = m.baro_alt_m;
    cfg.start_lat       = m.start_lat;
    cfg.start_lon       = m.start_lon;
    cfg.speed_mps       = m.speed_mps;
    cfg.azimuth_deg     = m.azimuth_deg;
    cfg.search_radius_m = m.search_radius_m;
    cfg.min_profile     = m.min_profile;
}

static void print_usage(const char* prog) {
    std::cerr << "Использование:\n"
              << "  " << prog << " --dem <path.tif> --nmea <path.nmea>\n"
              << "               [--baro 1500] [--lat 55.0] [--lon 60.0]\n"
              << "               [--speed 80] [--azimuth 120]\n"
              << "               [--radius 20000] [--min-profile 40]\n"
              << "               [--az-step 1] [--out results/] [--gui] [--sliding]\n"
              << "               [--source test_source]  # формат 3-го чекпоинта\n"
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
        else if (a == "--cpf")            cfg.use_cpf         = true;
        else if (a == "--cpf-n")          cfg.cpf_n           = std::stoi(next());
        else if (a == "--cpf-sigma")      cfg.cpf_sigma       = std::stod(next());
        else if (a == "--cpf-kp-window")  cfg.cpf_kp_window   = std::stoi(next());
        else if (a == "--lnav")           cfg.use_lnav        = true;
        else if (a == "--combined")       cfg.use_combined    = true;
        else if (a == "--lnav-gauss")     cfg.lnav_gauss      = std::stof(next());
        else if (a == "--lnav-thresh")    cfg.lnav_thresh     = std::stof(next());
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
        else if (a == "--source")       cfg.source_dir      = next();
        else if (a == "--help" || a == "-h") { print_usage(argv[0]); std::exit(0); }
        else { std::cerr << "Unknown arg: " << a << "\n"; print_usage(argv[0]); std::exit(1); }
    }
    if (!cfg.source_dir.empty()) {
        CheckpointManifest m;
        if (!load_checkpoint_source(cfg.source_dir, m)) std::exit(1);
        apply_checkpoint(m, cfg);
    } else if (cfg.dem_path.empty() || cfg.nmea_path.empty()) {
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
    csv_out << "timestamp_s,found_lat,found_lon,x_m,y_m,"
               "speed_mps,heading_deg,ncc\n";

    for (size_t i = 0; i < fixes.size(); ++i) {
        const auto& fix = fixes[i];
        bool gps_ok = (fix.gps_quality > 0);
        kf.predict(dt_s);

        auto res = correlator.add_measurement(fix.radio_alt_m);

        // GPS-точки всегда рисуем из реальных координат
        if (gps_ok) {
            TrajectoryPoint tp;
            tp.lat        = fix.lat;
            tp.lon        = fix.lon;
            tp.speed_mps  = cfg.speed_mps;
            tp.heading_deg = cfg.azimuth_deg;
            tp.ncc        = 1.0;
            tp.gps_denied = false;
            trajectory.push_back(tp);
            if (res) last_result = res;
            continue;
        }

        // GPS-denied: нужна корреляция
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
        tp.gps_denied = true;
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
        float offset_step_m = static_cast<float>(cfg.speed_mps * dt_s);
        Visualizer::save_correlation_heatmap(
            last_result->corr_map,
            cfg.out_dir + "/correlation_heatmap.png",
            last_result->best_azimuth_deg / cfg.az_step,
            static_cast<int>(last_result->best_offset_m / offset_step_m),
            offset_step_m,
            static_cast<float>(last_result->best_ncc));
    }

    if (!trajectory.empty()) {
        double faz  = last_result ? last_result->best_azimuth_deg : -1.0;
        double fdist = last_result ? last_result->best_offset_m   : 0.0;
        Visualizer::save_trajectory_on_dem(dem, trajectory,
                                           cfg.start_lat, cfg.start_lon,
                                           cfg.out_dir + "/trajectory.png",
                                           cfg.jammer_zone,
                                           faz, fdist, cfg.speed_mps);
    }

    if (trajectory.empty()) {
        std::cerr << "[TERCOM] Нет данных для корреляции.\n";
        return 2;
    }

    auto& last = trajectory.back();
    int n_corr = static_cast<int>(trajectory.size());
    double speed_kmh = last.speed_mps * 3.6;

    // Расстояние от старта до последней точки через координаты
    double dlat_km = (last.lat - cfg.start_lat) * 111.32;
    double dlon_km = (last.lon - cfg.start_lon) * 111.32
                   * std::cos(cfg.start_lat * M_PI / 180.0);
    double total_dist_km = std::sqrt(dlat_km*dlat_km + dlon_km*dlon_km);

    std::cerr << std::fixed << std::setprecision(6)
              << "\n=== РЕЗУЛЬТАТ (TERCOM) ===\n"
              << "Найденная позиция:  " << last.lat << "° N, " << last.lon << "° E\n";
    if (last_result) {
        std::cerr << std::fixed << std::setprecision(1)
                  << "Путевой угол:       " << last_result->best_azimuth_deg << "°\n"
                  << "Путевая скорость:   " << last.speed_mps << " м/с"
                  << "  (" << speed_kmh << " км/ч)\n"
                  << std::setprecision(4)
                  << "NCC (качество):     " << last_result->best_ncc << "\n"
                  << std::fixed << std::setprecision(2)
                  << "Пройденное расст.:  " << total_dist_km << " км\n";
    }
    std::cerr << "Корреляций всего:   " << n_corr << "\n"
              << "Файлы:              " << cfg.out_dir << "/\n";
    return 0;
}

// ─────────────────── Ветка Line Navigator (хребты и впадины) ─────────────────

static int run_line_navigator(const Config& cfg, const DemData& dem,
                               const std::vector<NmeaFix>& fixes, double dt_s) {
    std::cerr << "[LNAV] Вычисление карты хребтов/впадин (Hessian, NMS)...\n";
    RidgeMap rm(dem, cfg.lnav_gauss, cfg.lnav_thresh);
    std::cerr << "[LNAV] Карта готова: " << rm.size() << " точек\n";

    std::vector<TrajectoryPoint> trajectory;
    trajectory.reserve(fixes.size());

    std::ofstream csv_out(cfg.out_dir + "/trn_estimate.csv");
    csv_out << "timestamp_s,found_lat,found_lon,heading_deg,speed_mps,fix\n";

    double last_gps_lat = cfg.start_lat, last_gps_lon = cfg.start_lon;
    double last_gps_hdg = cfg.azimuth_deg;

    std::unique_ptr<LineNavigator> lnav;
    bool lnav_ready = false;

    for (size_t fi = 0; fi < fixes.size(); ++fi) {
        const auto& fix = fixes[fi];
        bool gps_ok = (fix.gps_quality > 0);

        TrajectoryPoint tp;
        tp.gps_denied = !gps_ok;

        if (gps_ok) {
            tp.lat         = fix.lat;
            tp.lon         = fix.lon;
            tp.speed_mps   = cfg.speed_mps;
            tp.heading_deg = last_gps_hdg;
            tp.ncc         = 1.0;
            if (fi > 0 && fixes[fi-1].gps_quality > 0) {
                double dlat = fix.lat - fixes[fi-1].lat;
                double dlon = (fix.lon - fixes[fi-1].lon)
                              * std::cos(fix.lat * M_PI / 180.0);
                if (std::abs(dlat) + std::abs(dlon) > 1e-9)
                    last_gps_hdg = std::atan2(dlon, dlat) * 180.0 / M_PI;
            }
            last_gps_lat = fix.lat;
            last_gps_lon = fix.lon;
            lnav_ready   = false;
        } else {
            if (!lnav_ready) {
                LNavConfig lcfg;
                lcfg.elev_tol_m = static_cast<float>(cfg.cpf_sigma) * 2.0f;
                lnav = std::make_unique<LineNavigator>(
                        rm, last_gps_lat, last_gps_lon,
                        last_gps_hdg, cfg.speed_mps, lcfg);
                lnav_ready = true;
            }
            auto est = lnav->step(fix.radio_alt_m, cfg.baro_alt_m, dt_s,
                                  cfg.pf_hdg_noise, 1.0);
            tp.lat         = est.lat;
            tp.lon         = est.lon;
            tp.heading_deg = est.heading_deg;
            tp.speed_mps   = est.speed_mps;
            tp.ncc         = est.fix_applied ? 1.0 : 0.5;
        }
        trajectory.push_back(tp);

        csv_out << std::fixed << std::setprecision(6)
                << fix.timestamp_s << ","
                << tp.lat << "," << tp.lon << ","
                << std::setprecision(2)
                << tp.heading_deg << "," << tp.speed_mps << ","
                << (tp.ncc > 0.9 ? 1 : 0) << "\n";
    }
    csv_out.close();

    if (!trajectory.empty()) {
        Visualizer::save_trajectory_on_dem(dem, trajectory,
                                           cfg.start_lat, cfg.start_lon,
                                           cfg.out_dir + "/trajectory.png",
                                           cfg.jammer_zone);
    }

    auto& last = trajectory.back();
    int total = lnav ? lnav->current().total_fixes : 0;
    std::cerr << "\n=== РЕЗУЛЬТАТ (LNAV) ===\n"
              << "Позиция:   " << last.lat << "° N, " << last.lon << "° E\n"
              << "Вектор:    " << last.speed_mps << " м/с, "
              << last.heading_deg << "° (азимут)\n"
              << "Всего fix: " << total << "\n"
              << "Карта:     " << cfg.out_dir << "/trajectory.png\n";
    return 0;
}

// ─────────────────── Ветка Contour PF (ключевые точки + изолинии) ────────────

static int run_contour_pf(const Config& cfg, const DemData& dem,
                           const std::vector<NmeaFix>& fixes, double dt_s) {
    std::cerr << "[CPF] Загрузка ключевых точек рельефа (W=" << cfg.cpf_kp_window << ")...\n";
    KeypointDB kpdb(dem, cfg.cpf_kp_window, 30.0f);
    std::cerr << "[CPF] Запуск ContourPF: N=" << cfg.cpf_n
              << " σ=" << cfg.cpf_sigma << " м"
              << " az0=" << cfg.azimuth_deg << "°\n";

    std::unique_ptr<ContourPF> cpf;

    std::vector<TrajectoryPoint> trajectory;
    trajectory.reserve(fixes.size());

    std::ofstream csv_out(cfg.out_dir + "/trn_estimate.csv");
    csv_out << "timestamp_s,found_lat,found_lon,heading_deg,speed_mps,neff\n";

    double last_gps_lat = cfg.start_lat, last_gps_lon = cfg.start_lon;
    double last_gps_hdg = cfg.azimuth_deg;
    bool   cpf_ready    = false;

    for (size_t fi = 0; fi < fixes.size(); ++fi) {
        const auto& fix = fixes[fi];
        bool gps_ok = (fix.gps_quality > 0);

        TrajectoryPoint tp;
        tp.gps_denied = !gps_ok;

        if (gps_ok) {
            tp.lat         = fix.lat;
            tp.lon         = fix.lon;
            tp.speed_mps   = cfg.speed_mps;
            tp.heading_deg = last_gps_hdg;
            tp.ncc         = 1.0;
            if (fi > 0 && fixes[fi-1].gps_quality > 0) {
                double dlat = fix.lat - fixes[fi-1].lat;
                double dlon = (fix.lon - fixes[fi-1].lon)
                              * std::cos(fix.lat * M_PI / 180.0);
                if (std::abs(dlat) + std::abs(dlon) > 1e-9)
                    last_gps_hdg = std::atan2(dlon, dlat) * 180.0 / M_PI;
            }
            last_gps_lat = fix.lat;
            last_gps_lon = fix.lon;
            cpf_ready    = false;
        } else {
            if (!cpf_ready) {
                cpf = std::make_unique<ContourPF>(
                        dem, kpdb,
                        last_gps_lat, last_gps_lon,
                        last_gps_hdg, cfg.speed_mps,
                        cfg.cpf_n,
                        /*pos_sigma_m=*/10.0,
                        /*hdg_sigma_deg=*/5.0);
                cpf_ready = true;
            }
            auto est = cpf->step(fix.radio_alt_m, cfg.baro_alt_m, dt_s,
                                 cfg.pf_hdg_noise, 1.0, cfg.cpf_sigma);
            tp.lat         = est.lat;
            tp.lon         = est.lon;
            tp.heading_deg = est.heading_deg;
            tp.speed_mps   = est.speed_mps;
            tp.ncc         = est.neff / cfg.cpf_n;
        }
        trajectory.push_back(tp);

        csv_out << std::fixed << std::setprecision(6)
                << fix.timestamp_s << ","
                << tp.lat << "," << tp.lon << ","
                << std::setprecision(2)
                << tp.heading_deg << "," << tp.speed_mps << ","
                << static_cast<int>(tp.ncc * cfg.cpf_n) << "\n";
    }
    csv_out.close();

    if (!trajectory.empty()) {
        Visualizer::save_trajectory_on_dem(dem, trajectory,
                                           cfg.start_lat, cfg.start_lon,
                                           cfg.out_dir + "/trajectory.png",
                                           cfg.jammer_zone);
    }

    auto& last = trajectory.back();
    std::cerr << "\n=== РЕЗУЛЬТАТ (CPF) ===\n"
              << "Позиция:   " << last.lat << "° N, " << last.lon << "° E\n"
              << "Вектор:    " << last.speed_mps << " м/с, "
              << last.heading_deg << "° (азимут)\n"
              << "Neff/N:    " << last.ncc << "\n"
              << "CSV:       " << cfg.out_dir << "/trn_estimate.csv\n"
              << "Траектория:" << cfg.out_dir << "/trajectory.png\n";
    return 0;
}

// ─────────────────── Ветка Combined (CPF + RidgeMap ridge fixes) ─────────────

static int run_combined(const Config& cfg, const DemData& dem,
                         const std::vector<NmeaFix>& fixes, double dt_s) {
    std::cerr << "[COMBINED] Построение RidgeMap (Hessian, gauss="
              << cfg.lnav_gauss << ", thresh=" << cfg.lnav_thresh << ")...\n";
    RidgeMap rm(dem, cfg.lnav_gauss, cfg.lnav_thresh);
    std::cerr << "[COMBINED] RidgeMap готов: " << rm.size() << " точек\n";

    std::cerr << "[COMBINED] Загрузка KeypointDB (W=" << cfg.cpf_kp_window << ")...\n";
    KeypointDB kpdb(dem, cfg.cpf_kp_window, 30.0f);

    std::cerr << "[COMBINED] ContourPF N=" << cfg.cpf_n
              << " σ=" << cfg.cpf_sigma << " м (с ridge-fixes)\n";

    std::unique_ptr<ContourPF> cpf;

    std::vector<TrajectoryPoint> trajectory;
    trajectory.reserve(fixes.size());

    std::ofstream csv_out(cfg.out_dir + "/trn_estimate.csv");
    csv_out << "timestamp_s,found_lat,found_lon,x_m,y_m,heading_deg,speed_mps,neff\n";

    double last_gps_lat = cfg.start_lat, last_gps_lon = cfg.start_lon;
    double last_gps_hdg = cfg.azimuth_deg;
    bool   cpf_ready    = false;

    std::deque<double> gps_hdg_history;
    constexpr int HDG_WIN = 10;

    double active_cpf_sigma = cfg.cpf_sigma;
    double active_hdg_noise = cfg.pf_hdg_noise;

    // Dead-reckoning state for flat-terrain fallback
    double dr_lat = 0.0, dr_lon = 0.0, dr_hdg = 0.0, dr_speed = 0.0;
    double dr_yaw_rate = 0.0;  // °/s from GPS history — kept fixed once set

    for (size_t fi = 0; fi < fixes.size(); ++fi) {
        const auto& fix = fixes[fi];
        bool gps_ok = (fix.gps_quality > 0);

        TrajectoryPoint tp;
        tp.gps_denied = !gps_ok;

        if (gps_ok) {
            tp.lat         = fix.lat;
            tp.lon         = fix.lon;
            tp.speed_mps   = cfg.speed_mps;
            tp.heading_deg = last_gps_hdg;
            tp.ncc         = 1.0;
            if (fi > 0 && fixes[fi-1].gps_quality > 0) {
                double dlat = fix.lat - fixes[fi-1].lat;
                double dlon = (fix.lon - fixes[fi-1].lon)
                              * std::cos(fix.lat * M_PI / 180.0);
                if (std::abs(dlat) + std::abs(dlon) > 1e-9)
                    last_gps_hdg = std::atan2(dlon, dlat) * 180.0 / M_PI;
            }
            gps_hdg_history.push_back(last_gps_hdg);
            if ((int)gps_hdg_history.size() > HDG_WIN) gps_hdg_history.pop_front();
            last_gps_lat = fix.lat;
            last_gps_lon = fix.lon;
            cpf_ready    = false;
        } else {
            if (!cpf_ready) {
                // Compute yaw_rate from GPS heading history
                double yaw_rate = 0.0;
                int nh = (int)gps_hdg_history.size();
                if (nh >= 2) {
                    double sum_dhdg = 0;
                    for (int k = 1; k < nh; k++) {
                        double d = gps_hdg_history[k] - gps_hdg_history[k-1];
                        while (d >  180) d -= 360;
                        while (d < -180) d += 360;
                        sum_dhdg += d;
                    }
                    yaw_rate = sum_dhdg / ((nh - 1) * dt_s);
                }
                dr_yaw_rate = yaw_rate;

                // Detect flat terrain by sampling DEM roughness around entry point.
                // Variance of pre-jammer AGL readings is unreliable (slow approach
                // to mountains also looks flat). DEM sampling is ground truth.
                double local_sigma = 0.0;
                {
                    double dlat = 3000.0 / 111320.0;
                    double cos_l = std::cos(last_gps_lat * M_PI / 180.0);
                    double dlon  = 3000.0 / (111320.0 * cos_l);
                    constexpr int S = 8;
                    double sum = 0, sum2 = 0;
                    int cnt = 0;
                    for (int i = -S; i <= S; i++) {
                        for (int j = -S; j <= S; j++) {
                            float h = DemLoader::sample(dem,
                                (float)(last_gps_lat + i * dlat / S),
                                (float)(last_gps_lon + j * dlon / S), -9999.0f);
                            if (h > -9000.0f) { sum += h; sum2 += h * h; cnt++; }
                        }
                    }
                    if (cnt > 4) {
                        double m = sum / cnt;
                        local_sigma = std::sqrt(std::max(0.0, sum2 / cnt - m * m));
                    }
                }

                constexpr double FLAT_THRESH = 50.0;  // m RMS
                bool flat = (local_sigma < FLAT_THRESH);

                double eff_hdg_sigma = flat ? 2.0  : 5.0;
                double eff_yaw_sigma = flat ? 0.02 : 0.2;
                active_cpf_sigma     = flat ? 50.0 : cfg.cpf_sigma;
                active_hdg_noise     = flat ? 0.3  : cfg.pf_hdg_noise;

                std::cerr << "[COMBINED] local_sigma=" << local_sigma
                          << (flat ? " → FLAT (guided DR)" : " → MOUNTAIN")
                          << " yaw_rate=" << yaw_rate << "°/s\n";

                cpf = std::make_unique<ContourPF>(
                        dem, kpdb,
                        last_gps_lat, last_gps_lon,
                        last_gps_hdg, cfg.speed_mps,
                        cfg.cpf_n,
                        /*pos_sigma_m=*/   10.0,
                        /*hdg_sigma_deg=*/ eff_hdg_sigma,
                        /*rm=*/            &rm,
                        /*yaw_rate_dps=*/  yaw_rate,
                        /*yaw_sigma_dps=*/ eff_yaw_sigma);
                cpf_ready = true;
            }
            auto est = cpf->step(fix.radio_alt_m, cfg.baro_alt_m, dt_s,
                                 active_hdg_noise, 1.0, active_cpf_sigma);

            if (active_cpf_sigma > cfg.cpf_sigma) {
                // Flat terrain: Neff-weighted PF/DR blend + step-consistency filter.
                double neff_ratio = std::min(1.0, est.neff / cfg.cpf_n);

                if (trajectory.empty() || !trajectory.back().gps_denied) {
                    // Seed DR from last GPS fix
                    dr_lat   = last_gps_lat;
                    dr_lon   = last_gps_lon;
                    dr_hdg   = last_gps_hdg;
                    dr_speed = (est.speed_mps > 1.0) ? est.speed_mps : cfg.speed_mps;
                }

                // Advance DR heading: base from GPS yaw_rate + weak PF correction (α=0.1).
                // Pure GPS anchor prevents PF from dragging the course off; small PF
                // blend lets the DR gradually follow real curvature and keeps the step
                // filter from rejecting valid corrections near the end of the zone.
                {
                    double base_hdg = dr_hdg + dr_yaw_rate * dt_s;
                    double pf_diff  = est.heading_deg - base_hdg;
                    while (pf_diff >  180.0) pf_diff -= 360.0;
                    while (pf_diff < -180.0) pf_diff += 360.0;
                    constexpr double PF_HDG_ALPHA = 0.3;
                    dr_hdg = base_hdg + PF_HDG_ALPHA * neff_ratio * pf_diff;
                }
                {
                    double hdg_rad = dr_hdg * M_PI / 180.0;
                    double cos_l   = std::cos(dr_lat * M_PI / 180.0);
                    dr_lat += dr_speed * dt_s * std::cos(hdg_rad) / 111320.0;
                    dr_lon += dr_speed * dt_s * std::sin(hdg_rad) / (111320.0 * cos_l);
                }

                // Candidate position from Neff-weighted blend
                double cand_lat = neff_ratio * est.lat + (1.0 - neff_ratio) * dr_lat;
                double cand_lon = neff_ratio * est.lon + (1.0 - neff_ratio) * dr_lon;

                // Step-consistency filter: reject candidate if heading or distance
                // implied by prev→cand deviates too far from the DR heading.
                bool step_ok = true;
                if (!trajectory.empty() && trajectory.back().gps_denied) {
                    const auto& prev = trajectory.back();
                    double cos_l  = std::cos(prev.lat * M_PI / 180.0);
                    double dlat   = (cand_lat - prev.lat) * 111320.0;
                    double dlon   = (cand_lon - prev.lon) * 111320.0 * cos_l;
                    double step_m = std::sqrt(dlat * dlat + dlon * dlon);

                    double step_hdg = std::atan2(dlon, dlat) * 180.0 / M_PI;
                    if (step_hdg < 0) step_hdg += 360.0;
                    double dhdg = step_hdg - dr_hdg;
                    while (dhdg >  180.0) dhdg -= 360.0;
                    while (dhdg < -180.0) dhdg += 360.0;

                    double max_step_m      = dr_speed * dt_s * 2.5;
                    constexpr double MAX_HDG_JUMP = 25.0;  // °/step

                    if (step_m > max_step_m || std::abs(dhdg) > MAX_HDG_JUMP) {
                        step_ok = false;
                    }
                }

                tp.lat = step_ok ? cand_lat : dr_lat;
                tp.lon = step_ok ? cand_lon : dr_lon;

                if (neff_ratio > 0.4 && step_ok) {
                    dr_speed = est.speed_mps;
                }
            } else {
                tp.lat = est.lat;
                tp.lon = est.lon;
            }
            tp.heading_deg = est.heading_deg;
            tp.speed_mps   = est.speed_mps;
            tp.ncc         = est.neff / cfg.cpf_n;
        }
        trajectory.push_back(tp);

        {
            double cos_lat = std::cos(cfg.start_lat * M_PI / 180.0);
            double x_m = (tp.lon - cfg.start_lon) * 111320.0 * cos_lat;
            double y_m = (tp.lat - cfg.start_lat) * 111320.0;
            csv_out << std::fixed << std::setprecision(6)
                    << fix.timestamp_s << ","
                    << tp.lat << "," << tp.lon << ","
                    << std::setprecision(2)
                    << x_m << "," << y_m << ","
                    << tp.heading_deg << "," << tp.speed_mps << ","
                    << static_cast<int>(tp.ncc * cfg.cpf_n) << "\n";
        }
    }
    csv_out.close();

    // Средняя скорость и среднее направление (start→end)
    double mean_speed = 0.0; int n_denied = 0;
    for (auto& pt : trajectory) {
        if (pt.gps_denied) { mean_speed += pt.speed_mps; n_denied++; }
    }
    if (n_denied > 0) mean_speed /= n_denied;

    auto& last = trajectory.back();
    double cos_lat = std::cos(cfg.start_lat * M_PI / 180.0);
    double end_x_m = (last.lon - cfg.start_lon) * 111320.0 * cos_lat;
    double end_y_m = (last.lat - cfg.start_lat) * 111320.0;
    double dist_m  = std::sqrt(end_x_m * end_x_m + end_y_m * end_y_m);
    double mean_hdg = std::atan2(end_x_m, end_y_m) * 180.0 / M_PI;
    if (mean_hdg < 0) mean_hdg += 360.0;

    if (!trajectory.empty()) {
        Visualizer::save_trajectory_on_dem(dem, trajectory,
                                           cfg.start_lat, cfg.start_lon,
                                           cfg.out_dir + "/trajectory.png",
                                           cfg.jammer_zone,
                                           mean_hdg, 0.0,
                                           n_denied > 0 ? mean_speed : cfg.speed_mps);
    }

    std::cerr << "\n=== РЕЗУЛЬТАТ (COMBINED) ===\n"
              << "Позиция:   " << last.lat << "° N, " << last.lon << "° E\n"
              << "Вектор:    " << last.speed_mps << " м/с, "
              << last.heading_deg << "° (азимут)\n"
              << "Neff/N:    " << last.ncc << "\n"
              << "CSV:       " << cfg.out_dir << "/trn_estimate.csv\n"
              << "Траектория:" << cfg.out_dir << "/trajectory.png\n"
              << "\n=== ИТОГ ===\n"
              << "Старт:    lat=" << cfg.start_lat << " lon=" << cfg.start_lon
              << "  (x=0 м, y=0 м)\n"
              << "Финиш:    lat=" << last.lat << " lon=" << last.lon
              << std::fixed << std::setprecision(1)
              << "  (x=" << end_x_m << " м, y=" << end_y_m << " м)\n"
              << "Пройдено:  " << dist_m << " м\n"
              << "Скорость:  " << (n_denied > 0 ? mean_speed : cfg.speed_mps)
              << " м/с (оценка PF)\n"
              << "Азимут:    " << mean_hdg << "°\n";
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

    if (cfg.baro_alt_m < 0.0) {
        float dem_h = DemLoader::sample(dem, cfg.start_lat, cfg.start_lon, -9999.f);
        if (!fixes.empty() && dem_h > -9000.f) {
            cfg.baro_alt_m = dem_h + fixes.front().radio_alt_m;
            std::cerr << "[MAIN] baro авто: DEM(" << dem_h << ") + AGL("
                      << fixes.front().radio_alt_m << ") = " << cfg.baro_alt_m << " м\n";
        } else {
            std::cerr << "ERROR: --baro не задан и DEM в стартовой точке пустой\n";
            return 1;
        }
    }

    if (cfg.use_lnav)
        return run_line_navigator(cfg, dem, fixes, dt_s);
    else if (cfg.use_combined) {
        // В режиме чекпоинта (--source) сначала строим тепловую карту через TERCOM,
        // затем Combined перезаписывает trajectory.png и CSV точными оценками PF.
        if (!cfg.source_dir.empty()) {
            std::cerr << "[MAIN] Генерация тепловой карты корреляций (TERCOM)...\n";
            run_tercom(cfg, dem, fixes, dt_s);
            std::cerr << "[MAIN] correlation_heatmap.png сохранён. Запуск Combined...\n";
        }
        return run_combined(cfg, dem, fixes, dt_s);
    }
    else if (cfg.use_cpf)
        return run_contour_pf(cfg, dem, fixes, dt_s);
    else if (cfg.use_pf)
        return run_particle_filter(cfg, dem, fixes, dt_s);
    else
        return run_tercom(cfg, dem, fixes, dt_s);
}
