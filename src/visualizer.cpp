#include "visualizer.hpp"
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <cmath>
#include <iostream>
#include <algorithm>

// ── Вспомогательные форматировщики ──────────────────────────────────────────
static std::string deg_sym() { return "deg"; }

static std::string format_km(int off_idx, float step_m) {
    if (step_m <= 0.0f) return std::to_string(off_idx);
    float d_km = off_idx * step_m / 1000.0f;
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.1f km", d_km);
    return buf;
}

static std::string format_km_label(float dist_m) {
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%.1f", dist_m / 1000.0f);
    return buf;
}

static std::string format_ncc(float ncc) {
    char buf[16];
    std::snprintf(buf, sizeof(buf), "%.3f", ncc);
    return buf;
}

void Visualizer::save_correlation_heatmap(const cv::Mat& corr_map,
                                           const std::string& path,
                                           int best_az, int best_off,
                                           float offset_step_m, float best_ncc)
{
    if (corr_map.empty()) return;

    // corr_map: rows=azimuths, cols=offsets.
    // Транспонируем → result rows=offsets, cols=azimuths.
    // Итог: X-ось = азимут (0..360°), Y-ось = расстояние (0..N км).
    cv::Mat tr;
    cv::transpose(corr_map, tr);  // tr: rows=offsets, cols=azimuths

    // Нормализовать [-1,1] → [0,255] → JET
    cv::Mat norm;
    cv::normalize(tr, norm, 0, 255, cv::NORM_MINMAX, CV_8U);
    cv::Mat colored;
    cv::applyColorMap(norm, colored, cv::COLORMAP_JET);

    // Целевой размер: ~1200×600
    const int TARGET_W = 1200, TARGET_H = 600;
    const int MARGIN_B = 40, MARGIN_T = 30, MARGIN_L = 65, MARGIN_R = 20;
    int plot_w = TARGET_W - MARGIN_L - MARGIN_R;
    int plot_h = TARGET_H - MARGIN_T - MARGIN_B;
    cv::resize(colored, colored, cv::Size(plot_w, plot_h), 0, 0, cv::INTER_LINEAR);

    // Холст с полями
    cv::Mat canvas(TARGET_H, TARGET_W, CV_8UC3, cv::Scalar(30, 30, 30));
    colored.copyTo(canvas(cv::Rect(MARGIN_L, MARGIN_T, plot_w, plot_h)));

    // ── X-ось: азимут 0°–360° ──────────────────────────────────────────────
    int n_az_ticks = corr_map.rows;  // = n_azimuths (до транспонирования)
    for (int az : {0, 90, 180, 270, 360}) {
        if (n_az_ticks <= 0) break;
        int safe_az = std::min(az, n_az_ticks - 1);
        int px = MARGIN_L + static_cast<int>(static_cast<float>(safe_az) / n_az_ticks * plot_w);
        int py_bottom = MARGIN_T + plot_h;
        cv::line(canvas, {px, py_bottom}, {px, py_bottom + 8},
                 cv::Scalar(220,220,220), 2);
        // grid line
        cv::line(canvas, {px, MARGIN_T}, {px, py_bottom},
                 cv::Scalar(60,60,60), 1);
        std::string lbl = std::to_string(az);
        int lw = (az == 360) ? 16 : 12;
        cv::putText(canvas, lbl, {px - lw, py_bottom + 22},
                    cv::FONT_HERSHEY_SIMPLEX, 0.55, cv::Scalar(230,230,230), 1);
    }
    cv::putText(canvas, "Azimuth",
                {MARGIN_L + plot_w/2 - 35, TARGET_H - 4},
                cv::FONT_HERSHEY_SIMPLEX, 0.55, cv::Scalar(220,220,220), 1);

    // ── Y-ось: расстояние ──────────────────────────────────────────────────
    int n_off_rows = corr_map.cols;  // = n_offsets (до транспонирования)
    if (n_off_rows > 0 && offset_step_m > 0.0f) {
        float total_dist_m = n_off_rows * offset_step_m;
        float step_km = (total_dist_m < 5000.0f) ? 0.5f : 1.0f;
        float step_m  = step_km * 1000.0f;
        for (float d = 0.0f; d <= total_dist_m + 1.0f; d += step_m) {
            int py = MARGIN_T + static_cast<int>(d / total_dist_m * plot_h);
            if (py > MARGIN_T + plot_h) break;
            cv::line(canvas, {MARGIN_L - 8, py}, {MARGIN_L, py},
                     cv::Scalar(220,220,220), 2);
            // grid line
            cv::line(canvas, {MARGIN_L, py}, {MARGIN_L + plot_w, py},
                     cv::Scalar(60,60,60), 1);
            std::string lbl = format_km_label(d) + " km";
            cv::putText(canvas, lbl, {2, py + 5},
                        cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(230,230,230), 1);
        }
    }

    // ── Отметить найденный максимум ─────────────────────────────────────────
    if (best_az >= 0 && best_off >= 0) {
        int n_az  = tr.cols;
        int n_off = tr.rows;
        float px_f = MARGIN_L + static_cast<float>(best_az)  / n_az  * plot_w;
        float py_f = MARGIN_T + static_cast<float>(best_off) / n_off * plot_h;
        cv::Point marker{static_cast<int>(px_f), static_cast<int>(py_f)};

        // Белый крест с чёрной обводкой
        cv::drawMarker(canvas, marker, cv::Scalar(0,0,0),
                       cv::MARKER_CROSS, 30, 4, cv::LINE_AA);
        cv::drawMarker(canvas, marker, cv::Scalar(255,255,255),
                       cv::MARKER_CROSS, 28, 2, cv::LINE_AA);

        // Подпись с непрозрачным фоном
        std::string label = "az=" + std::to_string(best_az) + deg_sym()
                          + "   NCC=" + format_ncc(best_ncc);
        int tx = marker.x + 14, ty = marker.y - 10;
        if (tx + 260 > TARGET_W) tx = marker.x - 270;
        if (ty < MARGIN_T + 18)  ty = marker.y + 22;
        // фон-прямоугольник
        int baseline = 0;
        cv::Size ts = cv::getTextSize(label, cv::FONT_HERSHEY_SIMPLEX, 0.65, 2, &baseline);
        cv::rectangle(canvas, {tx - 4, ty - ts.height - 4}, {tx + ts.width + 4, ty + baseline + 4},
                      cv::Scalar(10,10,10), cv::FILLED);
        cv::putText(canvas, label, {tx, ty}, cv::FONT_HERSHEY_SIMPLEX,
                    0.65, cv::Scalar(255,255,100), 2, cv::LINE_AA);
    }

    cv::imwrite(path, canvas);
    std::cerr << "[VIS] Тепловая карта сохранена: " << path << "\n";
}

void Visualizer::save_trajectory_on_dem(const DemData& dem,
                                         const std::vector<TrajectoryPoint>& traj,
                                         double start_lat, double start_lon,
                                         const std::string& path,
                                         const std::vector<double>& jammer_zone,
                                         double found_az_deg,
                                         double found_dist_m,
                                         double found_speed_mps)
{
    if (dem.elev.empty() || traj.empty()) return;

    int W = dem.elev.cols, H = dem.elev.rows;

    // --- Цветная карта высот ---
    double emin, emax;
    cv::minMaxLoc(dem.elev, &emin, &emax);

    // LUT: 256 оттенков green→yellow→brown (3 отдельных канала BGR)
    cv::Mat lut_b(1, 256, CV_8U), lut_g(1, 256, CV_8U), lut_r(1, 256, CV_8U);
    for (int i = 0; i < 256; ++i) {
        double t = i / 255.0;
        uchar r, g, b;
        if (t < 0.35) {
            double s = t / 0.35;
            r = static_cast<uchar>(60  + s * 80);
            g = static_cast<uchar>(130 + s * 50);
            b = static_cast<uchar>(50  - s * 20);
        } else if (t < 0.65) {
            double s = (t - 0.35) / 0.30;
            r = static_cast<uchar>(140 + s * 80);
            g = static_cast<uchar>(180 - s * 30);
            b = static_cast<uchar>(30  + s * 20);
        } else {
            double s = (t - 0.65) / 0.35;
            r = static_cast<uchar>(220 + s * 35);
            g = static_cast<uchar>(150 - s * 30);
            b = static_cast<uchar>(50  + s * 70);
        }
        lut_b.at<uchar>(0, i) = b;
        lut_g.at<uchar>(0, i) = g;
        lut_r.at<uchar>(0, i) = r;
    }

    cv::Mat elev_norm;
    cv::normalize(dem.elev, elev_norm, 0, 255, cv::NORM_MINMAX, CV_8U);
    cv::Mat ch_b, ch_g, ch_r;
    cv::LUT(elev_norm, lut_b, ch_b);
    cv::LUT(elev_norm, lut_g, ch_g);
    cv::LUT(elev_norm, lut_r, ch_r);
    cv::Mat vis;
    cv::merge(std::vector<cv::Mat>{ch_b, ch_g, ch_r}, vis);

    // --- Hillshade ---
    cv::Mat dx_m, dy_m;
    cv::Sobel(dem.elev, dx_m, CV_32F, 1, 0, 3);
    cv::Sobel(dem.elev, dy_m, CV_32F, 0, 1, 3);
    float mpp = static_cast<float>(std::abs(dem.gt.pixel_h) * 111320.0);
    dx_m /= (8.0f * mpp);
    dy_m /= (8.0f * mpp);

    constexpr float SUN_AZ  = static_cast<float>(315.0 * M_PI / 180.0);
    constexpr float SUN_ALT = static_cast<float>(45.0  * M_PI / 180.0);
    cv::Mat shade(H, W, CV_32F);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float dxv = dx_m.at<float>(y, x);
            float dyv = dy_m.at<float>(y, x);
            float slope  = std::atan(std::sqrt(dxv*dxv + dyv*dyv));
            float aspect = std::atan2(-dxv, dyv);
            float s = std::sin(SUN_ALT) * std::cos(slope)
                    + std::cos(SUN_ALT) * std::sin(slope) * std::cos(SUN_AZ - aspect);
            shade.at<float>(y, x) = std::max(0.0f, s);
        }
    }

    // Наложить тень: пиксели *= (0.4 + 0.6 * shade)
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float s = 0.4f + 0.6f * shade.at<float>(y, x);
            cv::Vec3b& px = vis.at<cv::Vec3b>(y, x);
            px[0] = cv::saturate_cast<uchar>(px[0] * s);
            px[1] = cv::saturate_cast<uchar>(px[1] * s);
            px[2] = cv::saturate_cast<uchar>(px[2] * s);
        }
    }

    // --- Зона глушения GPS (полупрозрачный красный прямоугольник) ---
    if (jammer_zone.size() == 4) {
        auto [jx1, jy1] = DemLoader::geo_to_pixel(dem.gt, jammer_zone[0], jammer_zone[1]);
        auto [jx2, jy2] = DemLoader::geo_to_pixel(dem.gt, jammer_zone[2], jammer_zone[3]);
        cv::Point jp1{static_cast<int>(std::round(std::min(jx1,jx2))),
                      static_cast<int>(std::round(std::min(jy1,jy2)))};
        cv::Point jp2{static_cast<int>(std::round(std::max(jx1,jx2))),
                      static_cast<int>(std::round(std::max(jy1,jy2)))};
        // Полупрозрачная заливка (blend)
        cv::Mat overlay = vis.clone();
        cv::rectangle(overlay, jp1, jp2, cv::Scalar(0, 0, 180), -1);
        cv::addWeighted(overlay, 0.30, vis, 0.70, 0, vis);
        // Рамка
        cv::rectangle(vis, jp1, jp2, cv::Scalar(0, 0, 255), 2);
        // Подпись
        cv::putText(vis, "GPS DENIED",
                    jp1 + cv::Point(4, 18),
                    cv::FONT_HERSHEY_SIMPLEX, 0.55, cv::Scalar(0,0,0), 3);
        cv::putText(vis, "GPS DENIED",
                    jp1 + cv::Point(4, 18),
                    cv::FONT_HERSHEY_SIMPLEX, 0.55, cv::Scalar(60,60,255), 1);
    }

    // --- Пиксельные координаты точек траектории ---
    std::vector<cv::Point> pts;
    for (const auto& p : traj) {
        auto [px, py] = DemLoader::geo_to_pixel(dem.gt, p.lat, p.lon);
        pts.push_back({static_cast<int>(std::round(px)),
                       static_cast<int>(std::round(py))});
    }

    // Пиксель старта (нужен раньше для линии накопления)
    auto [sx, sy] = DemLoader::geo_to_pixel(dem.gt, start_lat, start_lon);
    cv::Point sp{static_cast<int>(std::round(sx)), static_cast<int>(std::round(sy))};

    // Линия от старта до первой TRN-точки (период накопления профиля, ~40 шагов)
    if (!pts.empty()) {
        cv::line(vis, sp, pts[0], cv::Scalar(255, 255, 255), 5, cv::LINE_AA);
        cv::line(vis, sp, pts[0], cv::Scalar(0, 220, 255), 3, cv::LINE_AA);  // жёлтый = вне зоны
    }

    // Предвычислить bbox зоны глушения (в координатах) для точного окрашивания
    double jlat_min = 0, jlat_max = 0, jlon_min = 0, jlon_max = 0;
    bool has_jammer = (jammer_zone.size() == 4);
    if (has_jammer) {
        jlat_min = std::min(jammer_zone[0], jammer_zone[2]);
        jlat_max = std::max(jammer_zone[0], jammer_zone[2]);
        jlon_min = std::min(jammer_zone[1], jammer_zone[3]);
        jlon_max = std::max(jammer_zone[1], jammer_zone[3]);
    }

    // Линия маршрута: красная если оценённая позиция внутри зоны глушения
    for (size_t i = 1; i < pts.size(); ++i) {
        bool in_jammer = has_jammer &&
            traj[i].lat >= jlat_min && traj[i].lat <= jlat_max &&
            traj[i].lon >= jlon_min && traj[i].lon <= jlon_max;
        cv::Scalar col = in_jammer
                         ? cv::Scalar(30, 30, 220)   // красный внутри зоны
                         : cv::Scalar(0, 220, 255);  // жёлтый вне зоны
        cv::line(vis, pts[i-1], pts[i], cv::Scalar(255, 255, 255), 5, cv::LINE_AA);
        cv::line(vis, pts[i-1], pts[i], col, 3, cv::LINE_AA);
    }

    // Финиш — синий заливной квадрат с белой рамкой
    if (!pts.empty()) {
        cv::Point ep = pts.back();
        cv::rectangle(vis, ep - cv::Point(10,10), ep + cv::Point(10,10),
                      cv::Scalar(255,255,255), -1);
        cv::rectangle(vis, ep - cv::Point(8,8),  ep + cv::Point(8,8),
                      cv::Scalar(220, 80, 30), -1);
    }

    // Старт — зелёный заливной круг с белой рамкой (рисуем ПОСЛЕ линии, поверх)
    cv::circle(vis, sp, 11, cv::Scalar(255, 255, 255), -1);
    cv::circle(vis, sp, 9,  cv::Scalar(30, 210, 30),  -1);

    // Подписи
    cv::putText(vis, "Start", sp + cv::Point(13, 4),
                cv::FONT_HERSHEY_SIMPLEX, 0.55, cv::Scalar(0,0,0),   3);
    cv::putText(vis, "Start", sp + cv::Point(13, 4),
                cv::FONT_HERSHEY_SIMPLEX, 0.55, cv::Scalar(40,220,40), 1);
    if (!pts.empty()) {
        cv::putText(vis, "End", pts.back() + cv::Point(13, 4),
                    cv::FONT_HERSHEY_SIMPLEX, 0.55, cv::Scalar(0,0,0),   3);
        cv::putText(vis, "End", pts.back() + cv::Point(13, 4),
                    cv::FONT_HERSHEY_SIMPLEX, 0.55, cv::Scalar(220,100,40), 1);
    }

    // --- Стрелка найденного азимута (TERCOM) ---
    if (found_az_deg >= 0.0 && found_dist_m > 0.0) {
        double az_rad = found_az_deg * M_PI / 180.0;
        double cos_lat = std::cos(start_lat * M_PI / 180.0);
        double end_lat = start_lat + found_dist_m * std::cos(az_rad) / 111320.0;
        double end_lon = start_lon + found_dist_m * std::sin(az_rad) / (111320.0 * cos_lat);
        auto [ex, ey] = DemLoader::geo_to_pixel(dem.gt, end_lat, end_lon);
        cv::Point ep_arrow{static_cast<int>(std::round(ex)),
                           static_cast<int>(std::round(ey))};
        // Белая обводка + жёлтая стрелка
        cv::arrowedLine(vis, sp, ep_arrow, cv::Scalar(255,255,255), 6, cv::LINE_AA, 0, 0.05);
        cv::arrowedLine(vis, sp, ep_arrow, cv::Scalar(0, 220, 255), 3, cv::LINE_AA, 0, 0.05);
        // Подпись у острия
        std::string lbl;
        {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "az=%.0fdeg  v=%.0f m/s",
                          found_az_deg, found_speed_mps);
            lbl = buf;
        }
        cv::putText(vis, lbl, ep_arrow + cv::Point(10, -6),
                    cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0,0,0), 3, cv::LINE_AA);
        cv::putText(vis, lbl, ep_arrow + cv::Point(10, -6),
                    cv::FONT_HERSHEY_SIMPLEX, 0.5, cv::Scalar(0, 255, 220), 1, cv::LINE_AA);
    }

    // --- Обрезать до области маршрута + margin ---
    if (!pts.empty()) {
        // Pixel bbox всех точек (включая старт)
        int xmin = sp.x, xmax = sp.x, ymin = sp.y, ymax = sp.y;
        for (const auto& p : pts) {
            xmin = std::min(xmin, p.x); xmax = std::max(xmax, p.x);
            ymin = std::min(ymin, p.y); ymax = std::max(ymax, p.y);
        }
        // Также включить зону глушения в bbox
        if (jammer_zone.size() == 4) {
            auto [jx1, jy1] = DemLoader::geo_to_pixel(dem.gt, jammer_zone[0], jammer_zone[1]);
            auto [jx2, jy2] = DemLoader::geo_to_pixel(dem.gt, jammer_zone[2], jammer_zone[3]);
            xmin = std::min(xmin, (int)std::round(std::min(jx1,jx2)));
            xmax = std::max(xmax, (int)std::round(std::max(jx1,jx2)));
            ymin = std::min(ymin, (int)std::round(std::min(jy1,jy2)));
            ymax = std::max(ymax, (int)std::round(std::max(jy1,jy2)));
        }
        // Добавить margin ~4 км в пикселях (+30% к прежнему 3 км)
        double mpp_y = DemLoader::meters_per_pixel_y(dem.gt);
        int margin_px = static_cast<int>(4000.0 / mpp_y);
        int cx1 = std::max(0, xmin - margin_px);
        int cy1 = std::max(0, ymin - margin_px);
        int cx2 = std::min(W - 1, xmax + margin_px);
        int cy2 = std::min(H - 1, ymax + margin_px);
        cv::Rect roi(cx1, cy1, cx2 - cx1, cy2 - cy1);
        if (roi.width > 0 && roi.height > 0)
            vis = vis(roi).clone();
    }

    // Ограничить размер выходного изображения
    constexpr int MAX_DIM = 2048;
    if (vis.rows > MAX_DIM || vis.cols > MAX_DIM) {
        double scale = static_cast<double>(MAX_DIM) / std::max(vis.rows, vis.cols);
        cv::resize(vis, vis, cv::Size(), scale, scale, cv::INTER_AREA);
    }

    cv::imwrite(path, vis);
    std::cerr << "[VIS] Траектория сохранена: " << path << "\n";
}

void Visualizer::show(const std::string& title, const cv::Mat& img, int wait_ms) {
    try {
        cv::imshow(title, img);
        cv::waitKey(wait_ms);
    } catch (...) {
        // Нет дисплея — не падаем
    }
}
