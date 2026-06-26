#include "visualizer.hpp"
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <cmath>
#include <iostream>
#include <algorithm>

void Visualizer::save_correlation_heatmap(const cv::Mat& corr_map,
                                           const std::string& path,
                                           int best_az, int best_off)
{
    if (corr_map.empty()) return;

    // Нормализовать [-1, 1] → [0, 255]
    cv::Mat norm;
    cv::normalize(corr_map, norm, 0, 255, cv::NORM_MINMAX, CV_8U);

    cv::Mat colored;
    cv::applyColorMap(norm, colored, cv::COLORMAP_JET);

    // Масштабировать для читабельности
    int scale_x = std::max(1, 1200 / corr_map.cols);
    int scale_y = std::max(1, 720  / corr_map.rows);
    int scale = std::min(scale_x, scale_y);
    if (scale > 1)
        cv::resize(colored, colored, cv::Size(), scale, scale, cv::INTER_NEAREST);

    // Отметить найденный максимум
    if (best_az >= 0 && best_off >= 0) {
        int px = best_off * scale + scale / 2;
        int py = best_az  * scale + scale / 2;
        cv::drawMarker(colored, {px, py}, cv::Scalar(255, 255, 255),
                       cv::MARKER_CROSS, 15, 2);
    }

    // Подписи осей
    cv::putText(colored, "Offset (m) ->",
                {10, colored.rows - 10}, cv::FONT_HERSHEY_SIMPLEX,
                0.5, cv::Scalar(255,255,255), 1);
    cv::putText(colored, "Azimuth (deg)",
                {10, 20}, cv::FONT_HERSHEY_SIMPLEX,
                0.5, cv::Scalar(255,255,255), 1);

    cv::imwrite(path, colored);
    std::cerr << "[VIS] Тепловая карта сохранена: " << path << "\n";
}

void Visualizer::save_trajectory_on_dem(const DemData& dem,
                                         const std::vector<TrajectoryPoint>& traj,
                                         double start_lat, double start_lon,
                                         const std::string& path)
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

    // --- Пиксельные координаты точек траектории ---
    std::vector<cv::Point> pts;
    for (const auto& p : traj) {
        auto [px, py] = DemLoader::geo_to_pixel(dem.gt, p.lat, p.lon);
        pts.push_back({static_cast<int>(std::round(px)),
                       static_cast<int>(std::round(py))});
    }

    // Линия маршрута: белая обводка + красная
    for (size_t i = 1; i < pts.size(); ++i) {
        cv::line(vis, pts[i-1], pts[i], cv::Scalar(255, 255, 255), 5, cv::LINE_AA);
        cv::line(vis, pts[i-1], pts[i], cv::Scalar(30, 30, 220),   3, cv::LINE_AA);
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
    auto [sx, sy] = DemLoader::geo_to_pixel(dem.gt, start_lat, start_lon);
    cv::Point sp{static_cast<int>(std::round(sx)), static_cast<int>(std::round(sy))};
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
