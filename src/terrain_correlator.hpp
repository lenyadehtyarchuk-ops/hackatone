#pragma once
#include "dem_loader.hpp"
#include <vector>
#include <optional>
#include <opencv2/core.hpp>

struct CorrelationResult {
    int    best_azimuth_deg;   // 0..359
    double best_offset_m;      // расстояние от центра поиска до НАЧАЛА профиля (м)
    double best_ncc;           // максимальный коэффициент корреляции [-1, 1]
    double start_lat;
    double start_lon;
    double found_lat;          // позиция НАЧАЛА профиля (точка d*)
    double found_lon;
    double current_lat;        // ТЕКУЩАЯ позиция ВС = found + (n-1)*шаг
    double current_lon;
    double ground_speed_mps;   // путевая скорость = step_m / dt_s
    int    profile_len;        // длина использованного профиля
    cv::Mat corr_map;          // float32, 360 строк × max_offset_px столбцов
    bool   low_confidence;     // true если NCC < порога или рельеф слишком плоский
    double profile_std;        // σ измеренного профиля (информативность рельефа)
};

struct CorrelatorConfig {
    double start_lat;
    double start_lon;
    double baro_alt_m;          // барометрическая высота MSL
    double speed_mps;           // предполагаемая скорость (м/с) для шага сэмплирования
    double sample_dt_s;         // период дискретизации (с)
    int    min_profile_len;     // минимальный размер профиля для запуска поиска
    double search_radius_m;     // максимальный радиус поиска (м)
    int    azimuth_step_deg;    // шаг по азимуту (1 — полный перебор)
    double gradient_weight;     // вес градиентного канала [0..1], 0 = только высоты
    // Секторный поиск по azimuth (для непрямой траектории)
    double prior_heading_deg;   // текущий heading от Калмана (-1 = отключено, полный 360°)
    double heading_search_deg;  // полуширина сектора поиска (°), активна если prior_heading>=0
};

class TerrainCorrelator {
public:
    explicit TerrainCorrelator(const DemData& dem, const CorrelatorConfig& cfg);

    // Добавить измерение радиовысотомера; если профиль достаточно длинный — выполнить поиск
    std::optional<CorrelationResult> add_measurement(double radio_alt_m);

    // Принудительно выполнить поиск по текущему накопленному профилю
    std::optional<CorrelationResult> run_search();

    int profile_size() const { return static_cast<int>(m_measured.size()); }
    const std::vector<double>& profile() const { return m_measured; }

private:
    const DemData*   m_dem;
    CorrelatorConfig m_cfg;
    std::vector<double> m_measured;  // профиль высот рельефа (h = baro - AGL)
    std::vector<double> m_grad;      // градиент профиля: m_measured[k] - m_measured[k-1]
};
