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
    double baro_alt_m;         // барометрическая высота MSL
    double speed_mps;          // предполагаемая скорость (м/с) для шага сэмплирования
    double sample_dt_s;        // период дискретизации (с)
    int    min_profile_len;    // минимальный размер профиля для запуска поиска
    double search_radius_m;    // максимальный радиус поиска (м)
    int    azimuth_step_deg;   // шаг по азимуту (1 — полный перебор)
};

class TerrainCorrelator {
public:
    explicit TerrainCorrelator(const DemData& dem, const CorrelatorConfig& cfg);

    // Добавить измерение радиовысотомера; если профиль достаточно длинный — выполнить поиск
    // Возвращает nullopt если профиль ещё слишком короткий
    std::optional<CorrelationResult> add_measurement(double radio_alt_m);

    // Принудительно выполнить поиск по текущему накопленному профилю
    std::optional<CorrelationResult> run_search();

    int profile_size() const { return static_cast<int>(m_measured.size()); }
    const std::vector<double>& profile() const { return m_measured; }

private:
    // Извлечь эталонный профиль из ЦМР вдоль луча в направлении azimuth_deg
    // длиной n_points с шагом step_px (пиксели)
    std::vector<float> extract_profile(double start_lat, double start_lon,
                                       double azimuth_deg, int n_points,
                                       double step_px_x, double step_px_y) const;

    // Взвешенная нормированная кросс-корреляция (веса по отклонению M от среднего)
    static std::vector<float> ncc_sliding_weighted(
        const std::vector<double>& measured,
        const std::vector<float>&  reference,
        const std::vector<double>& weights);

    const DemData* m_dem;
    CorrelatorConfig m_cfg;
    std::vector<double> m_measured;  // накопленный измеренный профиль высот рельефа
};
