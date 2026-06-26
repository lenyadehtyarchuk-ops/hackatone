#pragma once
#include <opencv2/video/tracking.hpp>

// Фильтр Калмана для навигационного состояния [x, y, vx, vy] в метрах (локальная СК).
// Модель: constant velocity. Измерения: позиция (x, y) от корреляционного поиска.
class TrnKalmanFilter {
public:
    // dt_s — номинальный шаг по времени (с)
    explicit TrnKalmanFilter(double dt_s = 1.0);

    // Инициализировать из первого измерения позиции
    void init(double x_m, double y_m, double vx_mps = 0.0, double vy_mps = 0.0);

    // Предсказание на dt_s вперёд
    void predict(double dt_s = -1.0);

    // Коррекция по наблюдению (x, y) из корреляционного поиска
    // measurement_noise_m — СКО позиции (можно задать по уровню NCC)
    void correct(double x_m, double y_m, double measurement_noise_m = 50.0);

    // Текущая оценка состояния
    double x()  const;
    double y()  const;
    double vx() const;
    double vy() const;
    double speed_mps() const;
    double heading_deg() const;

    bool is_initialized() const { return m_initialized; }

private:
    cv::KalmanFilter m_kf;
    double m_dt;
    bool m_initialized{false};
};
