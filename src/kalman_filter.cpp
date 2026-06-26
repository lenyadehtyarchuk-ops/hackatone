#include "kalman_filter.hpp"
#include <cmath>

TrnKalmanFilter::TrnKalmanFilter(double dt_s)
    : m_kf(4, 2, 0, CV_64F), m_dt(dt_s)
{
    // Матрица перехода F (constant velocity): x' = x + vx*dt
    //  | 1  0  dt  0 |
    //  | 0  1  0  dt |
    //  | 0  0  1   0 |
    //  | 0  0  0   1 |
    m_kf.transitionMatrix = (cv::Mat_<double>(4,4)
        << 1, 0, m_dt, 0,
           0, 1, 0,    m_dt,
           0, 0, 1,    0,
           0, 0, 0,    1);

    // Матрица измерений H: наблюдаем только x, y
    m_kf.measurementMatrix = (cv::Mat_<double>(2,4)
        << 1, 0, 0, 0,
           0, 1, 0, 0);

    // Шум процесса Q — малый (полёт относительно гладкий)
    double q_pos = 0.1;   // м²/с²
    double q_vel = 1.0;   // (м/с)²/с
    m_kf.processNoiseCov = (cv::Mat_<double>(4,4)
        << q_pos, 0,     0,     0,
           0,     q_pos, 0,     0,
           0,     0,     q_vel, 0,
           0,     0,     0,     q_vel);

    // Ковариация ошибки начального состояния
    cv::setIdentity(m_kf.errorCovPost, cv::Scalar(1000));
}

void TrnKalmanFilter::init(double x_m, double y_m, double vx_mps, double vy_mps) {
    m_kf.statePost = (cv::Mat_<double>(4,1) << x_m, y_m, vx_mps, vy_mps);
    cv::setIdentity(m_kf.errorCovPost, cv::Scalar(500));
    m_initialized = true;
}

void TrnKalmanFilter::predict(double dt_s) {
    if (!m_initialized) return;
    double dt = (dt_s > 0) ? dt_s : m_dt;
    m_kf.transitionMatrix.at<double>(0, 2) = dt;
    m_kf.transitionMatrix.at<double>(1, 3) = dt;
    m_kf.predict();
}

void TrnKalmanFilter::correct(double x_m, double y_m, double measurement_noise_m) {
    if (!m_initialized) {
        init(x_m, y_m);
        return;
    }
    double r = measurement_noise_m * measurement_noise_m;
    m_kf.measurementNoiseCov = (cv::Mat_<double>(2,2)
        << r, 0,
           0, r);
    cv::Mat meas = (cv::Mat_<double>(2,1) << x_m, y_m);
    m_kf.correct(meas);
}

double TrnKalmanFilter::x()  const { return m_kf.statePost.at<double>(0); }
double TrnKalmanFilter::y()  const { return m_kf.statePost.at<double>(1); }
double TrnKalmanFilter::vx() const { return m_kf.statePost.at<double>(2); }
double TrnKalmanFilter::vy() const { return m_kf.statePost.at<double>(3); }

double TrnKalmanFilter::speed_mps() const {
    double vxv = vx(), vyv = vy();
    return std::sqrt(vxv*vxv + vyv*vyv);
}

double TrnKalmanFilter::heading_deg() const {
    // Азимут: север = 0°, восток = 90°
    double az = std::atan2(vx(), vy()) * 180.0 / M_PI;
    return az < 0 ? az + 360.0 : az;
}
