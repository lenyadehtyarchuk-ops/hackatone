#pragma once
#include "dem_loader.hpp"
#include <vector>
#include <deque>
#include <random>

struct Particle {
    double lat, lon, heading_deg, speed_mps, weight;
    std::deque<double> hdg_history;  // [0]=курс на последнем шаге, [1]=предпоследнем, ...
};

struct PFEstimate {
    double lat, lon, heading_deg, speed_mps;
    double neff;   // эффективное число частиц
};

// Particle Filter для TRN с профильным обновлением.
// Каждая частица хранит текущий курс. При обновлении веса вычисляются
// сравнением накопленного профиля AGL с ЦМР вдоль пути частицы назад.
class ParticleFilter {
public:
    // profile_win — длина профиля в шагах (как min-profile в TERCOM).
    // Частицы получают нулевой вес вклад до накопления profile_win измерений.
    ParticleFilter(const DemData& dem,
                   double start_lat,
                   double start_lon,
                   double heading_deg,
                   double speed_mps,
                   int    n_particles    = 1000,
                   int    profile_win    = 40,
                   double pos_sigma_m   = 100.0,
                   double hdg_sigma_deg = 30.0);

    // Один шаг: predict → добавить измерение → profile update если буфер полон → estimate.
    PFEstimate step(double agl_m,
                    double baro_alt_m,
                    double dt_s,
                    double hdg_noise_deg = 3.0,
                    double spd_noise_mps = 1.0,
                    double meas_sigma_m  = 3.0);

    PFEstimate estimate() const;
    int n() const { return n_; }

private:
    std::vector<Particle> p_;
    const DemData&        dem_;
    std::mt19937          rng_;
    int                   n_;
    int                   win_;          // размер профильного окна
    std::deque<double>    agl_buf_;      // накопленные измерения AGL
    double                speed_nom_;    // номинальная скорость для реконструкции пути

    void predict(double dt_s, double hdg_noise, double spd_noise);
    void profile_update(double baro_alt_m, double dt_s, double meas_sigma_m);
    void resample();
};
