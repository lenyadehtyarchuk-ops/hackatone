# Terrain-Referenced Navigation (TRN)

## Задача
Определить координаты и вектор скорости БПЛА без ГНСС, используя радиовысотомер и цифровую модель рельефа (ЦМР).

## Физика
- Радиовысотомер → высота над рельефом (AGL)
- Барометрический высотомер → абсолютная высота (MSL)
- Высота рельефа = MSL − AGL

---

## Режимы работы (--tercom / --combined)

### TERCOM (базовый, соответствует ТЗ)
```
./build/trn --dem ... --nmea ... --baro ... --lat ... --lon ... --speed ...
```
1. Накапливает профиль измерений `M[n]` из NMEA-потока
2. Для каждого θ ∈ [0°, 360°) строит эталонный профиль из ЦМР
3. Двухканальный NCC: `0.6·NCC(высоты) + 0.4·NCC(градиент)` — устойчивость к барометрическому дрейфу
4. θ*, d* = argmax NCC → азимут и пройденное расстояние
5. Фильтр Калмана [x, y, vx, vy] сглаживает решения, физфильтр отклоняет скачки > 250 м/с
6. В GPS-зоне траектория строится из реальных координат NMEA, в зоне глушения — из TERCOM-оценок

**Артефакты:** `correlation_heatmap.png` (тепловая карта 360°×d), `trajectory.png` (стрелка найденного азимута на ЦМР)

### Combined (расширенный, ноу-хау)
```
./build/trn --dem ... --nmea ... --baro ... ... --combined \
            --lnav-thresh 0.08 --pf-hdg-noise 1.5 --cpf-sigma 20
```
Два компонента:
- **ContourPF** — 3000 частиц по изолинии `DEM(lat,lon) = baro − AGL`; вес = `exp(−residual²/2σ²)`
- **RidgeMap** — матрица Гессе H=[[∂²h/∂x², ∂²h/∂xy],[∂²h/∂xy, ∂²h/∂y²]] → собственные числа λ₁≤λ₂ → хребты (λ₁<<0, λ₂≈0) / впадины (λ₂>>0, λ₁≈0); NMS → однопиксельные линии; при пересечении хребта бустит веса частиц (радиус 150 м)

**Время работы:** ~1.7 с на весь полёт (vs 3–4 мин у TERCOM на 1095 измерений)

---

## Стек
- C++17, OpenCV 4.10, GDAL, cv::KalmanFilter
- Python 3 + numpy + GDAL для генератора NMEA и визуализации

## Зависимости
```
sudo apt install -y libopencv-dev libgdal-dev libeigen3-dev python3-numpy python3-gdal cmake build-essential
```

## Сборка
```bash
mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc)
```

---

## Демонстрационные запуски (Эльбрус, реальный Copernicus GLO-30)

### TERCOM demo — прямолинейный полёт, чистая демонстрация алгоритма
```bash
./build/trn \
  --dem data/dem/elbrus.tif \
  --nmea data/nmea/tercom_straight.nmea \
  --baro 6500 --lat 43.503370 --lon 42.492140 \
  --speed 40 --min-profile 35 --radius 3000 \
  --jammer-zone 43.4639,42.5139,43.4925,42.5690 \
  --out results/tercom_demo
# Результат: az=120° (истинный 120°), NCC=0.924
```

### Combined — реалистичный маршрут (переменная скорость, изгибы курса)
```bash
./build/trn \
  --dem data/dem/elbrus.tif \
  --nmea data/nmea/elbrus_realistic.nmea \
  --baro 6500 --lat 43.503370 --lon 42.492140 --azimuth 120 \
  --speed 40 --out results/combined_realistic \
  --jammer-zone 43.4639,42.5139,43.4925,42.5690 \
  --combined --lnav-thresh 0.08 --pf-hdg-noise 1.5 --cpf-sigma 20
```

### Combined — синусоидальный маршрут
```bash
./build/trn \
  --dem data/dem/elbrus.tif \
  --nmea data/nmea/elbrus_sine.nmea \
  --baro 6500 --lat 43.503370 --lon 42.492140 --azimuth 120 \
  --speed 40 --out results/combined_sine \
  --jammer-zone 43.4639,42.5139,43.4925,42.5690 \
  --combined --lnav-thresh 0.08 --pf-hdg-noise 1.5 --cpf-sigma 20
```

### Генерация NMEA
```bash
# Прямолинейный (для TERCOM demo)
python3 scripts/generate_nmea.py \
  --dem data/dem/elbrus.tif \
  --lat 43.503370 --lon 42.492140 --azimuth 120 --speed 40 \
  --duration 180 --baro 6500 --noise 2 \
  --jammer-zone 43.4639,42.5139,43.4925,42.5690 \
  --out data/nmea/tercom_straight.nmea

# Реалистичный (переменная скорость, изгибы, dt=0.2 с)
python3 scripts/generate_nmea.py \
  --dem data/dem/elbrus.tif \
  --lat 43.503370 --lon 42.492140 --azimuth 120 \
  --speed 40 --speed-min 25 --speed-max 60 --heading-sigma 3 \
  --duration 219 --baro 6500 --noise 2 \
  --jammer-zone 43.4639,42.5139,43.4925,42.5690 \
  --out data/nmea/elbrus_realistic.nmea
```

---

## Результаты верификации

| Сценарий | Алгоритм | Азимут истинный | Найден | NCC | Время |
|---|---|---|---|---|---|
| Прямолинейный, σ=2 м | TERCOM | 120° | 120° | 0.924 | 2 с |
| Реалистичный (изгибы), σ=2 м | Combined | ~120° | 151° | — | 1.7 с |
| Синусоидальный, σ=2 м | Combined | ~120° | 147° | — | 1.7 с |
| Реалистичный, σ=2 м | TERCOM | 120° | 125° | 0.801 | 4 мин |

Вывод: Combined в 100× быстрее TERCOM и точнее на нелинейных траекториях.

---

## Структура проекта
```
src/
  main.cpp                    — CLI: --dem --nmea --baro --lat --lon --speed
                                      --min-profile --radius --az-step --sliding
                                      --combined --lnav-thresh --pf-hdg-noise --cpf-sigma
  nmea_parser.hpp/.cpp        — разбор GPGGA (поле 9=AGL, поле 11=MSL, поле 6=quality)
  dem_loader.hpp/.cpp         — GDAL → cv::Mat float32 + GeoTransform
  terrain_correlator.hpp/.cpp — двухканальный NCC по 360 азимутам, sliding window
  kalman_filter.hpp/.cpp      — cv::KalmanFilter [x, y, vx, vy]
  contour_pf.hpp/.cpp         — Particle Filter по изолинии высоты + ridge fixes
  ridge_map.hpp/.cpp          — Гессиан ЦМР → хребты/впадины → NMS → пространственный индекс
  keypoint_db.hpp/.cpp        — БД пиков рельефа (fallback для CPF без RidgeMap)
  line_navigator.hpp/.cpp     — dead reckoning + дискретные фиксы от хребтов (LNAV-режим)
  particle_filter.hpp/.cpp    — базовый PF (предшественник ContourPF)
  visualizer.hpp/.cpp         — тепловая карта NCC (X=азимут, Y=расстояние) +
                                  цветная карта высот hillshade + стрелка TERCOM

scripts/
  generate_nmea.py            — синтетический NMEA из реального ЦМР
                                  (Орнштейн-Уленбек для скорости, heading-sigma для курса,
                                   --jammer-zone → quality=0 внутри зоны)
  fetch_satellite.py          — скачивает ESRI World Imagery тайлы, накладывает GT+TRN треки

data/dem/
  elbrus.tif                  — Copernicus GLO-30, Эльбрус N43E042 (основной)
  orenburg.tif                — вырезка Оренбург, степь (для демо на плоском рельефе)
  steppe.tif                  — синтетический ЦМР степи

data/nmea/
  elbrus_realistic.nmea       — реалистичный маршрут, dt=0.2 с, 1095 измерений
  elbrus_realistic_gt.csv     — ground truth [t, lat, lon, heading, speed, terrain_h, gps_ok]
  elbrus_sine.nmea            — синусоидальный маршрут, dt=0.2 с, 1095 измерений
  elbrus_sine_gt.csv          — ground truth
  tercom_straight.nmea        — прямолинейный, dt=1 с, 180 измерений (для TERCOM demo)
  tercom_straight_gt.csv      — ground truth
  orenburg_flight.nmea        — полёт над степью (для демо на Оренбурге)

results/
  tercom_demo/                — TERCOM: correlation_heatmap.png, trajectory.png, trn_estimate.csv
  combined_realistic/         — Combined: trajectory.png, trajectory_satellite.png,
                                           correlation_heatmap.png, trn_estimate.csv
  combined_sine/              — Combined: то же для синусоидального маршрута
```

## Скачать ЦМР (Copernicus GLO-30, бесплатно)
```bash
# Эльбрус (N43 E042)
wget "https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N43_00_E042_00_DEM/Copernicus_DSM_COG_10_N43_00_E042_00_DEM.tif" -O data/dem/elbrus.tif
```
