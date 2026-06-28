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

### Адаптивный режим: FLAT vs MOUNTAIN
При входе в зону глушения алгоритм сэмплирует ЦМР в радиусе 3 км и вычисляет `local_sigma` (СКО высот):
- `local_sigma < 50 м` → **FLAT mode** (степь): DR-якорь с мягкой PF-коррекцией
- `local_sigma ≥ 50 м` → **MOUNTAIN mode**: стандартный ContourPF

**Почему PF плохо работает на степи:** изолиния `baro−AGL = const` на плоском рельефе покрывает полосу шириной несколько км. Частицы имеют одинаковый вес везде на этой полосе → среднее по облаку = случайная точка. Neff остаётся высоким, но это ложная уверенность.

**FLAT mode — DR/PF blend со step-filter:**
1. Dead Reckoning строится от последней GPS-точки по курсу и `yaw_rate` из GPS-истории
2. DR-курс обновляется: `hdg = hdg_gps_base + yaw_rate·dt + 0.3·neff_ratio·(pf_hdg − hdg_base)` — PF влияет слабо (30%)
3. Позиция: `pos = neff_ratio·est_pf + (1−neff_ratio)·dr` — Neff-взвешенный blend
4. **Step-filter**: если предложенный шаг > 2.5× ожидаемого ИЛИ отклонение курса > 25° → отбросить PF, взять DR
5. Скорость обновляется из PF когда `neff_ratio > 0.4` и шаг прошёл фильтр

**Итог:** на степи траектория следует GPS-курсу; PF вносит мягкую коррекцию не захватывая управление.

---

## Формат входных данных (чекпоинт жюри)

```
./build/trn --source <dir> --combined ...
```
Директория `<dir>` должна содержать:
- `manifest.ini` — параметры полёта
- `heights.txt` — AGL-высоты в метрах, по одной на строку
- `map.tif` — ЦМР в формате GeoTIFF

**manifest.ini:**
```ini
[flight]
origin_lat = 43.503370
origin_lon = 42.492140
heading_deg = 120
speed_mps = 40
; baro_alt_m = 6500   ; опционально — если не задан, авто: DEM(start) + AGL[0]
```

**Автооценка baro:** если `baro_alt_m` не задан в manifest, вычисляется как `DEM(start_lat, start_lon) + AGL[0]`. Работает при условии что дрон в момент старта находится в начальной точке.

---

## Формат вывода CSV

Оба режима (TERCOM и Combined) пишут `trn_estimate.csv`:

```
timestamp_s,found_lat,found_lon,x_m,y_m,heading_deg,speed_mps,neff
```
- `x_m, y_m` — локальные координаты в метрах от стартовой точки (восток/север)
- `x_m = (lon − start_lon) · 111320 · cos(start_lat)`
- `y_m = (lat − start_lat) · 111320`

**Итоговый вывод `=== ИТОГ ===`:**
```
Старт:    lat=... lon=...  (x=0 м, y=0 м)
Финиш:    lat=... lon=...  (x=XXXX м, y=XXXX м)
Пройдено:  XXXX м
Скорость:  XX.X м/с (оценка PF)
Азимут:    XXX.X°
```

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

## Демонстрационные запуски

### Combined — реалистичный маршрут (горы, переменная скорость, изгибы курса)
```bash
./build/trn \
  --dem data/dem/elbrus.tif \
  --nmea data/nmea/elbrus_realistic.nmea \
  --baro 6500 --lat 43.503370 --lon 42.492140 --azimuth 120 \
  --speed 40 --out results/combined_realistic \
  --jammer-zone 43.4639,42.5139,43.4925,42.5690 \
  --combined --lnav-thresh 0.08 --pf-hdg-noise 1.5 --cpf-sigma 20
```

### Combined — синусоидальный маршрут (горы)
```bash
./build/trn \
  --dem data/dem/elbrus.tif \
  --nmea data/nmea/elbrus_sine.nmea \
  --baro 6500 --lat 43.503370 --lon 42.492140 --azimuth 120 \
  --speed 40 --out results/combined_sine \
  --jammer-zone 43.4639,42.5139,43.4925,42.5690 \
  --combined --lnav-thresh 0.08 --pf-hdg-noise 1.5 --cpf-sigma 20
```

### Combined — mostly jammed (94.7% зона глушения, горы)
```bash
./build/trn \
  --dem data/dem/elbrus.tif \
  --nmea data/nmea/elbrus_mostly_jammed.nmea \
  --baro 6500 --lat 43.503370 --lon 42.492140 --azimuth 120 \
  --speed 40 --out results/mostly_jammed \
  --jammer-zone 43.4639,42.5139,43.4925,42.5690 \
  --combined --lnav-thresh 0.08 --pf-hdg-noise 1.5 --cpf-sigma 20
```

### Combined — Оренбург (степь, FLAT mode)
```bash
./build/trn \
  --dem data/dem/orenburg.tif \
  --nmea data/nmea/orenburg_demo.nmea \
  --baro 300 --lat 51.537 --lon 55.033 --speed 40 \
  --out results/orenburg_demo \
  --jammer-zone 51.510,55.056,51.565,55.093 \
  --combined --lnav-thresh 0.08 --pf-hdg-noise 1.5 --cpf-sigma 20
# FLAT mode (local_sigma=22 м) → DR/PF blend, az=90°
```

### Спутниковые изображения
```bash
python3 scripts/fetch_satellite.py \
  --dem data/dem/elbrus.tif \
  --gt data/nmea/elbrus_realistic_gt.csv \
  --trn results/combined_realistic/trn_estimate.csv \
  --jammer-zone 43.4639,42.5139,43.4925,42.5690 \
  --out results/combined_realistic/trajectory_satellite.png --zoom 13
```

### Генерация NMEA
```bash
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

| Сценарий | Алгоритм | Режим | Азимут истинный | Найден | Время |
|---|---|---|---|---|---|
| Реалистичный (изгибы) | Combined | MOUNTAIN | ~120° | 127° | 1.7 с |
| Синусоидальный | Combined | MOUNTAIN | ~120° | 120° | 1.7 с |
| Mostly jammed (94.7%) | Combined | MOUNTAIN | 120° | 120° | 1.7 с |
| Оренбург (степь) | Combined | FLAT | 90° | 90°, Δy=−0.3 м | 1.7 с |

Вывод: Combined в 100× быстрее TERCOM; FLAT mode корректно обрабатывает степной рельеф через DR-якорь.

---

## Структура проекта
```
src/
  main.cpp                    — CLI: --dem --nmea --baro --lat --lon --speed
                                      --min-profile --radius --az-step --sliding
                                      --combined --lnav-thresh --pf-hdg-noise --cpf-sigma
                                      --source (формат чекпоинта жюри)
  nmea_parser.hpp/.cpp        — разбор GPGGA (поле 9=AGL, поле 11=MSL, поле 6=quality)
  dem_loader.hpp/.cpp         — GDAL → cv::Mat float32 + GeoTransform
  terrain_correlator.hpp/.cpp — двухканальный NCC по 360 азимутам, sliding window
  kalman_filter.hpp/.cpp      — cv::KalmanFilter [x, y, vx, vy]
  contour_pf.hpp/.cpp         — Particle Filter по изолинии высоты + ridge fixes
  ridge_map.hpp/.cpp          — Гессиан ЦМР → хребты/впадины → NMS → пространственный индекс
  keypoint_db.hpp/.cpp        — БД пиков рельефа (fallback для CPF без RidgeMap)
  line_navigator.hpp/.cpp     — dead reckoning + дискретные фиксы от хребтов (LNAV-режим)
  particle_filter.hpp/.cpp    — базовый PF (предшественник ContourPF)
  checkpoint_loader.hpp/.cpp  — чтение manifest.ini + heights.txt → NmeaFix[]
  visualizer.hpp/.cpp         — тепловая карта NCC + цветная карта высот hillshade +
                                  стрелка TERCOM + az/v подпись на всех изображениях

scripts/
  generate_nmea.py            — синтетический NMEA из реального ЦМР
                                  (Орнштейн-Уленбек для скорости, heading-sigma для курса,
                                   --jammer-zone → quality=0 внутри зоны)
  fetch_satellite.py          — скачивает ESRI World Imagery тайлы, накладывает GT+TRN треки,
                                  az/v подпись рядом со стартом

data/dem/
  elbrus.tif                  — Copernicus GLO-30, Эльбрус N43E042 (основной)
  orenburg.tif                — вырезка Оренбург, степь (FLAT mode demo)
  steppe.tif                  — синтетический ЦМР степи

data/nmea/
  elbrus_realistic.nmea       — реалистичный маршрут, dt=0.2 с, 1095 измерений
  elbrus_realistic_gt.csv     — ground truth [t, lat, lon, heading, speed, terrain_h, gps_ok]
  elbrus_sine.nmea            — синусоидальный маршрут, dt=0.2 с, 1095 измерений
  elbrus_sine_gt.csv          — ground truth
  elbrus_mostly_jammed.nmea   — 94.7% GPS-denied, dt=0.2 с
  elbrus_mostly_jammed_gt.csv — ground truth
  orenburg_demo.nmea          — полёт над степью, az=90°, FLAT mode demo
  orenburg_demo_gt.csv        — ground truth

results/
  combined_realistic/         — Combined MOUNTAIN: trajectory.png, trajectory_satellite.png,
                                                    trn_estimate.csv
  combined_sine/              — Combined MOUNTAIN: то же для синусоидального маршрута
  mostly_jammed/              — Combined MOUNTAIN: 94.7% GPS-denied demo
  orenburg_demo/              — Combined FLAT: степной рельеф, DR/PF blend demo
```

## Скачать ЦМР (Copernicus GLO-30, бесплатно)
```bash
# Эльбрус (N43 E042)
wget "https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N43_00_E042_00_DEM/Copernicus_DSM_COG_10_N43_00_E042_00_DEM.tif" -O data/dem/elbrus.tif
```
