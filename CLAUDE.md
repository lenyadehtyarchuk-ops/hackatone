# Terrain-Referenced Navigation (TRN)

## Задача
Определить координаты и вектор скорости БПЛА без ГНСС, используя радиовысотомер и цифровую модель рельефа (ЦМР). Алгоритм: TERCOM-подобная корреляция профиля рельефа.

## Физика
- Радиовысотомер → высота над рельефом (AGL)
- Барометрический высотомер → абсолютная высота (MSL)
- Высота рельефа = MSL − AGL

## Алгоритм
1. Накапливаем профиль измерений `M[n]` из NMEA-потока
2. Для каждого азимута θ ∈ [0°, 360°) строим эталонный профиль `R[n]` из ЦМР
3. NCC(θ, d) — нормированная кросс-корреляция с поиском смещения d
4. θ*, d* = argmax NCC → направление полёта и пройденное расстояние
5. Фильтр Калмана [x,y,vx,vy] сглаживает последовательные решения
6. Физический фильтр: отклоняем решения со скачком > 250 м/с

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

## Запуск (пример — Оренбургская область)
```bash
# Генерация тестовых данных (ЦМР берётся из data/dem/)
python3 scripts/generate_nmea.py --dem data/dem/orenburg.tif \
    --lat 51.52 --lon 55.04 --azimuth 60 --speed 80 \
    --duration 45 --baro 600 --noise 1 --out data/nmea/orenburg_flight.nmea

# Запуск TRN (центр поиска = стартовая точка дрона)
./build/trn --dem data/dem/orenburg.tif \
            --nmea data/nmea/orenburg_flight.nmea \
            --baro 600 --lat 51.52 --lon 55.04 \
            --speed 80 --min-profile 35 --radius 2000 \
            --out results/
```

## Результаты верификации
| Местность | Δh | Шум | Азимут истинный | Найден | NCC |
|---|---|---|---|---|---|
| Эльбрус (горы) | ~1000 м | σ=7 м | 120° | 120° | 0.9997 |
| Оренбург (степь, реальный ЦМР) | 154 м | σ=1 м | 60° | 60° | 0.9960 |
| Синтетическая степь | 70 м | σ=4 м | 75° | **195°** | 0.94 |

Вывод: алгоритм надёжен при σ(рельефа) > 30 м. На плоской местности возникает 180°-амбигуальность.

## Структура
```
src/
  main.cpp                    — CLI (--dem, --nmea, --baro, --lat, --lon, --speed,
                                      --min-profile, --radius, --az-step, --sliding)
  nmea_parser.hpp/.cpp        — разбор GPGGA (поле 9=AGL, поле 11=MSL)
  dem_loader.hpp/.cpp         — GDAL → cv::Mat float32 + GeoTransform
  terrain_correlator.hpp/.cpp — NCC по 360 азимутам, sliding window
  kalman_filter.hpp/.cpp      — cv::KalmanFilter [x,y,vx,vy]
  visualizer.hpp/.cpp         — тепловая карта NCC + цветная карта высот с hillshading
scripts/
  generate_nmea.py            — синтетический NMEA из реального ЦМР
results/
  trajectory.png              — траектория на цветной карте высот (hillshade)
  correlation_heatmap.png     — тепловая карта NCC (360° × смещение)
  terrain_color.png           — спутниковый снимок ESRI / цветной рельеф с маршрутом
data/dem/                     — .tif ЦМР (не в git, скачивать отдельно)
data/nmea/                    — NMEA логи полётов
```

## Скачать ЦМР (Copernicus GLO-30, бесплатно)
```bash
# Эльбрус (N43 E042)
wget "https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N43_00_E042_00_DEM/Copernicus_DSM_COG_10_N43_00_E042_00_DEM.tif" -O data/dem/elbrus.tif

# Оренбург (N51 E055) — вырезать 10×10 км участок скриптом ниже
wget "https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_N51_00_E055_00_DEM/Copernicus_DSM_COG_10_N51_00_E055_00_DEM.tif" -O data/dem/orenburg_full.tif
python3 scripts/crop_dem.py  # см. скрипт обрезки
```
