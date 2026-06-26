#!/usr/bin/env python3
"""
Генератор синтетических NMEA GPGGA данных для тестирования TRN-алгоритма.
Симулирует полёт над ЦМР с переменной скоростью и непрямой траекторией.

Пример:
    python3 generate_nmea.py --dem ../data/dem/elbrus.tif \
        --lat 43.50014 --lon 42.49986 --azimuth 120 \
        --speed 40 --speed-min 25 --speed-max 60 \
        --heading-sigma 3 --duration 150 --baro 6500 --noise 2.0 \
        --out ../data/nmea/elbrus_realistic.nmea
"""

import argparse
import csv
import math
import random
import struct
import sys
from datetime import datetime, timedelta

try:
    from osgeo import gdal
    import numpy as np
    gdal.UseExceptions()
except ImportError:
    print("ERROR: требуется osgeo (gdal) и numpy. Установите: sudo apt install python3-gdal python3-numpy", file=sys.stderr)
    sys.exit(1)


def nmea_checksum(sentence: str) -> str:
    cs = 0
    for ch in sentence:
        cs ^= ord(ch)
    return f"{cs:02X}"


def format_gpgga(dt: datetime, radio_alt_m: float, baro_alt_m: float,
                 quality: int = 1,
                 lat: float = 0.0, lon: float = 0.0) -> str:
    """
    Формирует GPGGA строку.
    quality: 1 = GPS доступен, 0 = GPS отсутствует (зона глушения).
    При quality>0 вписывает реальные координаты lat/lon.
    """
    time_str = dt.strftime("%H%M%S.") + f"{dt.microsecond // 1000:03d}"
    if quality > 0:
        # Перевести десятичные градусы → DDMM.MMMMM
        def deg2nmea(d, is_lon=False):
            d = abs(d)
            deg = int(d)
            minutes = (d - deg) * 60.0
            width = 3 if is_lon else 2
            return f"{deg:0{width}d}{minutes:08.5f}"
        lat_s = deg2nmea(lat)
        lat_h = "N" if lat >= 0 else "S"
        lon_s = deg2nmea(lon, is_lon=True)
        lon_h = "E" if lon >= 0 else "W"
        body = (f"GPGGA,{time_str},{lat_s},{lat_h},{lon_s},{lon_h},"
                f"{quality},08,1.0,{radio_alt_m:.1f},M,{baro_alt_m:.1f},M,,")
    else:
        body = (f"GPGGA,{time_str},,,,,{quality},,,{radio_alt_m:.1f},M,{baro_alt_m:.1f},M,,")
    cs = nmea_checksum(body)
    return f"${body}*{cs}"


def sample_dem(band, gt, lat: float, lon: float, default: float = 0.0) -> float:
    """Билинейно сэмплировать ЦМР в точке (lat, lon)."""
    px = (lon - gt[0]) / gt[1]
    py = (lat - gt[3]) / gt[5]
    x0, y0 = int(px), int(py)
    nx, ny = band.XSize, band.YSize
    if x0 < 0 or y0 < 0 or x0 >= nx - 1 or y0 >= ny - 1:
        return default
    dx, dy = px - x0, py - y0
    row = band.ReadRaster(x0, y0, 2, 2, 2, 2)
    vals = struct.unpack('f' * 4, row)
    v00, v10, v01, v11 = vals
    return (v00 * (1-dx)*(1-dy) + v10*dx*(1-dy) +
            v01*(1-dx)*dy     + v11*dx*dy)


def move(lat: float, lon: float, azimuth_deg: float, distance_m: float):
    """Вычислить новые (lat, lon) при движении на distance_m в направлении azimuth_deg."""
    az = math.radians(azimuth_deg)
    dlat = (distance_m * math.cos(az)) / 111320.0
    dlon = (distance_m * math.sin(az)) / (111320.0 * math.cos(math.radians(lat)))
    return lat + dlat, lon + dlon


def inside_jammer(lat: float, lon: float, jammer_zone) -> bool:
    """Проверить, находится ли точка внутри зоны глушения."""
    if jammer_zone is None:
        return False
    lat1, lon1, lat2, lon2 = jammer_zone
    return (min(lat1, lat2) <= lat <= max(lat1, lat2) and
            min(lon1, lon2) <= lon <= max(lon1, lon2))


def main():
    parser = argparse.ArgumentParser(description="Генератор NMEA для TRN")
    parser.add_argument("--dem",          required=True,        help="Путь к GeoTIFF ЦМР")
    parser.add_argument("--lat",          type=float, required=True, help="Стартовая широта")
    parser.add_argument("--lon",          type=float, required=True, help="Стартовая долгота")
    parser.add_argument("--azimuth",      type=float, default=90.0,  help="Начальный азимут полёта (°, 0=север)")
    parser.add_argument("--speed",        type=float, default=40.0,  help="Номинальная путевая скорость (м/с)")
    parser.add_argument("--speed-min",    type=float, default=25.0,  help="Минимальная скорость (м/с)")
    parser.add_argument("--speed-max",    type=float, default=60.0,  help="Максимальная скорость (м/с)")
    parser.add_argument("--heading-sigma",type=float, default=0.0,   help="СКО поворота курса (°/шаг), 0 = прямо")
    parser.add_argument("--duration",     type=float, default=120.0, help="Длительность полёта (с)")
    parser.add_argument("--baro",         type=float, default=1500.0,help="Барометрическая высота MSL (м)")
    parser.add_argument("--hz",           type=float, default=1.0,   help="Частота обновлений (Гц)")
    parser.add_argument("--noise",        type=float, default=2.0,   help="СКО шума высотомера (м)")
    parser.add_argument("--jammer-zone",  type=str,   default=None,
                        help="Зона глушения GPS: lat1,lon1,lat2,lon2")
    parser.add_argument("--out",          default="flight.nmea",     help="Выходной файл NMEA")
    parser.add_argument("--seed",         type=int,   default=42,    help="Зерно генератора")
    parser.add_argument("--sin-exp-amplitude", type=float, default=0.0,
                        help="Амплитуда формы sin(x)+0.1e^x (м боковое отклонение). 0=выкл")
    args = parser.parse_args()

    random.seed(args.seed)

    jammer_zone = None
    if args.jammer_zone:
        parts = [float(x) for x in args.jammer_zone.split(',')]
        if len(parts) != 4:
            print("ERROR: --jammer-zone должен быть lat1,lon1,lat2,lon2", file=sys.stderr)
            sys.exit(1)
        jammer_zone = tuple(parts)

    ds = gdal.Open(args.dem)
    if ds is None:
        print(f"ERROR: не удалось открыть ЦМР: {args.dem}", file=sys.stderr)
        sys.exit(1)

    band = ds.GetRasterBand(1)
    gt = ds.GetGeoTransform()

    dt_start = datetime(2024, 6, 15, 10, 0, 0)
    step_s   = 1.0 / args.hz
    n_steps  = int(args.duration * args.hz)

    v_nom = args.speed
    v_min = args.speed_min
    v_max = args.speed_max
    speed = v_nom

    lat, lon     = args.lat, args.lon
    heading      = args.azimuth
    ou_theta     = 0.1   # скорость возврата к номиналу (Орнштейн–Уленбек)
    ou_sigma     = 2.0   # диффузия скорости (м/с / √с)

    nmea_lines = []
    gt_rows    = []   # ground truth

    print(f"Генерация {n_steps} точек: v_nom={v_nom} м/с, "
          f"v=[{v_min},{v_max}] м/с, heading_σ={args.heading_sigma}°, "
          f"baro={args.baro} м, шум={args.noise} м", file=sys.stderr)
    print(f"Старт: {lat:.6f}°N, {lon:.6f}°E → {args.out}", file=sys.stderr)

    for i in range(n_steps):
        # Рельеф и показание высотомера
        terrain_h = sample_dem(band, gt, lat, lon, default=0.0)
        radio_alt = args.baro - terrain_h + random.gauss(0, args.noise)
        radio_alt = max(0.0, radio_alt)

        # Флаг GPS (0 в зоне глушения)
        quality = 0 if inside_jammer(lat, lon, jammer_zone) else 1

        dt_cur = dt_start + timedelta(seconds=i * step_s)
        nmea_lines.append(format_gpgga(dt_cur, radio_alt, args.baro, quality, lat, lon))

        # Ground truth
        gt_rows.append({
            "t_s":     round(i * step_s, 3),
            "lat":     round(lat, 8),
            "lon":     round(lon, 8),
            "heading": round(heading, 2),
            "speed_mps": round(speed, 3),
            "terrain_h": round(terrain_h, 2),
            "gps_ok":  quality,
        })

        # Обновить скорость (Орнштейн–Уленбек → тянется к v_nom)
        speed += -ou_theta * (speed - v_nom) * step_s + \
                  ou_sigma * math.sqrt(step_s) * random.gauss(0, 1)
        speed = max(v_min, min(v_max, speed))

        # Расстояние шага нужно до обновления курса (для sin-exp)
        dist_step = speed * step_s

        # Обновить курс
        if args.sin_exp_amplitude > 0:
            # Предписанный курс из производной sin(x)+0.1*e^x
            # x пробегает [-1.5, 4.0] по мере продвижения маршрута
            x_i = -1.5 + 5.5 * (i / max(n_steps - 1, 1))
            dydx = math.cos(x_i) + 0.1 * math.exp(x_i)
            dx_per_step = 5.5 / max(n_steps - 1, 1)
            lateral_m = args.sin_exp_amplitude * dydx * dx_per_step  # м/шаг вбок
            heading = args.azimuth + math.degrees(math.atan2(lateral_m, dist_step))
            if args.heading_sigma > 0:
                heading += random.gauss(0, args.heading_sigma)
        else:
            if args.heading_sigma > 0:
                heading += random.gauss(0, args.heading_sigma)
        heading %= 360.0

        lat, lon = move(lat, lon, heading, dist_step)

    # Записать NMEA
    with open(args.out, "w") as f:
        f.write("\n".join(nmea_lines) + "\n")

    # Записать ground truth CSV
    gt_path = args.out.replace(".nmea", "_gt.csv")
    with open(gt_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=gt_rows[0].keys())
        writer.writeheader()
        writer.writerows(gt_rows)

    total_dist = sum(r["speed_mps"] * step_s for r in gt_rows)
    print(f"Записано {len(nmea_lines)} строк NMEA → {args.out}", file=sys.stderr)
    print(f"Ground truth      → {gt_path}", file=sys.stderr)
    print(f"Финиш: {lat:.6f}°N, {lon:.6f}°E", file=sys.stderr)
    print(f"Пройденное расстояние: {total_dist:.0f} м", file=sys.stderr)
    if jammer_zone:
        jammed = sum(1 for r in gt_rows if r["gps_ok"] == 0)
        print(f"Точек в зоне глушения: {jammed}/{n_steps}", file=sys.stderr)


if __name__ == "__main__":
    main()
