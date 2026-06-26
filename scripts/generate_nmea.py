#!/usr/bin/env python3
"""
Генератор синтетических NMEA GPGGA данных для тестирования TRN-алгоритма.
Симулирует полёт над ЦМР и генерирует показания радиовысотомера.

Пример:
    python3 generate_nmea.py --dem ../data/dem/region.tif \
        --lat 55.0 --lon 60.0 --azimuth 45 --speed 80 \
        --duration 120 --baro 1500 --hz 1 --noise 2.0 \
        --out ../data/nmea/flight.nmea
"""

import argparse
import math
import random
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
    """Вычислить XOR-контрольную сумму NMEA (без $ и *)."""
    cs = 0
    for ch in sentence:
        cs ^= ord(ch)
    return f"{cs:02X}"


def format_gpgga(dt: datetime, radio_alt_m: float, baro_alt_m: float) -> str:
    """
    Формирует GPGGA строку. Координаты не заполняются (алгоритм их не знает).
    radio_alt_m — высота над рельефом (поле 9, antenna altitude above MSL).
    baro_alt_m  — абсолютная высота полёта (поле геоид сепарации не используется напрямую).

    Используем поле antenna altitude = radio_alt_m (высота над рельефом).
    Барометрическую высоту выводим как undulation (поле 11).
    """
    time_str = dt.strftime("%H%M%S.") + f"{dt.microsecond // 1000:03d}"
    # lat/lon пустые — устройство не знает своих координат
    body = (
        f"GPGGA,{time_str},,,,,,,,{radio_alt_m:.1f},M,{baro_alt_m:.1f},M,,"
    )
    cs = nmea_checksum(body)
    return f"${body}*{cs}"


def geo_to_pixel(ds, lat: float, lon: float):
    """Перевести (lat, lon) в пиксельные координаты растра."""
    gt = ds.GetGeoTransform()
    # gt: [x_origin, pixel_w, 0, y_origin, 0, pixel_h]
    px = (lon - gt[0]) / gt[1]
    py = (lat - gt[3]) / gt[5]
    return px, py


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
    import struct
    vals = struct.unpack('f' * 4, row)  # float32
    v00, v10, v01, v11 = vals
    return (v00 * (1 - dx) * (1 - dy) +
            v10 * dx * (1 - dy) +
            v01 * (1 - dx) * dy +
            v11 * dx * dy)


def azimuth_to_dlat_dlon(azimuth_deg: float, distance_m: float, current_lat: float):
    """Вычислить (delta_lat, delta_lon) для движения на расстояние distance_m в направлении azimuth_deg."""
    az_rad = math.radians(azimuth_deg)
    dlat = (distance_m * math.cos(az_rad)) / 111320.0
    dlon = (distance_m * math.sin(az_rad)) / (111320.0 * math.cos(math.radians(current_lat)))
    return dlat, dlon


def main():
    parser = argparse.ArgumentParser(description="Генератор NMEA для TRN")
    parser.add_argument("--dem", required=True, help="Путь к GeoTIFF ЦМР")
    parser.add_argument("--lat", type=float, required=True, help="Стартовая широта")
    parser.add_argument("--lon", type=float, required=True, help="Стартовая долгота")
    parser.add_argument("--azimuth", type=float, default=90.0, help="Азимут полёта (градусы, 0=север)")
    parser.add_argument("--speed", type=float, default=80.0, help="Путевая скорость (м/с)")
    parser.add_argument("--duration", type=float, default=120.0, help="Длительность полёта (с)")
    parser.add_argument("--baro", type=float, default=1500.0, help="Барометрическая высота MSL (м)")
    parser.add_argument("--hz", type=float, default=1.0, help="Частота обновлений (Гц, 1-10)")
    parser.add_argument("--noise", type=float, default=1.5, help="СКО шума высотомера (м)")
    parser.add_argument("--out", default="flight.nmea", help="Выходной файл NMEA")
    parser.add_argument("--seed", type=int, default=42, help="Зерно генератора случайных чисел")
    args = parser.parse_args()

    random.seed(args.seed)

    ds = gdal.Open(args.dem)
    if ds is None:
        print(f"ERROR: не удалось открыть ЦМР: {args.dem}", file=sys.stderr)
        sys.exit(1)

    band = ds.GetRasterBand(1)
    band.SetNoDataValue(band.GetNoDataValue() or -9999)
    gt = ds.GetGeoTransform()

    dt_start = datetime(2024, 6, 15, 10, 0, 0)
    step_s = 1.0 / args.hz
    n_steps = int(args.duration * args.hz)
    dist_per_step = args.speed * step_s

    lat, lon = args.lat, args.lon

    print(f"Генерация {n_steps} точек: азимут={args.azimuth}°, скорость={args.speed} м/с, "
          f"высота={args.baro} м MSL, шум={args.noise} м")
    print(f"Старт: {lat:.6f}°N, {lon:.6f}°E → {args.out}")

    lines = []
    for i in range(n_steps):
        terrain_h = sample_dem(band, gt, lat, lon, default=0.0)
        radio_alt = args.baro - terrain_h + random.gauss(0, args.noise)
        radio_alt = max(0.0, radio_alt)

        dt_cur = dt_start + timedelta(seconds=i * step_s)
        line = format_gpgga(dt_cur, radio_alt, args.baro)
        lines.append(line)

        dlat, dlon = azimuth_to_dlat_dlon(args.azimuth, dist_per_step, lat)
        lat += dlat
        lon += dlon

    with open(args.out, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Записано {len(lines)} строк в {args.out}")
    print(f"Финиш: {lat:.6f}°N, {lon:.6f}°E")
    total_dist = args.speed * args.duration
    print(f"Пройденное расстояние: {total_dist:.0f} м")


if __name__ == "__main__":
    main()
