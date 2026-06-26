#!/usr/bin/env python3
"""
Скачивает ESRI World Imagery тайлы для области ЦМР,
сшивает их и рисует на них траекторию из ground-truth CSV.

Использование:
    python3 fetch_satellite.py \
        --dem ../data/dem/elbrus.tif \
        --gt  ../data/nmea/elbrus_realistic_gt.csv \
        --out ../results/trajectory_satellite.png \
        [--zoom 12]
"""

import argparse
import csv
import math
import os
import struct
import subprocess
import sys
import tempfile

try:
    from osgeo import gdal
    gdal.UseExceptions()
except ImportError:
    print("ERROR: требуется osgeo (gdal). sudo apt install python3-gdal", file=sys.stderr)
    sys.exit(1)

try:
    from PIL import Image, ImageDraw
except ImportError:
    print("ERROR: требуется Pillow. sudo apt install python3-pil", file=sys.stderr)
    sys.exit(1)


# --- Тайловая математика ---

def deg2tile(lat_deg: float, lon_deg: float, zoom: int):
    lat_r = math.radians(lat_deg)
    n = 2 ** zoom
    xtile = int((lon_deg + 180.0) / 360.0 * n)
    ytile = int((1.0 - math.log(math.tan(lat_r) + 1.0 / math.cos(lat_r)) / math.pi) / 2.0 * n)
    return xtile, ytile


def tile2deg(xtile: int, ytile: int, zoom: int):
    """Северо-западный угол тайла → (lat, lon)."""
    n = 2 ** zoom
    lon = xtile / n * 360.0 - 180.0
    lat_r = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
    lat = math.degrees(lat_r)
    return lat, lon


def download_tile(x: int, y: int, z: int, cache_dir: str) -> str | None:
    """Скачать тайл ESRI World Imagery через wget, вернуть путь к файлу."""
    url = (f"https://server.arcgisonline.com/ArcGIS/rest/services/"
           f"World_Imagery/MapServer/tile/{z}/{y}/{x}")
    path = os.path.join(cache_dir, f"{z}_{x}_{y}.jpg")
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        return path
    try:
        result = subprocess.run(
            ["wget", "-q", "-O", path, url],
            timeout=15, capture_output=True
        )
        if result.returncode == 0 and os.path.getsize(path) > 1000:
            return path
        else:
            os.remove(path) if os.path.exists(path) else None
            return None
    except Exception:
        return None


def latlon_to_pixel(lat: float, lon: float,
                    zoom: int, x0_tile: int, y0_tile: int,
                    tile_size: int = 256):
    """(lat, lon) → пиксель в сшитом изображении (x0_tile, y0_tile) верхний-левый тайл."""
    n = 2 ** zoom
    # Дробный тайловый индекс
    xtf = (lon + 180.0) / 360.0 * n
    lat_r = math.radians(lat)
    ytf = (1.0 - math.log(math.tan(lat_r) + 1.0 / math.cos(lat_r)) / math.pi) / 2.0 * n
    px = (xtf - x0_tile) * tile_size
    py = (ytf - y0_tile) * tile_size
    return int(px), int(py)


def get_dem_bbox(dem_path: str):
    """Вернуть (lat_min, lat_max, lon_min, lon_max) из GeoTIFF."""
    ds = gdal.Open(dem_path)
    if ds is None:
        raise RuntimeError(f"Cannot open DEM: {dem_path}")
    gt = ds.GetGeoTransform()
    w, h = ds.RasterXSize, ds.RasterYSize
    lon_min = gt[0]
    lon_max = gt[0] + w * gt[1]
    lat_max = gt[3]
    lat_min = gt[3] + h * gt[5]
    return lat_min, lat_max, lon_min, lon_max


def main():
    parser = argparse.ArgumentParser(description="Спутниковая подложка для TRN")
    parser.add_argument("--dem",  required=True, help="Путь к GeoTIFF ЦМР (bbox из него)")
    parser.add_argument("--gt",   required=True, help="Ground-truth CSV (t_s,lat,lon,...)")
    parser.add_argument("--out",  default="results/trajectory_satellite.png")
    parser.add_argument("--zoom", type=int, default=12, help="Zoom тайлов (10-14)")
    parser.add_argument("--cache", default="/tmp/sat_tiles", help="Папка кэша тайлов")
    parser.add_argument("--jammer-zone", default="", help="Зона глушения: lat1,lon1,lat2,lon2")
    parser.add_argument("--trn", default="", help="CSV с TRN-оценками (found_lat,found_lon из stdout trn)")
    args = parser.parse_args()

    os.makedirs(args.cache, exist_ok=True)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    # Читаем GT один раз: для bbox и для рисования
    gt_rows = []  # список (lat, lon, gps_ok)
    with open(args.gt, newline="") as f:
        for row in csv.DictReader(f):
            gt_rows.append((float(row["lat"]), float(row["lon"]),
                            int(row.get("gps_ok", 1))))

    if not gt_rows:
        raise RuntimeError("GT CSV пустой")

    gt_lats = [r[0] for r in gt_rows]
    gt_lons = [r[1] for r in gt_rows]

    if not gt_lats:
        raise RuntimeError("GT CSV пустой")

    MARGIN_KM = 4.0  # +30% к прежнему 3 км
    lat_margin = MARGIN_KM / 111.32
    cos_lat = math.cos(math.radians(sum(gt_lats) / len(gt_lats)))
    lon_margin = MARGIN_KM / (111.32 * cos_lat)

    lat_min = min(gt_lats) - lat_margin
    lat_max = max(gt_lats) + lat_margin
    lon_min = min(gt_lons) - lon_margin
    lon_max = max(gt_lons) + lon_margin
    print(f"Bbox маршрута+margin: lat=[{lat_min:.4f},{lat_max:.4f}] lon=[{lon_min:.4f},{lon_max:.4f}]",
          file=sys.stderr)

    # Тайловые индексы для bbox
    x_min_t, y_min_t = deg2tile(lat_max, lon_min, args.zoom)
    x_max_t, y_max_t = deg2tile(lat_min, lon_max, args.zoom)

    nx = x_max_t - x_min_t + 1
    ny = y_max_t - y_min_t + 1
    TILE = 256
    total = nx * ny
    print(f"Тайлы zoom={args.zoom}: {nx}×{ny} = {total} шт.", file=sys.stderr)

    # Скачать и сшить тайлы
    mosaic = Image.new("RGB", (nx * TILE, ny * TILE), (50, 50, 50))
    downloaded = 0
    for ty in range(y_min_t, y_max_t + 1):
        for tx in range(x_min_t, x_max_t + 1):
            path = download_tile(tx, ty, args.zoom, args.cache)
            if path:
                try:
                    tile_img = Image.open(path).convert("RGB")
                    px_off = (tx - x_min_t) * TILE
                    py_off = (ty - y_min_t) * TILE
                    mosaic.paste(tile_img, (px_off, py_off))
                    downloaded += 1
                except Exception:
                    pass

    print(f"Скачано {downloaded}/{total} тайлов", file=sys.stderr)

    # --- Зона глушения GPS (полупрозрачный красный прямоугольник) ---
    if args.jammer_zone:
        j = [float(x) for x in args.jammer_zone.split(",")]
        jx1, jy1 = latlon_to_pixel(j[0], j[1], args.zoom, x_min_t, y_min_t, TILE)
        jx2, jy2 = latlon_to_pixel(j[2], j[3], args.zoom, x_min_t, y_min_t, TILE)
        rx1, ry1 = min(jx1, jx2), min(jy1, jy2)
        rx2, ry2 = max(jx1, jx2), max(jy1, jy2)
        overlay = Image.new("RGBA", mosaic.size, (0, 0, 0, 0))
        draw_ov = ImageDraw.Draw(overlay)
        draw_ov.rectangle([rx1, ry1, rx2, ry2], fill=(255, 30, 30, 55))     # ~22% красный
        draw_ov.rectangle([rx1, ry1, rx2, ry2], outline=(255, 60, 60, 220), width=2)
        mosaic = Image.alpha_composite(mosaic.convert("RGBA"), overlay).convert("RGB")
        draw_txt = ImageDraw.Draw(mosaic)
        draw_txt.text((rx1 + 5, ry1 + 5), "GPS DENIED", fill=(255, 200, 200))

    # Строим пиксели из уже прочитанных GT-точек
    draw = ImageDraw.Draw(mosaic)
    pts = []
    gps_flags = []
    for lat, lon, gps_ok in gt_rows:
        px, py = latlon_to_pixel(lat, lon, args.zoom, x_min_t, y_min_t, TILE)
        pts.append((px, py))
        gps_flags.append(gps_ok)

    if len(pts) >= 2:
        for i in range(1, len(pts)):
            color = (255, 220, 0) if gps_flags[i] else (255, 50, 50)  # жёлтый вне зоны, красный внутри
            # Белая обводка
            draw.line([pts[i-1], pts[i]], fill=(255,255,255), width=5)
            draw.line([pts[i-1], pts[i]], fill=color, width=3)

    # Попытка загрузить шрифт побольше
    try:
        from PIL import ImageFont
        try:
            font_label = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 18)
            font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 14)
        except Exception:
            font_label = ImageFont.load_default()
            font_small = font_label
    except Exception:
        font_label = None
        font_small = None

    def draw_label(d, x, y, text, color, shadow=(0,0,0)):
        kw = {"font": font_label} if font_label else {}
        d.text((x+1, y+1), text, fill=shadow, **kw)
        d.text((x, y), text, fill=color, **kw)

    if pts:
        # Старт — зелёный
        sx, sy = pts[0]
        draw.ellipse([sx-10, sy-10, sx+10, sy+10], fill=(255,255,255))
        draw.ellipse([sx-8,  sy-8,  sx+8,  sy+8 ], fill=(30, 210, 30))
        draw_label(draw, sx + 14, sy - 10, "Start", (120, 255, 120))
        # Финиш — оранжевый квадрат
        ex, ey = pts[-1]
        draw.rectangle([ex-10, ey-10, ex+10, ey+10], fill=(255,255,255))
        draw.rectangle([ex-8,  ey-8,  ex+8,  ey+8 ], fill=(220, 80, 30))
        draw_label(draw, ex + 14, ey - 10, "End", (255, 180, 80))

    # --- TRN estimated track ---
    if args.trn and os.path.exists(args.trn):
        trn_rows = []
        with open(args.trn, newline="") as f:
            for row in csv.DictReader(f):
                trn_rows.append((float(row["found_lat"]), float(row["found_lon"])))
        if trn_rows:
            # Линия от старта GT до первого TRN фикса (период накопления)
            start_px = latlon_to_pixel(gt_rows[0][0], gt_rows[0][1], args.zoom, x_min_t, y_min_t, TILE)
            first_trn_px = latlon_to_pixel(trn_rows[0][0], trn_rows[0][1], args.zoom, x_min_t, y_min_t, TILE)
            draw.line([start_px, first_trn_px], fill=(80, 80, 255), width=5)
            draw.line([start_px, first_trn_px], fill=(120, 180, 255), width=3)
            # Соединяем TRN фиксы
            trn_pts = [first_trn_px] + [
                latlon_to_pixel(lat, lon, args.zoom, x_min_t, y_min_t, TILE)
                for lat, lon in trn_rows[1:]
            ]
            for i in range(1, len(trn_pts)):
                draw.line([trn_pts[i-1], trn_pts[i]], fill=(80, 80, 255), width=5)
                draw.line([trn_pts[i-1], trn_pts[i]], fill=(120, 180, 255), width=3)
            # Финиш TRN
            ex, ey = trn_pts[-1]
            draw.rectangle([ex-8, ey-8, ex+8, ey+8], fill=(255,255,255))
            draw.rectangle([ex-6, ey-6, ex+6, ey+6], fill=(60, 120, 255))
            draw_label(draw, ex + 12, ey - 10, "TRN end", (120, 180, 255))

    # Легенда
    try:
        kw_s = {"font": font_small} if font_small else {}
        draw.text((10, 10), "GPS TRACK (ground truth)", fill=(255,255,255), **kw_s)
        draw.text((10, 30), f"Yellow = GPS ok   Red = GPS denied   zoom={args.zoom}", fill=(200,200,200), **kw_s)
        draw.text((10, 50), "Blue = TRN estimated track", fill=(120, 180, 255), **kw_s)
    except Exception:
        pass

    # Ограничить размер
    MAX_DIM = 2048
    if mosaic.width > MAX_DIM or mosaic.height > MAX_DIM:
        scale = MAX_DIM / max(mosaic.width, mosaic.height)
        new_w = int(mosaic.width * scale)
        new_h = int(mosaic.height * scale)
        mosaic = mosaic.resize((new_w, new_h), Image.LANCZOS)

    mosaic.save(args.out)
    print(f"Сохранено: {args.out} ({mosaic.width}×{mosaic.height})", file=sys.stderr)


if __name__ == "__main__":
    main()
