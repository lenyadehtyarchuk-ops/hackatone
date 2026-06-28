#!/usr/bin/env python3
"""
Формат входных данных 3-го чекпоинта (как мы его понимаем):

  heights.txt   — показания радиовысотомера (AGL, м), по одному значению на строку.
                  Частота задаётся sample_hz в manifest.ini (1–10 Гц).
  map.tif       — фрагмент ЦМР в GeoTIFF (Copernicus / SRTM / ALOS).
  manifest.ini  — метаданные полёта:
      start_x_m, start_y_m  — локальные координаты старта (м, ось X=восток, Y=север)
                              относительно origin_lat / origin_lon
      heading_deg           — курс относительно севера (°)
      baro_alt_m            — барометрическая высота MSL (м)
      speed_mps             — предполагаемая путевая скорость (м/с)
  reference.csv — эталонная траектория (x_m, y_m, lat, lon) для оценки точности

Алгоритм TRN получает AGL + baro → профиль высот рельефа и ищет совпадение на ЦМР.
"""

from __future__ import annotations

import configparser
import csv
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

try:
    from osgeo import gdal
    gdal.UseExceptions()
    HAS_GDAL = True
except ImportError:
    gdal = None
    HAS_GDAL = False

try:
    import cv2
    HAS_CV2 = True
except ImportError:
    HAS_CV2 = False

from geotiff_minimal import read_elevation, read_geotransform


EARTH_R = 6378137.0


@dataclass
class CheckpointConfig:
    source_dir: Path
    map_path: Path
    heights_path: Path
    reference_path: Optional[Path]
    origin_lat: float
    origin_lon: float
    start_x_m: float
    start_y_m: float
    heading_deg: float
    baro_alt_m: float
    speed_mps: float
    sample_hz: float
    output_dir: Path
    min_profile: int = 35
    search_radius_m: float = 3000.0

    @property
    def start_lat(self) -> float:
        return meters_to_latlon(self.start_x_m, self.start_y_m,
                                self.origin_lat, self.origin_lon)[0]

    @property
    def start_lon(self) -> float:
        return meters_to_latlon(self.start_x_m, self.start_y_m,
                                self.origin_lat, self.origin_lon)[1]


def meters_to_latlon(x_m: float, y_m: float,
                     ref_lat: float, ref_lon: float) -> tuple[float, float]:
    lat = ref_lat + y_m / EARTH_R * (180.0 / math.pi)
    lon = ref_lon + x_m / (EARTH_R * math.cos(math.radians(ref_lat))) * (180.0 / math.pi)
    return lat, lon


def latlon_to_meters(lat: float, lon: float,
                     ref_lat: float, ref_lon: float) -> tuple[float, float]:
    x = (lon - ref_lon) * (math.pi / 180.0) * EARTH_R * math.cos(math.radians(ref_lat))
    y = (lat - ref_lat) * (math.pi / 180.0) * EARTH_R
    return x, y


def load_manifest(source_dir: Path) -> CheckpointConfig:
    source_dir = Path(source_dir).resolve()
    ini_path = source_dir / "manifest.ini"
    if not ini_path.exists():
        raise FileNotFoundError(f"Нет manifest.ini в {source_dir}")

    cp = configparser.ConfigParser()
    cp.read(ini_path, encoding="utf-8")
    if "flight" not in cp:
        raise ValueError("manifest.ini должен содержать секцию [flight]")

    f = cp["flight"]
    ref = f.get("reference", "").strip()
    out = f.get("output", "output").strip()

    return CheckpointConfig(
        source_dir=source_dir,
        map_path=source_dir / f.get("map", "map.tif"),
        heights_path=source_dir / f.get("heights", "heights.txt"),
        reference_path=(source_dir / ref) if ref else None,
        origin_lat=float(f.get("origin_lat", "0")),
        origin_lon=float(f.get("origin_lon", "0")),
        start_x_m=float(f.get("start_x_m", "0")),
        start_y_m=float(f.get("start_y_m", "0")),
        heading_deg=float(f.get("heading_deg", "0")),
        baro_alt_m=float(f.get("baro_alt_m", "1500")),
        speed_mps=float(f.get("speed_mps", "40")),
        sample_hz=float(f.get("sample_hz", "1")),
        output_dir=source_dir / out,
        min_profile=int(f.get("min_profile", "35")),
        search_radius_m=float(f.get("search_radius_m", "3000")),
    )


def read_heights(path: Path) -> list[float]:
    values = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        values.append(float(line.replace(",", ".")))
    if not values:
        raise ValueError(f"Файл высот пуст: {path}")
    return values


def heights_to_nmea(heights: list[float], baro_alt_m: float,
                    sample_hz: float, out_path: Path) -> None:
    """Конвертировать AGL-профиль в NMEA GPGGA (GPS quality=0)."""
    dt = 1.0 / sample_hz
    lines = []
    t0 = 36000.0  # 10:00:00 UTC
    for i, agl in enumerate(heights):
        ts = t0 + i * dt
        hh = int(ts // 3600) % 24
        mm = int((ts % 3600) // 60)
        ss = ts % 60
        time_str = f"{hh:02d}{mm:02d}{ss:06.3f}"
        body = f"GPGGA,{time_str},,,,,0,,,{agl:.1f},M,{baro_alt_m:.1f},M,,"
        cs = 0
        for ch in body:
            cs ^= ord(ch)
        lines.append(f"${body}*{cs:02X}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_reference_csv(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def load_estimate_csv(path: Path) -> list[dict]:
    return load_reference_csv(path)


def compute_errors(reference: list[dict], estimate: list[dict]) -> dict:
    """RMSE и финальное смещение по парам точек (по индексу)."""
    n = min(len(reference), len(estimate))
    if n == 0:
        return {"rmse_m": float("nan"), "final_error_m": float("nan"),
                "max_error_m": float("nan"), "mean_error_m": float("nan")}

    errors = []
    for i in range(n):
        ref = reference[i]
        est = estimate[i]
        if "x_m" in ref and "x_m" in est:
            dx = float(est["x_m"]) - float(ref["x_m"])
            dy = float(est["y_m"]) - float(ref["y_m"])
        else:
            rlat, rlon = float(ref["lat"]), float(ref["lon"])
            elat = float(est.get("found_lat", est.get("lat", 0)))
            elon = float(est.get("found_lon", est.get("lon", 0)))
            olat = float(ref.get("origin_lat", rlat))
            olon = float(ref.get("origin_lon", rlon))
            rx, ry = latlon_to_meters(rlat, rlon, olat, olon)
            ex, ey = latlon_to_meters(elat, elon, olat, olon)
            dx, dy = ex - rx, ey - ry
        errors.append(math.hypot(dx, dy))

    final = errors[-1]
    rmse = math.sqrt(sum(e * e for e in errors) / n)
    return {
        "rmse_m": rmse,
        "final_error_m": final,
        "max_error_m": max(errors),
        "mean_error_m": sum(errors) / n,
        "n_points": n,
    }


def run_trn(cfg: CheckpointConfig, trn_bin: Path) -> Path:
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    heights = read_heights(cfg.heights_path)
    nmea_path = cfg.output_dir / "flight.nmea"
    heights_to_nmea(heights, cfg.baro_alt_m, cfg.sample_hz, nmea_path)

    cmd = [
        str(trn_bin),
        "--dem", str(cfg.map_path),
        "--nmea", str(nmea_path),
        "--baro", str(cfg.baro_alt_m),
        "--lat", str(cfg.start_lat),
        "--lon", str(cfg.start_lon),
        "--speed", str(cfg.speed_mps),
        "--azimuth", str(cfg.heading_deg),
        "--min-profile", str(cfg.min_profile),
        "--radius", str(cfg.search_radius_m),
        "--out", str(cfg.output_dir),
    ]
    print("[checkpoint] Запуск:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)
    return cfg.output_dir / "trn_estimate.csv"


def render_dem_hillshade(map_path: Path):
    """BGR uint8 изображение ЦМР + geotransform."""
    import numpy as np

    if HAS_GDAL:
        ds = gdal.Open(str(map_path))
        if ds is None:
            raise FileNotFoundError(f"Не удалось открыть {map_path}")
        band = ds.GetRasterBand(1)
        elev = band.ReadAsArray().astype(np.float32)
        gt = ds.GetGeoTransform()
        ds = None
    else:
        elev, gt = read_elevation(map_path)

    nodata = -9999.0
    elev = np.where(elev < -9000, np.nan, elev)
    valid = np.nanmin(elev), np.nanmax(elev)
    t = (elev - valid[0]) / max(valid[1] - valid[0], 1e-6)
    t = np.nan_to_num(t, nan=0.0)

    r = np.clip(np.where(t < 0.35, 60 + t / 0.35 * 80,
                 np.where(t < 0.65, 140 + (t - 0.35) / 0.30 * 80,
                          220 + (t - 0.65) / 0.35 * 35)), 0, 255).astype(np.uint8)
    g = np.clip(np.where(t < 0.35, 130 + t / 0.35 * 50,
                 np.where(t < 0.65, 180 - (t - 0.35) / 0.30 * 30,
                          150 - (t - 0.65) / 0.35 * 30)), 0, 255).astype(np.uint8)
    b = np.clip(np.where(t < 0.35, 50 - t / 0.35 * 20,
                 np.where(t < 0.65, 30 + (t - 0.35) / 0.30 * 20,
                          50 + (t - 0.65) / 0.35 * 70)), 0, 255).astype(np.uint8)
    vis = np.stack([b, g, r], axis=-1)

    gy, gx = np.gradient(elev, abs(gt[5]) * 111320.0, abs(gt[1]) * 111320.0)
    slope = np.arctan(np.sqrt(gx * gx + gy * gy))
    aspect = np.arctan2(-gx, gy)
    sun_az, sun_alt = math.radians(315), math.radians(45)
    shade = np.sin(sun_alt) * np.cos(slope) + np.cos(sun_alt) * np.sin(slope) * np.cos(sun_az - aspect)
    shade = np.clip(0.4 + 0.6 * np.maximum(0, shade), 0, 1)
    vis = (vis.astype(np.float32) * shade[..., None]).astype(np.uint8)
    return vis, gt


def geo_to_pixel(gt, lat: float, lon: float) -> tuple[float, float]:
    px = (lon - gt[0]) / gt[1]
    py = (lat - gt[3]) / gt[5]
    return px, py


def render_trajectory_map(
    map_path: Path,
    reference: Optional[list[dict]],
    estimate: Optional[list[dict]],
    cfg: CheckpointConfig,
    *,
    show_reference: bool = True,
    show_estimate: bool = True,
    show_endpoints: bool = True,
):
    import numpy as np

    vis, gt = render_dem_hillshade(map_path)
    h, w = vis.shape[:2]

    def collect_pts(rows):
        pts = []
        for row in rows:
            lat = float(row.get("lat", row.get("found_lat", 0)))
            lon = float(row.get("lon", row.get("found_lon", 0)))
            px, py = geo_to_pixel(gt, lat, lon)
            if 0 <= px < w and 0 <= py < h:
                pts.append((int(round(px)), int(round(py))))
        return pts

    if HAS_CV2:
        import cv2

        def draw_poly(rows, color, thickness=2):
            pts = collect_pts(rows)
            if len(pts) >= 2:
                cv2.polylines(vis, [np.array(pts, dtype=np.int32)], False,
                              color, thickness, cv2.LINE_AA)

        if show_reference and reference:
            draw_poly(reference, (0, 255, 255), 2)
        if show_estimate and estimate:
            draw_poly(estimate, (0, 180, 255), 3)
        if show_endpoints:
            for rows, color in [(reference, (0, 255, 0)), (estimate, (0, 0, 255))]:
                if not rows:
                    continue
                for row, label in [(rows[0], "S"), (rows[-1], "E")]:
                    lat = float(row.get("lat", row.get("found_lat", 0)))
                    lon = float(row.get("lon", row.get("found_lon", 0)))
                    px, py = geo_to_pixel(gt, lat, lon)
                    if 0 <= px < w and 0 <= py < h:
                        cv2.circle(vis, (int(px), int(py)), 8, color, -1, cv2.LINE_AA)
                        cv2.putText(vis, label, (int(px) + 10, int(py) - 6),
                                    cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), 2, cv2.LINE_AA)
        slat, slon = cfg.start_lat, cfg.start_lon
        px, py = geo_to_pixel(gt, slat, slon)
        if 0 <= px < w and 0 <= py < h:
            cv2.drawMarker(vis, (int(px), int(py)), (255, 255, 255),
                           cv2.MARKER_STAR, 14, 2, cv2.LINE_AA)
    else:
        from PIL import Image, ImageDraw

        # BGR → RGB для PIL
        img = Image.fromarray(vis[:, :, ::-1])
        draw = ImageDraw.Draw(img)

        def draw_poly(rows, color, width=2):
            pts = collect_pts(rows)
            if len(pts) >= 2:
                draw.line(pts, fill=color, width=width)

        if show_reference and reference:
            draw_poly(reference, (255, 255, 0), 2)
        if show_estimate and estimate:
            draw_poly(estimate, (255, 140, 0), 3)
        if show_endpoints:
            for rows, color in [(reference, (0, 255, 0)), (estimate, (255, 0, 0))]:
                if not rows:
                    continue
                for row, label in [(rows[0], "S"), (rows[-1], "E")]:
                    lat = float(row.get("lat", row.get("found_lat", 0)))
                    lon = float(row.get("lon", row.get("found_lon", 0)))
                    px, py = geo_to_pixel(gt, lat, lon)
                    if 0 <= px < w and 0 <= py < h:
                        r = 8
                        draw.ellipse([px - r, py - r, px + r, py + r], fill=color)
                        draw.text((px + 10, py - 12), label, fill=(255, 255, 255))
        slat, slon = cfg.start_lat, cfg.start_lon
        px, py = geo_to_pixel(gt, slat, slon)
        if 0 <= px < w and 0 <= py < h:
            r = 10
            draw.polygon([
                (px, py - r), (px + r * 0.3, py + r * 0.3),
                (px - r * 0.3, py + r * 0.3),
            ], fill=(255, 255, 255))
        vis = np.array(img)[:, :, ::-1]

    return vis


def estimate_rows_with_local(estimate: list[dict], cfg: CheckpointConfig) -> list[dict]:
    out = []
    for row in estimate:
        lat = float(row["found_lat"])
        lon = float(row["found_lon"])
        x, y = latlon_to_meters(lat, lon, cfg.origin_lat, cfg.origin_lon)
        r = dict(row)
        r["lat"] = lat
        r["lon"] = lon
        r["x_m"] = x
        r["y_m"] = y
        out.append(r)
    return out
