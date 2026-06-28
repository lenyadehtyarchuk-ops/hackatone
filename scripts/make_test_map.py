#!/usr/bin/env python3
"""Создать test_source/map.tif, читаемый UI без python-gdal."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from checkpoint_io import load_manifest  # noqa: E402
from geotiff_minimal import crop_geotiff, read_elevation, write_geotiff  # noqa: E402

EARTH_R = 6378137.0
PIXEL_M = 30.0
DEFAULT_GT = ROOT / "data/nmea/tercom_straight_gt.csv"


def _flight_bounds(source_dir: Path, margin_m: float) -> tuple[float, float, float, float]:
    cfg = load_manifest(source_dir)
    lats, lons = [cfg.start_lat], [cfg.start_lon]
    ref_path = cfg.reference_path
    if ref_path and ref_path.exists():
        for row in csv.DictReader(ref_path.open(encoding="utf-8")):
            lats.append(float(row["lat"]))
            lons.append(float(row["lon"]))

    lat_min, lat_max = min(lats), max(lats)
    lon_min, lon_max = min(lons), max(lons)
    m = margin_m
    dlat = m / EARTH_R * (180.0 / math.pi)
    cos_lat = math.cos(math.radians((lat_min + lat_max) / 2.0))
    dlon = m / (EARTH_R * cos_lat) * (180.0 / math.pi)
    return lon_min - dlon, lat_min - dlat, lon_max + dlon, lat_max + dlat


def _synthetic_from_gt(gt_csv: Path, source_dir: Path) -> tuple[np.ndarray, float, float, float]:
    """ЦМР по terrain_h из ground truth (IDW на сетке 30 м)."""
    rows = list(csv.DictReader(gt_csv.open(encoding="utf-8")))
    if not rows or "terrain_h" not in rows[0]:
        raise ValueError("GT CSV без terrain_h")

    cfg = load_manifest(source_dir)
    lon_min, lat_min, lon_max, lat_max = _flight_bounds(source_dir, margin_m=3500.0)

    m_per_deg_lat = EARTH_R * math.pi / 180.0
    m_per_deg_lon = m_per_deg_lat * math.cos(math.radians(cfg.origin_lat))
    w = max(32, int(math.ceil((lon_max - lon_min) * m_per_deg_lon / PIXEL_M)))
    h = max(32, int(math.ceil((lat_max - lat_min) * m_per_deg_lat / PIXEL_M)))

    pts = np.array([[float(r["lon"]), float(r["lat"]), float(r["terrain_h"])] for r in rows])
    grid = np.zeros((h, w), dtype=np.float32)
    lon_step = (lon_max - lon_min) / max(w - 1, 1)
    lat_step = (lat_max - lat_min) / max(h - 1, 1)

    for j in range(h):
        lat = lat_max - j * lat_step
        for i in range(w):
            lon = lon_min + i * lon_step
            d = (pts[:, 0] - lon) ** 2 + (pts[:, 1] - lat) ** 2
            k = min(6, len(pts))
            idx = np.argpartition(d, k - 1)[:k]
            wts = 1.0 / (d[idx] + 1e-10)
            grid[j, i] = float(np.sum(pts[idx, 2] * wts) / np.sum(wts))

    return grid, lon_min, lat_max, PIXEL_M


def _crop_ok(path: Path) -> bool:
    try:
        elev, _ = read_elevation(path)
        valid = elev[(elev > 500) & (elev < 9000)]
        return valid.size > elev.size * 0.5
    except Exception:
        return False


def main() -> None:
    ap = argparse.ArgumentParser(description="Создать UI-совместимый test_source/map.tif")
    ap.add_argument("--source", default=str(ROOT / "test_source"))
    ap.add_argument("--from-dem", default="", help="Обрезать из Copernicus GeoTIFF")
    ap.add_argument("--gt", default=str(DEFAULT_GT))
    ap.add_argument("--margin-m", type=float, default=3500.0)
    ap.add_argument("--out", default="", help="По умолчанию <source>/map.tif")
    args = ap.parse_args()

    source_dir = Path(args.source)
    out_path = Path(args.out) if args.out else source_dir / "map.tif"

    if args.from_dem:
        src = Path(args.from_dem)
        if not src.exists():
            sys.exit(f"Нет файла: {src}")
        print(f"[make_test_map] Обрезка из {src}")
        bounds = _flight_bounds(source_dir, args.margin_m)
        w, h = crop_geotiff(src, out_path, *bounds)
        if not _crop_ok(out_path):
            sys.exit("Обрезка дала некорректные высоты — используйте --gt без --from-dem")
        print(f"[make_test_map] → {out_path} ({w}×{h})")
        return

    gt_path = Path(args.gt)
    if not gt_path.exists():
        sys.exit(f"Нет {gt_path}. Запустите bootstrap_test_source.py")

    print(f"[make_test_map] Синтетическая ЦМР из {gt_path.name}")
    grid, lon0, lat0, ps = _synthetic_from_gt(gt_path, source_dir)
    write_geotiff(out_path, grid, lon0, lat0, ps)
    print(f"[make_test_map] → {out_path} ({grid.shape[1]}×{grid.shape[0]})")


if __name__ == "__main__":
    main()
