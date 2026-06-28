#!/usr/bin/env python3
"""Подготовить test_source из существующих NMEA + ground truth."""

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from checkpoint_io import latlon_to_meters  # noqa: E402


def parse_nmea_agl(path: Path) -> list[float]:
    heights = []
    for line in path.read_text(encoding="utf-8").splitlines():
        m = re.search(r",(\d+\.?\d*),M,(\d+\.?\d*),M,", line)
        if m:
            heights.append(float(m.group(1)))
    return heights


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmea", default=str(ROOT / "data/nmea/tercom_straight.nmea"))
    ap.add_argument("--gt", default=str(ROOT / "data/nmea/tercom_straight_gt.csv"))
    ap.add_argument("--out", default=str(ROOT / "test_source"))
    ap.add_argument("--map", default="map.tif")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    gt_rows = list(csv.DictReader(open(args.gt, encoding="utf-8")))
    if not gt_rows:
        sys.exit("GT пуст")

    origin_lat = float(gt_rows[0]["lat"])
    origin_lon = float(gt_rows[0]["lon"])

    ref_rows = []
    for row in gt_rows:
        lat, lon = float(row["lat"]), float(row["lon"])
        x, y = latlon_to_meters(lat, lon, origin_lat, origin_lon)
        ref_rows.append({
            "t_s": row["t_s"],
            "x_m": f"{x:.3f}",
            "y_m": f"{y:.3f}",
            "lat": row["lat"],
            "lon": row["lon"],
            "heading_deg": row.get("heading", "120"),
            "speed_mps": row.get("speed_mps", "40"),
        })

    heights = parse_nmea_agl(Path(args.nmea))
    (out / "heights.txt").write_text(
        "\n".join(f"{h:.1f}" for h in heights) + "\n", encoding="utf-8"
    )

    with (out / "reference.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=ref_rows[0].keys())
        w.writeheader()
        w.writerows(ref_rows)

    heading = float(gt_rows[0].get("heading", 120))
    speed = float(gt_rows[0].get("speed_mps", 40))
    n = len(heights)
    hz = 1.0 if n <= 200 else n / 219.0

    manifest = f"""[flight]
# Формат 3-го чекпоинта: heights.txt + map.tif + manifest.ini
map = {args.map}
heights = heights.txt
reference = reference.csv
origin_lat = {origin_lat:.8f}
origin_lon = {origin_lon:.8f}
start_x_m = 0
start_y_m = 0
heading_deg = {heading:.1f}
baro_alt_m = 6500
speed_mps = {speed:.1f}
sample_hz = {hz:.2f}
min_profile = 35
search_radius_m = 3000
output = output
"""
    (out / "manifest.ini").write_text(manifest, encoding="utf-8")
    print(f"test_source готов: {out}")
    print(f"  heights: {len(heights)} значений AGL")
    print(f"  reference: {len(ref_rows)} точек")
    print(f"  origin: {origin_lat:.6f}, {origin_lon:.6f}")

    map_path = out / args.map
    if not map_path.exists():
        import subprocess
        r = subprocess.run(
            [sys.executable, str(ROOT / "scripts/make_test_map.py"), "--source", str(out)],
            cwd=str(ROOT),
        )
        if r.returncode != 0 or not map_path.exists():
            print("  Создайте карту: python3 scripts/make_test_map.py")


if __name__ == "__main__":
    main()
