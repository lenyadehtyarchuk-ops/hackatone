#!/usr/bin/env bash
# Запуск TRN из test_source (формат 3-го чекпоинта)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

if [[ ! -f test_source/manifest.ini ]]; then
  python3 scripts/bootstrap_test_source.py
fi

if [[ ! -f test_source/map.tif ]]; then
  echo "Скачайте GeoTIFF в test_source/map.tif (см. CLAUDE.md)"
  exit 1
fi

if [[ ! -x build/trn ]]; then
  echo "Соберите: mkdir -p build && cd build && cmake .. && make -j\$(nproc)"
  exit 1
fi

./build/trn --source test_source
echo ""
echo "Результаты: test_source/output/"
echo "  trn_estimate.csv  — локальные + глобальные координаты"
echo "  trajectory.png    — визуализация"
echo "  correlation_heatmap.png"
