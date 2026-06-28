#!/usr/bin/env bash
# Десктопный UI — только apt, без Streamlit
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

missing=()
python3 -c "import numpy" 2>/dev/null || missing+=(python3-numpy)
python3 -c "import PIL" 2>/dev/null || missing+=(python3-pil python3-pil.imagetk)
python3 -c "import tkinter" 2>/dev/null || missing+=(python3-tk)

if ((${#missing[@]})); then
  echo "Установите:"
  echo "  sudo apt install ${missing[*]}"
  exit 1
fi

if [[ ! -f test_source/manifest.ini ]]; then
  python3 scripts/bootstrap_test_source.py
fi

if [[ ! -f test_source/map.tif ]]; then
  python3 scripts/make_test_map.py || {
    echo "Нет test_source/map.tif."
    echo "  python3 scripts/make_test_map.py   # из map_full.tif или GT"
    echo "  или скачайте Copernicus GLO-30 → test_source/map_full.tif"
    exit 1
  }
fi

exec python3 app/trn_viewer_tk.py
