#!/usr/bin/env bash
# UI через venv + Streamlit (python3-streamlit в apt нет на Python 3.14)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

VENV="$ROOT/.venv"
PY="$VENV/bin/python"
PIP="$VENV/bin/pip"
STREAMLIT="$VENV/bin/streamlit"

if [[ ! -d "$VENV" ]]; then
  echo "Создаю виртуальное окружение .venv …"
  if ! python3 -m venv "$VENV"; then
    echo ""
    echo "Нужен пакет venv. Установите:"
    echo "  sudo apt install python3.14-venv python3-pip"
    echo ""
    echo "Или запустите десктопный UI без Streamlit:"
    echo "  sudo apt install python3-tk python3-pil python3-pil.imagetk python3-numpy"
    echo "  ./run_viewer_tk.sh"
    exit 1
  fi
fi

if [[ ! -f test_source/manifest.ini ]]; then
  python3 scripts/bootstrap_test_source.py
fi

if [[ ! -f test_source/map.tif ]]; then
  python3 scripts/make_test_map.py || {
    echo "Нет test_source/map.tif — python3 scripts/make_test_map.py"
    exit 1
  }
fi

echo "Устанавливаю зависимости UI в .venv …"
"$PIP" install -q --upgrade pip
"$PIP" install -q -r requirements-ui.txt

exec "$STREAMLIT" run app/trn_viewer.py --server.headless true
