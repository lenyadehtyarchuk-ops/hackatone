#!/usr/bin/env python3
"""Интерактивный просмотр результатов TRN для формата 3-го чекпоинта."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import streamlit as st

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from checkpoint_io import (  # noqa: E402
    CheckpointConfig,
    compute_errors,
    estimate_rows_with_local,
    load_estimate_csv,
    load_manifest,
    load_reference_csv,
    render_trajectory_map,
    run_trn,
)

st.set_page_config(page_title="TRN Checkpoint Viewer", layout="wide")
st.title("TRN — проверка 3-го чекпоинта")
st.caption(
    "Вход: heights.txt (AGL, м) · map.tif (GeoTIFF) · manifest.ini · reference.csv (эталон)"
)

DEFAULT_SOURCE = ROOT / "test_source"
TRN_BIN = ROOT / "build" / "trn"


def save_uploaded(uploaded, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(uploaded.getbuffer())


def run_pipeline(source_dir: Path) -> tuple[CheckpointConfig, list, list, dict, object]:
    cfg = load_manifest(source_dir)
    if not cfg.map_path.exists():
        raise FileNotFoundError(f"Нет карты: {cfg.map_path}")

    if TRN_BIN.exists():
        est_path = run_trn(cfg, TRN_BIN)
    elif (cfg.output_dir / "trn_estimate.csv").exists():
        st.warning("Бинарник build/trn не найден — используем сохранённый trn_estimate.csv")
        est_path = cfg.output_dir / "trn_estimate.csv"
    else:
        raise FileNotFoundError(
            "Соберите ./build/trn или положите готовый output/trn_estimate.csv"
        )

    estimate = load_estimate_csv(est_path)
    estimate = estimate_rows_with_local(estimate, cfg)

    reference = []
    if cfg.reference_path and cfg.reference_path.exists():
        reference = load_reference_csv(cfg.reference_path)

    errors = compute_errors(reference, estimate) if reference else {}

    img = render_trajectory_map(
        cfg.map_path, reference, estimate, cfg,
        show_reference=True, show_estimate=True, show_endpoints=True,
    )
    return cfg, reference, estimate, errors, img


# ── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Источник данных")
    mode = st.radio("Режим", ["test_source (по умолчанию)", "Загрузить свои файлы"])

    work_dir = DEFAULT_SOURCE

    if mode == "Загрузить свои файлы":
        st.subheader("Файлы")
        up_manifest = st.file_uploader("manifest.ini", type=["ini", "txt"])
        up_heights = st.file_uploader("heights.txt", type=["txt"])
        up_map = st.file_uploader("map.tif (GeoTIFF)", type=["tif", "tiff"])
        up_ref = st.file_uploader("reference.csv (эталон, опционально)", type=["csv"])

        if st.button("Применить загрузку", type="primary"):
            tmp = Path(tempfile.mkdtemp(prefix="trn_upload_"))
            if up_manifest:
                save_uploaded(up_manifest, tmp / "manifest.ini")
            if up_heights:
                save_uploaded(up_heights, tmp / "heights.txt")
            if up_map:
                save_uploaded(up_map, tmp / "map.tif")
            if up_ref:
                save_uploaded(up_ref, tmp / "reference.csv")
            st.session_state["upload_dir"] = str(tmp)
            st.success(f"Данные сохранены во временный каталог")

        if "upload_dir" in st.session_state:
            work_dir = Path(st.session_state["upload_dir"])
    else:
        if not DEFAULT_SOURCE.exists():
            st.error("Нет test_source/. Запустите: python3 scripts/bootstrap_test_source.py")
        else:
            st.info(str(DEFAULT_SOURCE))

    st.divider()
    st.header("Слои на карте")
    show_ref = st.checkbox("Эталонная траектория", value=True)
    show_est = st.checkbox("TRN (предсказанная)", value=True)
    show_ep = st.checkbox("Начало / конец", value=True)

    st.divider()
    run_clicked = st.button("Запустить TRN", type="primary", use_container_width=True)

# ── Main ─────────────────────────────────────────────────────────────────────
if run_clicked:
    try:
        with st.spinner("Обработка…"):
            cfg, reference, estimate, errors, _ = run_pipeline(work_dir)
            st.session_state["last_result"] = {
                "cfg": cfg, "reference": reference, "estimate": estimate,
                "errors": errors, "work_dir": str(work_dir),
            }
    except Exception as e:
        st.error(str(e))

elif "last_result" not in st.session_state and DEFAULT_SOURCE.exists() and mode.startswith("test_source"):
    # Автозагрузка метрик из готового output без перезапуска TRN
    try:
        cfg = load_manifest(DEFAULT_SOURCE)
        est_path = cfg.output_dir / "trn_estimate.csv"
        if est_path.exists() and cfg.map_path.exists():
            estimate = estimate_rows_with_local(load_estimate_csv(est_path), cfg)
            reference = (load_reference_csv(cfg.reference_path)
                         if cfg.reference_path and cfg.reference_path.exists() else [])
            errors = compute_errors(reference, estimate) if reference else {}
            st.session_state["last_result"] = {
                "cfg": cfg, "reference": reference, "estimate": estimate,
                "errors": errors, "work_dir": str(DEFAULT_SOURCE),
            }
    except Exception:
        pass

if "last_result" in st.session_state:
    data = st.session_state["last_result"]
    cfg: CheckpointConfig = data["cfg"]
    reference = data["reference"]
    estimate = data["estimate"]
    errors = data["errors"]

    # Перерисовать с текущими слоями без повторного TRN
    img = render_trajectory_map(
        cfg.map_path, reference, estimate, cfg,
        show_reference=show_ref, show_estimate=show_est, show_endpoints=show_ep,
    )

    col_map, col_metrics = st.columns([3, 1])

    with col_map:
        st.subheader("Карта траектории")
        st.image(img, channels="BGR", use_container_width=True)
        st.caption(
            "⭐ белая звезда — старт (manifest) · "
            "жёлтая линия — эталон · оранжевая — TRN · "
            "S/E — начало/конец"
        )
        if (cfg.output_dir / "correlation_heatmap.png").exists():
            with st.expander("Тепловая карта корреляции"):
                st.image(str(cfg.output_dir / "correlation_heatmap.png"), use_container_width=True)

    with col_metrics:
        st.subheader("Параметры полёта")
        st.write(f"**Старт (лок.):** ({cfg.start_x_m:.0f}, {cfg.start_y_m:.0f}) м")
        st.write(f"**Старт (глоб.):** {cfg.start_lat:.6f}°N, {cfg.start_lon:.6f}°E")
        st.write(f"**Курс:** {cfg.heading_deg:.1f}°")
        st.write(f"**Baro MSL:** {cfg.baro_alt_m:.0f} м")
        st.write(f"**Скорость:** {cfg.speed_mps:.1f} м/с")
        st.write(f"**Измерений:** {len(estimate)}")

        st.divider()
        st.subheader("Отклонение от эталона")
        if errors:
            st.metric("RMSE", f"{errors['rmse_m']:.1f} м")
            st.metric("Финальное смещение", f"{errors['final_error_m']:.1f} м")
            st.metric("Среднее", f"{errors['mean_error_m']:.1f} м")
            st.metric("Максимум", f"{errors['max_error_m']:.1f} м")
        else:
            st.info("reference.csv не задан — метрики недоступны")

        if estimate:
            last = estimate[-1]
            st.divider()
            st.subheader("Финальная точка TRN")
            st.write(f"**Lat/Lon:** {float(last['lat']):.6f}, {float(last['lon']):.6f}")
            st.write(f"**Лок.:** ({float(last['x_m']):.1f}, {float(last['y_m']):.1f}) м")
            st.write(f"**Курс:** {float(last.get('heading_deg', 0)):.1f}°")
            st.write(f"**NCC:** {float(last.get('ncc', 0)):.3f}")

    with st.expander("Как мы понимаем входные данные судей"):
        st.markdown("""
**heights.txt** — последовательность показаний *радиовысотомера* (AGL, метры над рельефом),
по одному числу на строку. Частота задаётся `sample_hz` в manifest.ini.

**map.tif** — фрагмент цифровой модели рельефа (GeoTIFF). Алгоритм строит эталонные
профили высот рельефа: `terrain = baro_alt_m − AGL`.

**manifest.ini** — метаданные:
- `start_x_m`, `start_y_m` — локальные координаты старта (метры, X=восток, Y=север)
  относительно `origin_lat` / `origin_lon`;
- `heading_deg` — направление движения от севера;
- `baro_alt_m` — абсолютная высота полёта MSL;
- `speed_mps` — предполагаемая скорость для TERCOM.

**reference.csv** — эталонная траектория (для оценки точности на демо; судьи используют свой эталон).

**Выход:** локальные (`x_m`, `y_m`) и глобальные (`lat`, `lon`) координаты + PNG траектории.
        """)

else:
    st.info("Нажмите «Запустить TRN» в боковой панели или положите данные в test_source/")
