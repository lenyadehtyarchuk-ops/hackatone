#!/usr/bin/env python3
"""Десктопный просмотрщик TRN (Tkinter) — без Streamlit, только apt-пакеты."""

from __future__ import annotations

import shutil
import sys
import tempfile
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

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

try:
    from PIL import Image, ImageTk
except ImportError:
    print("Установите: sudo apt install python3-pil python3-pil.imagetk python3-tk python3-numpy")
    sys.exit(1)

import numpy as np

DEFAULT_SOURCE = ROOT / "test_source"
TRN_BIN = ROOT / "build" / "trn"


class TrnViewerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("TRN Checkpoint Viewer")
        self.geometry("1200x800")
        self.source_dir = DEFAULT_SOURCE
        self.cfg: CheckpointConfig | None = None
        self.reference: list = []
        self.estimate: list = []
        self.errors: dict = {}

        self.var_ref = tk.BooleanVar(value=True)
        self.var_est = tk.BooleanVar(value=True)
        self.var_ep = tk.BooleanVar(value=True)

        self._build_ui()
        self._load_default()

    def _build_ui(self):
        top = ttk.Frame(self, padding=8)
        top.pack(fill=tk.X)
        ttk.Button(top, text="Открыть каталог…", command=self._pick_dir).pack(side=tk.LEFT, padx=4)
        ttk.Button(top, text="Загрузить файлы…", command=self._pick_files).pack(side=tk.LEFT, padx=4)
        ttk.Button(top, text="Запустить TRN", command=self._run_trn).pack(side=tk.LEFT, padx=4)
        self.lbl_dir = ttk.Label(top, text=str(self.source_dir))
        self.lbl_dir.pack(side=tk.LEFT, padx=12)

        opts = ttk.Frame(self, padding=(8, 0))
        opts.pack(fill=tk.X)
        ttk.Checkbutton(opts, text="Эталонная траектория", variable=self.var_ref,
                        command=self._redraw).pack(side=tk.LEFT, padx=8)
        ttk.Checkbutton(opts, text="TRN (предсказанная)", variable=self.var_est,
                        command=self._redraw).pack(side=tk.LEFT, padx=8)
        ttk.Checkbutton(opts, text="Начало / конец", variable=self.var_ep,
                        command=self._redraw).pack(side=tk.LEFT, padx=8)

        body = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        body.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        self.canvas_frame = ttk.Frame(body)
        self.metrics_frame = ttk.Frame(body, width=280)
        body.add(self.canvas_frame, weight=3)
        body.add(self.metrics_frame, weight=1)

        self.canvas = tk.Canvas(self.canvas_frame, bg="#1e1e1e", highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        self._photo = None

        self.metrics = tk.Text(self.metrics_frame, wrap=tk.WORD, width=36, font=("Sans", 10))
        self.metrics.pack(fill=tk.BOTH, expand=True)

    def _pick_dir(self):
        d = filedialog.askdirectory(initialdir=str(self.source_dir))
        if d:
            self.source_dir = Path(d)
            self.lbl_dir.config(text=str(self.source_dir))
            self._load_data()

    def _pick_files(self):
        tmp = Path(tempfile.mkdtemp(prefix="trn_upload_"))
        for name, patterns in [
            ("manifest.ini", [("INI", "*.ini"), ("All", "*.*")]),
            ("heights.txt", [("Text", "*.txt"), ("All", "*.*")]),
            ("map.tif", [("GeoTIFF", "*.tif"), ("All", "*.*")]),
        ]:
            p = filedialog.askopenfilename(title=f"Выберите {name}", filetypes=patterns)
            if p:
                shutil.copy(p, tmp / name)
        ref = filedialog.askopenfilename(title="reference.csv (опционально)",
                                         filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if ref:
            shutil.copy(ref, tmp / "reference.csv")
        if (tmp / "manifest.ini").exists() and (tmp / "heights.txt").exists():
            self.source_dir = tmp
            self.lbl_dir.config(text=str(tmp))
            self._load_data()
        else:
            messagebox.showwarning("Недостаточно файлов",
                                   "Нужны как минимум manifest.ini и heights.txt")

    def _load_default(self):
        if DEFAULT_SOURCE.exists():
            self._load_data()

    def _load_data(self):
        try:
            self.cfg = load_manifest(self.source_dir)
            est_path = self.cfg.output_dir / "trn_estimate.csv"
            if not est_path.exists():
                messagebox.showinfo("TRN", "Нет output/trn_estimate.csv — нажмите «Запустить TRN»")
                self.estimate = []
            else:
                self.estimate = estimate_rows_with_local(load_estimate_csv(est_path), self.cfg)
            self.reference = []
            if self.cfg.reference_path and self.cfg.reference_path.exists():
                self.reference = load_reference_csv(self.cfg.reference_path)
            self.errors = compute_errors(self.reference, self.estimate) if self.reference else {}
            self._redraw()
            self._update_metrics()
        except Exception as e:
            messagebox.showerror("Ошибка", str(e))

    def _run_trn(self):
        if not TRN_BIN.exists():
            messagebox.showerror("TRN", f"Соберите бинарник: {TRN_BIN}")
            return
        try:
            run_trn(self.cfg or load_manifest(self.source_dir), TRN_BIN)
            self._load_data()
        except Exception as e:
            messagebox.showerror("TRN", str(e))

    def _redraw(self):
        if not self.cfg or not self.cfg.map_path.exists():
            return
        bgr = render_trajectory_map(
            self.cfg.map_path, self.reference, self.estimate, self.cfg,
            show_reference=self.var_ref.get(),
            show_estimate=self.var_est.get(),
            show_endpoints=self.var_ep.get(),
        )
        rgb = bgr[:, :, ::-1]
        img = Image.fromarray(rgb)
        cw = max(self.canvas.winfo_width(), 400)
        ch = max(self.canvas.winfo_height(), 300)
        img.thumbnail((cw, ch), Image.Resampling.LANCZOS)
        self._photo = ImageTk.PhotoImage(img)
        self.canvas.delete("all")
        self.canvas.create_image(cw // 2, ch // 2, image=self._photo, anchor=tk.CENTER)

    def _update_metrics(self):
        if not self.cfg:
            return
        lines = [
            "=== Параметры ===",
            f"Старт (лок.): ({self.cfg.start_x_m:.0f}, {self.cfg.start_y_m:.0f}) м",
            f"Старт: {self.cfg.start_lat:.6f}°N, {self.cfg.start_lon:.6f}°E",
            f"Курс: {self.cfg.heading_deg:.1f}°",
            f"Baro: {self.cfg.baro_alt_m:.0f} м",
            f"Скорость: {self.cfg.speed_mps:.1f} м/с",
            f"Точек TRN: {len(self.estimate)}",
            "",
            "=== Отклонение от эталона ===",
        ]
        if self.errors:
            lines += [
                f"RMSE: {self.errors['rmse_m']:.1f} м",
                f"Финальное: {self.errors['final_error_m']:.1f} м",
                f"Среднее: {self.errors['mean_error_m']:.1f} м",
                f"Максимум: {self.errors['max_error_m']:.1f} м",
            ]
        else:
            lines.append("(reference.csv не задан)")
        if self.estimate:
            last = self.estimate[-1]
            lines += [
                "",
                "=== Финал TRN ===",
                f"Lat/Lon: {float(last['lat']):.6f}, {float(last['lon']):.6f}",
                f"Лок.: ({float(last['x_m']):.1f}, {float(last['y_m']):.1f}) м",
                f"NCC: {float(last.get('ncc', 0)):.3f}",
            ]
        self.metrics.delete("1.0", tk.END)
        self.metrics.insert(tk.END, "\n".join(lines))


def main():
    app = TrnViewerApp()
    app.bind("<Configure>", lambda _: app._redraw())
    app.mainloop()


if __name__ == "__main__":
    main()
