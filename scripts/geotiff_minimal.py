"""Минимальное чтение/запись GeoTIFF: geotransform + растр (без GDAL)."""

from __future__ import annotations

import struct
import zlib
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

# GeoTIFF tags
GT_MODEL_PIXEL_SCALE = 33550
GT_MODEL_TIEPOINT = 33922
SAMPLE_FORMAT = 339
BITS_PER_SAMPLE = 258
IMAGE_WIDTH = 256
IMAGE_LENGTH = 257
STRIP_OFFSETS = 273
STRIP_BYTE_COUNTS = 279
TILE_WIDTH = 322
TILE_LENGTH = 323
TILE_OFFSETS = 324
TILE_BYTE_COUNTS = 325
COMPRESSION = 259
PREDICTOR = 317
PHOTOMETRIC = 262
SAMPLES_PER_PIXEL = 277
PLANAR_CONFIG = 284


def _read_ifd(data: bytes, offset: int) -> tuple[list[tuple[int, int, int, int, int]], int]:
    if offset + 2 > len(data):
        return [], 0
    n = struct.unpack_from("<H", data, offset)[0]
    entries = []
    pos = offset + 2
    for _ in range(n):
        tag, typ, cnt, val = struct.unpack_from("<HHII", data, pos)
        entries.append((tag, typ, cnt, val, pos))
        pos += 12
    next_ifd = struct.unpack_from("<I", data, pos)[0]
    return entries, next_ifd


def _tag_value(data: bytes, typ: int, cnt: int, val: int) -> bytes:
    # TIFF: 1=BYTE, 2=ASCII, 3=SHORT, 4=LONG, 5=RATIONAL, 12=DOUBLE
    sizes = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 12: 8}
    need = sizes.get(typ, 1) * cnt
    if need <= 4:
        return struct.pack("<I", val)[:need]
    return data[val: val + need]


def _tag_uint_array(data: bytes, typ: int, cnt: int, val: int) -> tuple[int, ...]:
    raw = _tag_value(data, typ, cnt, val)
    if typ == 3:
        return struct.unpack(f"<{cnt}H", raw[: 2 * cnt])
    if typ == 4:
        return struct.unpack(f"<{cnt}I", raw[: 4 * cnt])
    raise ValueError(f"Неожиданный тип массива TIFF: {typ}")


def _parse_ifd_tags(data: bytes, entries) -> dict:
    tags: dict = {}
    for tag, typ, cnt, val, _ in entries:
        raw = _tag_value(data, typ, cnt, val)
        if tag in (IMAGE_WIDTH, IMAGE_LENGTH, BITS_PER_SAMPLE, SAMPLE_FORMAT,
                   COMPRESSION, PREDICTOR, TILE_WIDTH, TILE_LENGTH):
            tags[tag] = struct.unpack("<H", raw[:2])[0] if typ == 3 else struct.unpack("<I", raw[:4])[0]
        elif tag in (STRIP_OFFSETS, STRIP_BYTE_COUNTS, TILE_OFFSETS, TILE_BYTE_COUNTS):
            tags[tag] = _tag_uint_array(data, typ, cnt, val)
        elif tag == GT_MODEL_PIXEL_SCALE and cnt >= 3:
            tags[tag] = struct.unpack("<3d", raw[:24])
        elif tag == GT_MODEL_TIEPOINT and cnt >= 6:
            tags[tag] = struct.unpack("<6d", raw[:48])
    return tags


def _geotransform_from_tags(tags: dict) -> Tuple[float, float, float, float, float, float]:
    scale = tags.get(GT_MODEL_PIXEL_SCALE)
    tie = tags.get(GT_MODEL_TIEPOINT)
    if scale and tie:
        return (tie[3], scale[0], 0.0, tie[4], 0.0, -abs(scale[1]))
    return (0.0, 1.0, 0.0, 0.0, 0.0, -1.0)


def _decode_float32(raw: bytes, w: int, h: int, predictor: int) -> np.ndarray:
    arr = np.frombuffer(raw, dtype=np.float32).reshape(h, w).copy()
    if predictor == 3:
        # TIFF floating-point predictor: горизонтальное differencing по строкам
        arr = np.cumsum(arr, axis=1)
    return arr


def _decompress_strip(raw: bytes, compression: int) -> bytes:
    if compression == 1:
        return raw
    if compression == 8:
        return zlib.decompress(raw)
    raise ValueError(f"Неподдерживаемое сжатие TIFF: {compression}")


def read_geotransform(path: Path) -> Optional[Tuple[float, float, float, float, float, float]]:
    data = path.read_bytes()
    if data[:2] not in (b"II", b"MM"):
        return None
    fmt_i = "<I" if data[:2] == b"II" else ">I"
    ifd_off = struct.unpack(fmt_i, data[4:8])[0]
    entries, _ = _read_ifd(data, ifd_off)
    tags = _parse_ifd_tags(data, entries)
    gt = _geotransform_from_tags(tags)
    if GT_MODEL_PIXEL_SCALE not in tags:
        return None
    return gt


def read_elevation(path: Path) -> Tuple[np.ndarray, Tuple[float, ...]]:
    """Вернуть (elev float32, geotransform). Strip или tiled, DEFLATE + predictor 3."""
    data = path.read_bytes()
    if data[:2] not in (b"II", b"MM"):
        raise ValueError(f"Не TIFF: {path}")
    fmt_i = "<I" if data[:2] == b"II" else ">I"
    ifd_off = struct.unpack(fmt_i, data[4:8])[0]
    entries, _ = _read_ifd(data, ifd_off)
    tags = _parse_ifd_tags(data, entries)

    w = tags.get(IMAGE_WIDTH, 0)
    h = tags.get(IMAGE_LENGTH, 0)
    bps = tags.get(BITS_PER_SAMPLE, 32)
    fmt_code = tags.get(SAMPLE_FORMAT, 3)
    compression = tags.get(COMPRESSION, 1)
    predictor = tags.get(PREDICTOR, 1)
    gt = _geotransform_from_tags(tags)

    if fmt_code != 3 or bps != 32:
        raise ValueError("Ожидается float32 GeoTIFF")

    if TILE_OFFSETS in tags:
        tile_w = tags[TILE_WIDTH]
        tile_h = tags[TILE_LENGTH]
        offsets = tags[TILE_OFFSETS]
        counts = tags[TILE_BYTE_COUNTS]
        n_tx = (w + tile_w - 1) // tile_w
        arr = np.zeros((h, w), dtype=np.float32)
        for idx, (off, cnt) in enumerate(zip(offsets, counts)):
            ty = idx // n_tx
            tx = idx % n_tx
            y0, x0 = ty * tile_h, tx * tile_w
            th = min(tile_h, h - y0)
            tw = min(tile_w, w - x0)
            decoded = _decompress_strip(data[off: off + cnt], compression)
            # Copernicus хранит тайлы фиксированного размера tile_w×tile_h (с паддингом)
            tile = _decode_float32(decoded, tile_w, tile_h, predictor)
            arr[y0: y0 + th, x0: x0 + tw] = tile[:th, :tw]
    elif STRIP_OFFSETS in tags and STRIP_BYTE_COUNTS in tags:
        strip_off = tags[STRIP_OFFSETS]
        strip_cnt = tags[STRIP_BYTE_COUNTS]
        off = strip_off[0] if isinstance(strip_off, tuple) else strip_off
        cnt = strip_cnt[0] if isinstance(strip_cnt, tuple) else strip_cnt
        decoded = _decompress_strip(data[off: off + cnt], compression)
        arr = _decode_float32(decoded, w, h, predictor)
    else:
        raise ValueError(f"Не удалось прочитать растр из {path}")

    arr = np.where(arr < -9000, np.nan, arr)
    return arr, gt


def write_geotiff(
    path: Path,
    elev: np.ndarray,
    origin_lon: float,
    origin_lat: float,
    pixel_size: float,
) -> None:
    """Записать float32 strip GeoTIFF (без сжатия) — читается UI без GDAL."""
    path = Path(path)
    h, w = elev.shape
    pixel_w = pixel_size
    pixel_h = -abs(pixel_size)
    raw = np.ascontiguousarray(elev, dtype=np.float32).tobytes()

    scale = struct.pack("<3d", pixel_w, abs(pixel_h), 0.0)
    tie = struct.pack("<6d", 0.0, 0.0, 0.0, origin_lon, origin_lat, 0.0)

    # IFD layout: header 8 + IFD + values + raster
    n_tags = 12
    ifd_size = 2 + n_tags * 12 + 4
    header_size = 8
    scale_off = header_size + ifd_size
    tie_off = scale_off + len(scale)
    strip_off = tie_off + len(tie)
    data_start = strip_off + 4  # strip offset value stored inline, raster follows

    # Recompute: strip data at data_start
    strip_data_off = header_size + ifd_size + len(scale) + len(tie)

    buf = bytearray()
    buf += b"II"
    buf += struct.pack("<H", 42)
    buf += struct.pack("<I", header_size)

    def entry(tag: int, typ: int, cnt: int, val: int) -> bytes:
        return struct.pack("<HHII", tag, typ, cnt, val)

    ifd = struct.pack("<H", n_tags)
    ifd += entry(IMAGE_WIDTH, 3, 1, w)
    ifd += entry(IMAGE_LENGTH, 3, 1, h)
    ifd += entry(BITS_PER_SAMPLE, 3, 1, 32)
    ifd += entry(COMPRESSION, 3, 1, 1)
    ifd += entry(PHOTOMETRIC, 3, 1, 1)
    ifd += entry(SAMPLES_PER_PIXEL, 3, 1, 1)
    ifd += entry(SAMPLE_FORMAT, 3, 1, 3)
    ifd += entry(PLANAR_CONFIG, 3, 1, 1)
    ifd += entry(STRIP_OFFSETS, 4, 1, strip_data_off)
    ifd += entry(STRIP_BYTE_COUNTS, 4, 1, len(raw))
    ifd += entry(GT_MODEL_PIXEL_SCALE, 12, 3, header_size + ifd_size)
    ifd += entry(GT_MODEL_TIEPOINT, 12, 6, header_size + ifd_size + len(scale))
    ifd += struct.pack("<I", 0)

    buf += ifd
    buf += scale
    buf += tie
    buf += raw

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(buf)


def crop_geotiff(
    src: Path,
    dst: Path,
    lon_min: float,
    lat_min: float,
    lon_max: float,
    lat_max: float,
) -> Tuple[int, int]:
    """Вырезать фрагмент по WGS84 и сохранить UI-совместимый GeoTIFF."""
    elev, gt = read_elevation(src)
    px0 = max(0, int((lon_min - gt[0]) / gt[1]))
    px1 = min(elev.shape[1], int(np.ceil((lon_max - gt[0]) / gt[1])))
    py0 = max(0, int((lat_max - gt[3]) / gt[5]))
    py1 = min(elev.shape[0], int(np.ceil((lat_min - gt[3]) / gt[5])))
    if px1 <= px0 or py1 <= py0:
        raise ValueError("Пустая область обрезки")

    crop = elev[py0:py1, px0:px1].copy()
    crop = np.where(np.isnan(crop), -9999.0, crop)
    origin_lon = gt[0] + px0 * gt[1]
    origin_lat = gt[3] + py0 * gt[5]
    write_geotiff(dst, crop, origin_lon, origin_lat, abs(gt[1]))
    return crop.shape[1], crop.shape[0]
