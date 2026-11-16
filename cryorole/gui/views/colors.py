"""SLD colour mapping shared by the 3D, 2D and legend widgets."""

from __future__ import annotations

from functools import lru_cache

import numpy as np


@lru_cache(maxsize=8)
def colormap_lut(name: str, size: int = 256) -> np.ndarray:
    """``(size, 4)`` uint8 RGBA lookup table of a Matplotlib colormap."""

    from matplotlib import colormaps

    return (colormaps[name](np.linspace(0.0, 1.0, size)) * 255).round().astype(np.ubyte)


def values_to_rgba(values: np.ndarray, vmin: float, vmax: float, lut: np.ndarray) -> np.ndarray:
    """Map values to RGBA through ``lut`` with clipping to ``[vmin, vmax]``; NaN becomes transparent."""

    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    scaled = np.zeros(values.shape, dtype=float)
    scaled[finite] = (values[finite] - vmin) / (vmax - vmin)
    index = np.clip(scaled * (len(lut) - 1), 0, len(lut) - 1).round().astype(np.intp)
    rgba = lut[index]
    rgba[~finite] = 0
    return rgba
