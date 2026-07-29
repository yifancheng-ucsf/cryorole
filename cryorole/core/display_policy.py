"""Shared display-only policy for static visualization and animation."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

import numpy as np


EULER_PROJECTIONS = (
    ("alpha_beta", 0, 1),
    ("beta_gamma", 1, 2),
    ("alpha_gamma", 0, 2),
)
EULER_AXIS_NAMES = ("alpha", "beta", "gamma")
LEGACY_DISPLAY_STYLE = {
    "color_map": "rainbow_r",
    "point_size": 1.0,
    "point_alpha": 1.0,
    "figure_size_2d": (20.0, 9.0),
    "figure_size_3d": (8.0, 7.0),
    "colorbar_position": "bottom",
    "aspect": "equal",
    "sort_points_by_color": "ascending",
    "axis_limits": {
        "alpha": (-180.0, 180.0),
        "beta": (-180.0, 180.0),
        "gamma": (-180.0, 180.0),
        "x": (-3.14, 3.14),
        "y": (-3.14, 3.14),
        "z": (-3.14, 3.14),
    },
}


@dataclass(frozen=True)
class ColorScale:
    """Resolved display color bounds and their provenance."""

    vmin: float | None
    vmax: float | None
    vmin_source: str
    vmax_source: str


def resolve_display_indices(
    values,
    *,
    threshold: float | None = None,
    top_fraction: float | None = None,
    coordinates: np.ndarray | None = None,
    axis_names: tuple[str, str, str] | None = None,
    range_bounds: Mapping[str, tuple[float | None, float | None] | None] | None = None,
    wraparound: bool = False,
) -> np.ndarray:
    """Resolve display rows with inclusive threshold and stable top-fraction ties."""

    density = np.asarray(values, dtype=float)
    if density.ndim != 1 or not np.isfinite(density).all():
        raise ValueError("display density values must be a finite one-dimensional array")
    if threshold is not None and top_fraction is not None:
        raise ValueError("threshold and top_fraction are mutually exclusive")
    if threshold is not None and not np.isfinite(threshold):
        raise ValueError("threshold must be finite")
    if top_fraction is not None and (
        not np.isfinite(top_fraction) or not 0 < top_fraction <= 1
    ):
        raise ValueError("top_fraction must be in (0, 1]")

    selected = np.arange(len(density), dtype=int)
    if threshold is not None:
        selected = selected[density >= float(threshold)]
    if top_fraction is not None and selected.size:
        count = max(1, int(math.ceil(selected.size * float(top_fraction))))
        order = np.argsort(density[selected], kind="mergesort")
        selected = np.sort(selected[order[-count:]])

    bounds = dict(range_bounds or {})
    for axis, axis_bounds in bounds.items():
        if axis_bounds is None:
            continue
        for bound in axis_bounds:
            if bound is not None and not np.isfinite(bound):
                raise ValueError(f"display range for {axis} must use finite bounds")
    if bounds and selected.size:
        if coordinates is None or axis_names is None:
            raise ValueError("range_bounds require coordinates and axis_names")
        selected = selected[
            coordinate_range_mask(
                np.asarray(coordinates, dtype=float)[selected],
                axis_names=axis_names,
                range_bounds=bounds,
                wraparound=wraparound,
            )
        ]
    return np.asarray(selected, dtype=int)


def coordinate_range_mask(
    coordinates: np.ndarray,
    *,
    axis_names: tuple[str, str, str],
    range_bounds: Mapping[str, tuple[float | None, float | None] | None],
    wraparound: bool,
) -> np.ndarray:
    """Return the inclusive display-only coordinate range mask."""

    values = np.asarray(coordinates, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(axis_names):
        raise ValueError("display coordinates must have one column per axis name")
    mask = np.ones(values.shape[0], dtype=bool)
    for axis_index, axis_name in enumerate(axis_names):
        bounds = range_bounds.get(axis_name)
        if bounds is None:
            continue
        lower, upper = bounds
        axis_values = values[:, axis_index]
        if lower is not None and upper is not None and wraparound and lower > upper:
            mask &= (axis_values >= lower) | (axis_values <= upper)
            continue
        if lower is not None:
            mask &= axis_values >= lower
        if upper is not None:
            mask &= axis_values <= upper
    return mask


def resolve_color_scale(
    values,
    *,
    display_outlier_mask: np.ndarray | None,
    visual_style: str,
    color_vmin: float | None,
    color_vmax: float | None,
    display_threshold: float | None,
) -> ColorScale:
    """Resolve display color bounds, excluding flagged tail-jump outliers."""

    if color_vmin is not None and not np.isfinite(color_vmin):
        raise ValueError("color_vmin must be finite when provided")
    if color_vmax is not None and not np.isfinite(color_vmax):
        raise ValueError("color_vmax must be finite when provided")
    if color_vmin is not None and color_vmax is not None and color_vmin >= color_vmax:
        raise ValueError("color_vmin must be less than color_vmax")
    resolved_vmin = color_vmin
    resolved_vmax = color_vmax
    vmin_source = "explicit" if color_vmin is not None else "auto"
    vmax_source = "explicit" if color_vmax is not None else "auto"
    raw_values = np.asarray(values, dtype=float)
    finite_values = raw_values[np.isfinite(raw_values)]
    scale_values = finite_values
    if display_outlier_mask is not None and np.asarray(display_outlier_mask).shape == raw_values.shape:
        outliers = np.asarray(display_outlier_mask, dtype=bool)
        scale_values = raw_values[np.isfinite(raw_values) & ~outliers]
        if scale_values.size and outliers.any() and resolved_vmax is None:
            resolved_vmax = float(np.max(scale_values))
            vmax_source = "tail_jump_display_outliers"
    if visual_style == "legacy" and resolved_vmin is None and display_threshold is not None:
        resolved_vmin = float(display_threshold)
        vmin_source = "legacy_display_threshold"
    if visual_style == "legacy" and resolved_vmax is None and scale_values.size:
        resolved_vmax = math.ceil(float(scale_values.max()) * 10.0) / 10.0
        vmax_source = "legacy_color_max_ceiling"
    return ColorScale(resolved_vmin, resolved_vmax, vmin_source, vmax_source)


def sort_display_indices(indices, values, *, order: str) -> np.ndarray:
    """Sort display indices stably by color value, or preserve their order."""

    selected = np.asarray(indices, dtype=int)
    if order == "none" or selected.size == 0:
        return selected
    density = np.asarray(values, dtype=float)
    ranked = np.argsort(density[selected], kind="mergesort")
    if order == "descending":
        ranked = ranked[::-1]
    elif order != "ascending":
        raise ValueError("display sort order must be none, ascending, or descending")
    return selected[ranked]
