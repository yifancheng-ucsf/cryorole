"""Array-native landscape model for production NPZ workflows."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LandscapeArrays:
    """Compact landscape arrays without pandas object-coordinate columns."""

    particle_key: np.ndarray
    coordinates_analysis: np.ndarray
    coordinates_display: np.ndarray | None
    sld_unfloored: np.ndarray
    sld_raw: np.ndarray
    sld_display: np.ndarray
    sld_was_floored: np.ndarray
    sld_local_k_mean: np.ndarray
    sld_effective_local_k_mean: np.ndarray
    sld_distance_floor: np.ndarray
    sld_display_is_outlier: np.ndarray | None = None
    ref_source_row_id: np.ndarray | None = None
    mov_source_row_id: np.ndarray | None = None
    coordinates_canonical: np.ndarray | None = None
    canonical_transform: np.ndarray | None = None

    def __post_init__(self) -> None:
        particle_key = np.asarray(self.particle_key).astype(str)
        n_points = len(particle_key)
        object.__setattr__(self, "particle_key", particle_key)
        object.__setattr__(
            self,
            "coordinates_analysis",
            _coordinate_array(self.coordinates_analysis, "coordinates_analysis", n_points),
        )
        if self.coordinates_display is not None:
            object.__setattr__(
                self,
                "coordinates_display",
                _coordinate_array(
                    self.coordinates_display,
                    "coordinates_display",
                    n_points,
                ),
            )
        if self.coordinates_canonical is not None:
            object.__setattr__(
                self,
                "coordinates_canonical",
                _coordinate_array(
                    self.coordinates_canonical,
                    "coordinates_canonical",
                    n_points,
                ),
            )
        for field_name in (
            "sld_unfloored",
            "sld_raw",
            "sld_display",
            "sld_local_k_mean",
            "sld_effective_local_k_mean",
            "sld_distance_floor",
        ):
            object.__setattr__(
                self,
                field_name,
                _one_dimensional_array(getattr(self, field_name), field_name, n_points, float),
            )
        object.__setattr__(
            self,
            "sld_was_floored",
            _one_dimensional_array(self.sld_was_floored, "sld_was_floored", n_points, bool),
        )
        display_outlier = (
            np.zeros(n_points, dtype=bool)
            if self.sld_display_is_outlier is None
            else self.sld_display_is_outlier
        )
        object.__setattr__(
            self,
            "sld_display_is_outlier",
            _one_dimensional_array(
                display_outlier,
                "sld_display_is_outlier",
                n_points,
                bool,
            ),
        )
        for field_name in ("ref_source_row_id", "mov_source_row_id"):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(
                    self,
                    field_name,
                    _one_dimensional_array(value, field_name, n_points, np.int64),
                )
        if self.canonical_transform is not None:
            transform = np.asarray(self.canonical_transform, dtype=float)
            if transform.shape != (3, 3) or not np.isfinite(transform).all():
                raise ValueError("canonical_transform must be a finite shape (3, 3) array")
            object.__setattr__(self, "canonical_transform", transform)

    @property
    def n_points(self) -> int:
        """Return particle count."""

        return int(len(self.particle_key))


def _coordinate_array(value, field_name: str, n_points: int) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != (n_points, 3):
        raise ValueError(f"{field_name} must have shape ({n_points}, 3)")
    if not np.isfinite(array).all():
        raise ValueError(f"{field_name} contains non-finite values")
    return array


def _one_dimensional_array(value, field_name: str, n_points: int, dtype) -> np.ndarray:
    array = np.asarray(value, dtype=dtype)
    if array.shape != (n_points,):
        raise ValueError(f"{field_name} must have shape ({n_points},)")
    if dtype is not bool and not np.isfinite(array).all():
        raise ValueError(f"{field_name} contains non-finite values")
    return array
