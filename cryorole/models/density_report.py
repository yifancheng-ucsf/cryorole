"""Density report model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class DensityReport:
    """Diagnostic report for SLD density computation."""

    n_points: int
    requested_k_neighbors: int
    effective_k_neighbors: int
    global_local_k_mean: float
    distance_floor: float
    distance_floor_fraction: float
    n_floored_points: int
    fraction_floored_points: float
    n_inf_sld_unfloored: int
    n_inf_sld_raw: int
    max_sld_unfloored: float
    max_sld_raw: float
    p99_sld_unfloored: float
    p99_sld_raw: float
    max_over_p99_unfloored: float
    max_over_p99_raw: float
    floored_particle_keys: tuple[Any, ...]
    floored_rows: tuple[Mapping[str, Any], ...]
    p99_5_sld_raw: float = 0.0
    sld_display_mode: str = "identity"
    sld_display_outlier_mode: str = "tail_jump"
    sld_tail_search_fraction: float = 0.01
    sld_tail_jump_factor: float = 5.0
    sld_max_display_outlier_fraction: float = 0.002
    sld_display_outlier_threshold: float | None = None
    sld_display_color_vmax: float | None = None
    largest_tail_jump_ratio: float = 1.0
    n_sld_display_outliers: int = 0
    fraction_sld_display_outliers: float = 0.0
    max_over_display_vmax: float = 1.0
    near_identity_ro_tolerance_rad: float = 1e-8
    n_near_identity_ro: int = 0
    fraction_near_identity_ro: float = 0.0
    near_duplicate_coordinate_tolerance_rad: float = 1e-8
    n_near_duplicate_coordinate_points: int = 0
    largest_near_duplicate_coordinate_cluster: int = 0
    warnings: tuple[str, ...] = ()
