"""SLD density computation in analysis space.

SLD is kNN-scaled Local Density. ``sld_raw`` is the floor-stabilized production
density value; ``sld_unfloored`` is retained for diagnostics; ``sld_display`` is
a display-only view and must not overwrite raw truth.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from cryorole.models.density_report import DensityReport
from cryorole.models.landscape import Landscape
from cryorole.models.policies import DensityPolicy
from cryorole.models.ro_result import ROResult


def _require_ro_result(ro_result: ROResult) -> None:
    if not isinstance(ro_result, ROResult):
        raise TypeError("density computation requires an ROResult")
    if ro_result.data.empty:
        raise ValueError("ROResult must contain at least one row")


def _analysis_coordinates_from_ro(
    ro_result: ROResult,
    *,
    coordinate_source: str,
) -> np.ndarray:
    if coordinate_source not in ro_result.data.columns:
        raise ValueError(f"ROResult does not contain coordinate source {coordinate_source!r}")
    coords = np.stack(ro_result.data[coordinate_source].map(lambda value: np.asarray(value, dtype=float)))
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"Analysis coordinates must have shape (n, 3), got {coords.shape}")
    return coords


def _validate_density_policy(policy: DensityPolicy) -> None:
    if policy.density_metric != "sld":
        raise ValueError("cryoROLE 2.0 density_metric must be 'sld'")
    if policy.k_neighbors < 1:
        raise ValueError("k_neighbors must be at least 1")
    if policy.distance_floor_mode != "relative_to_global_local_k_mean":
        raise ValueError(
            "Only distance_floor_mode='relative_to_global_local_k_mean' is supported"
        )
    if policy.distance_floor_fraction < 0:
        raise ValueError("distance_floor_fraction must be non-negative")
    if policy.duplicate_policy != "report_only":
        raise ValueError("Only duplicate_policy='report_only' is supported")
    if policy.display_normalization_mode not in {
        "identity",
        "minmax",
    }:
        raise ValueError("display_normalization_mode must be identity or minmax")
    if policy.display_outlier_mode not in {"none", "tail_jump"}:
        raise ValueError("display_outlier_mode must be none or tail_jump")
    if not 0.0 < policy.tail_search_fraction <= 1.0:
        raise ValueError("tail_search_fraction must be > 0 and <= 1")
    if policy.tail_jump_factor <= 1.0 or not np.isfinite(policy.tail_jump_factor):
        raise ValueError("tail_jump_factor must be finite and > 1")
    if not 0.0 <= policy.max_display_outlier_fraction <= 1.0:
        raise ValueError("max_display_outlier_fraction must be >= 0 and <= 1")
    if policy.near_identity_ro_tolerance_rad < 0:
        raise ValueError("near_identity_ro_tolerance_rad must be non-negative")
    if policy.near_duplicate_coordinate_tolerance_rad < 0:
        raise ValueError("near_duplicate_coordinate_tolerance_rad must be non-negative")


def _ratio(numerator: float, denominator: float) -> float:
    if np.isinf(numerator) and np.isinf(denominator):
        return 1.0
    if np.isinf(numerator):
        return float("inf")
    if np.isinf(denominator):
        return 0.0
    if np.isclose(denominator, 0.0):
        return float("inf") if numerator > 0 else 1.0
    return float(numerator / denominator)


def compute_sld_values(
    coordinates_analysis: np.ndarray,
    *,
    k_neighbors: int,
    distance_floor_fraction: float,
) -> dict[str, np.ndarray | float]:
    """Compute SLD values and floor diagnostics in analysis space."""

    coords = np.asarray(coordinates_analysis, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"coordinates_analysis must have shape (n, 3), got {coords.shape}")
    n_points = len(coords)
    if n_points == 0:
        raise ValueError("coordinates_analysis must contain at least one point")
    if n_points == 1:
        raise ValueError("SLD requires at least two analysis coordinates")

    k = min(k_neighbors, n_points - 1)

    tree = cKDTree(coords)
    distances, _ = tree.query(coords, k=k + 1)
    neighbor_distances = distances[:, 1:]
    local_k_mean = np.mean(neighbor_distances, axis=1)
    global_local_k_mean = float(np.mean(local_k_mean))
    if np.isclose(global_local_k_mean, 0.0):
        raise ValueError(
            "SLD is undefined when all analysis coordinates are duplicated or collapsed; "
            "check particle identity matching or upstream pose metadata."
        )
    distance_floor = float(distance_floor_fraction * global_local_k_mean)
    effective_local_k_mean = np.maximum(local_k_mean, distance_floor)

    with np.errstate(divide="ignore", invalid="ignore"):
        sld_unfloored = global_local_k_mean / local_k_mean
        sld_raw = global_local_k_mean / effective_local_k_mean
    sld_unfloored = np.where(local_k_mean == 0.0, np.inf, sld_unfloored)
    sld_raw = np.where(effective_local_k_mean == 0.0, np.inf, sld_raw)

    return {
        "sld_local_k_mean": local_k_mean,
        "global_local_k_mean": global_local_k_mean,
        "sld_distance_floor": distance_floor,
        "sld_effective_local_k_mean": effective_local_k_mean,
        "sld_unfloored": sld_unfloored,
        "sld_raw": sld_raw,
        "sld_was_floored": local_k_mean < distance_floor,
    }


def normalize_sld_for_display(
    sld_raw: np.ndarray,
    *,
    mode: str,
) -> np.ndarray:
    """Return a display-only normalization of raw SLD values."""

    result = compute_sld_display_values(
        sld_raw,
        mode=mode,
    )
    return np.asarray(result["sld_display"], dtype=float)


def compute_sld_display_values(
    sld_raw: np.ndarray,
    *,
    mode: str,
    outlier_mode: str = "tail_jump",
    tail_search_fraction: float = 0.01,
    tail_jump_factor: float = 5.0,
    max_display_outlier_fraction: float = 0.002,
) -> dict[str, np.ndarray | float | int | str | None]:
    """Return display-only SLD values and tail-jump outlier diagnostics."""

    values = np.asarray(sld_raw, dtype=float)
    if mode == "identity":
        display = values.copy()
    if mode == "minmax":
        minimum = float(np.min(values))
        maximum = float(np.max(values))
        if np.isclose(minimum, maximum):
            display = np.zeros_like(values)
        else:
            display = (values - minimum) / (maximum - minimum)
    elif mode != "identity":
        raise ValueError(f"Unsupported display normalization mode: {mode}")

    outlier = compute_tail_jump_display_outliers(
        values,
        mode=outlier_mode,
        tail_search_fraction=tail_search_fraction,
        tail_jump_factor=tail_jump_factor,
        max_display_outlier_fraction=max_display_outlier_fraction,
    )
    outlier_mask = np.asarray(outlier["sld_display_is_outlier"], dtype=bool)
    non_outlier_display = display[~outlier_mask & np.isfinite(display)]
    if non_outlier_display.size:
        outlier["sld_display_color_vmax"] = float(np.max(non_outlier_display))
    return {
        "sld_display": display,
        "sld_display_is_outlier": outlier_mask,
        "sld_display_mode": mode,
        **outlier,
    }


def compute_tail_jump_display_outliers(
    sld_raw: np.ndarray,
    *,
    mode: str,
    tail_search_fraction: float,
    tail_jump_factor: float,
    max_display_outlier_fraction: float,
) -> dict[str, np.ndarray | float | int | str | None]:
    """Mark display-only SLD outliers from sudden jumps in the high-density tail."""

    values = np.asarray(sld_raw, dtype=float)
    finite_mask = np.isfinite(values)
    finite_values = values[finite_mask]
    outliers = np.zeros(values.shape, dtype=bool)
    finite_max = float(np.max(finite_values)) if finite_values.size else None
    base = {
        "sld_display_outlier_mode": mode,
        "sld_tail_search_fraction": float(tail_search_fraction),
        "sld_tail_jump_factor": float(tail_jump_factor),
        "sld_max_display_outlier_fraction": float(max_display_outlier_fraction),
        "sld_display_outlier_threshold": None,
        "sld_display_color_vmax": finite_max,
        "largest_tail_jump_ratio": 1.0,
        "n_sld_display_outliers": 0,
        "fraction_sld_display_outliers": 0.0,
    }
    if mode == "none" or finite_values.size < 2:
        return {"sld_display_is_outlier": outliers, **base}
    if mode != "tail_jump":
        raise ValueError(f"Unsupported display outlier mode: {mode}")

    order = np.argsort(finite_values, kind="mergesort")
    sorted_values = finite_values[order]
    n_finite = sorted_values.size
    n_total = values.size
    n_tail = max(2, int(np.ceil(n_finite * tail_search_fraction)))
    tail_start = max(0, n_finite - n_tail)
    best_index = None
    best_ratio = 1.0
    largest_ratio = 1.0
    for index in range(tail_start, n_finite - 1):
        before = float(sorted_values[index])
        after = float(sorted_values[index + 1])
        if before <= 0.0:
            ratio = float("inf") if after > 0.0 else 1.0
        else:
            ratio = after / before
        largest_ratio = max(largest_ratio, ratio)
        n_after = n_finite - (index + 1)
        if ratio >= tail_jump_factor and n_after / n_total <= max_display_outlier_fraction:
            if ratio > best_ratio:
                best_ratio = ratio
                best_index = index

    base["largest_tail_jump_ratio"] = float(largest_ratio)
    if best_index is None:
        return {"sld_display_is_outlier": outliers, **base}

    threshold = float(sorted_values[best_index])
    outliers = values > threshold
    n_outliers = int(outliers.sum())
    base.update(
        {
            "sld_display_outlier_threshold": threshold,
            "sld_display_color_vmax": threshold,
            "n_sld_display_outliers": n_outliers,
            "fraction_sld_display_outliers": (
                float(n_outliers / n_total) if n_total else 0.0
            ),
        }
    )
    return {"sld_display_is_outlier": outliers, **base}


def _finite_percentile(values: np.ndarray, percentile: float) -> float:
    clean = np.asarray(values, dtype=float)
    clean = clean[np.isfinite(clean)]
    if len(clean) == 0:
        return float("nan")
    return float(np.percentile(clean, percentile))


def _finite_max(values: np.ndarray) -> float:
    clean = np.asarray(values, dtype=float)
    clean = clean[~np.isnan(clean)]
    if len(clean) == 0:
        return float("nan")
    return float(np.max(clean))


def _build_density_report(
    *,
    data: pd.DataFrame,
    coordinates_analysis: np.ndarray,
    requested_k_neighbors: int,
    effective_k_neighbors: int,
    global_local_k_mean: float,
    distance_floor: float,
    distance_floor_fraction: float,
    display_diagnostics: dict[str, np.ndarray | float | int | str | None],
    near_identity_ro_tolerance_rad: float,
    near_duplicate_coordinate_tolerance_rad: float,
) -> DensityReport:
    floored = data[data["sld_was_floored"]]
    sld_unfloored_values = data["sld_unfloored"].to_numpy(dtype=float)
    sld_raw_values = data["sld_raw"].to_numpy(dtype=float)
    max_sld_unfloored = _finite_max(data["sld_unfloored"].to_numpy(dtype=float))
    max_sld_raw = _finite_max(data["sld_raw"].to_numpy(dtype=float))
    p99_sld_unfloored = _finite_percentile(sld_unfloored_values, 99)
    p99_sld_raw = _finite_percentile(sld_raw_values, 99)
    p99_5_sld_raw = _finite_percentile(sld_raw_values, 99.5)
    near_identity_count = _near_identity_ro_count(
        coordinates_analysis,
        tolerance_rad=near_identity_ro_tolerance_rad,
    )
    duplicate_stats = _near_duplicate_coordinate_stats(
        coordinates_analysis,
        tolerance_rad=near_duplicate_coordinate_tolerance_rad,
    )
    warnings = _density_warnings(
        n_points=len(data),
        near_identity_count=near_identity_count,
        duplicate_point_count=duplicate_stats["point_count"],
        largest_duplicate_cluster=duplicate_stats["largest_cluster"],
    )
    floored_rows = []
    for index, row in floored.iterrows():
        floored_rows.append(
            {
                "particle_key": row["particle_key"],
                "row": int(index),
                "sld_local_k_mean": float(row["sld_local_k_mean"]),
                "sld_effective_local_k_mean": float(row["sld_effective_local_k_mean"]),
                "sld_unfloored": float(row["sld_unfloored"]),
                "sld_raw": float(row["sld_raw"]),
                "reason": "local_k_mean_below_distance_floor",
            }
        )

    n_points = len(data)
    n_floored = len(floored)
    return DensityReport(
        n_points=n_points,
        requested_k_neighbors=requested_k_neighbors,
        effective_k_neighbors=effective_k_neighbors,
        global_local_k_mean=global_local_k_mean,
        distance_floor=distance_floor,
        distance_floor_fraction=distance_floor_fraction,
        n_floored_points=n_floored,
        fraction_floored_points=float(n_floored / n_points) if n_points else 0.0,
        n_inf_sld_unfloored=int(np.isinf(sld_unfloored_values).sum()),
        n_inf_sld_raw=int(np.isinf(sld_raw_values).sum()),
        max_sld_unfloored=max_sld_unfloored,
        max_sld_raw=max_sld_raw,
        p99_sld_unfloored=p99_sld_unfloored,
        p99_sld_raw=p99_sld_raw,
        max_over_p99_unfloored=_ratio(max_sld_unfloored, p99_sld_unfloored),
        max_over_p99_raw=_ratio(max_sld_raw, p99_sld_raw),
        floored_particle_keys=tuple(floored["particle_key"].tolist()),
        floored_rows=tuple(floored_rows),
        sld_display_mode=str(display_diagnostics["sld_display_mode"]),
        p99_5_sld_raw=p99_5_sld_raw,
        sld_display_outlier_mode=str(display_diagnostics["sld_display_outlier_mode"]),
        sld_tail_search_fraction=float(display_diagnostics["sld_tail_search_fraction"]),
        sld_tail_jump_factor=float(display_diagnostics["sld_tail_jump_factor"]),
        sld_max_display_outlier_fraction=float(
            display_diagnostics["sld_max_display_outlier_fraction"]
        ),
        sld_display_outlier_threshold=display_diagnostics["sld_display_outlier_threshold"],
        sld_display_color_vmax=display_diagnostics["sld_display_color_vmax"],
        largest_tail_jump_ratio=float(display_diagnostics["largest_tail_jump_ratio"]),
        n_sld_display_outliers=int(display_diagnostics["n_sld_display_outliers"]),
        fraction_sld_display_outliers=float(
            display_diagnostics["fraction_sld_display_outliers"]
        ),
        max_over_display_vmax=_ratio(
            max_sld_raw,
            float(display_diagnostics["sld_display_color_vmax"])
            if display_diagnostics["sld_display_color_vmax"] is not None
            else max_sld_raw,
        ),
        near_identity_ro_tolerance_rad=near_identity_ro_tolerance_rad,
        n_near_identity_ro=near_identity_count,
        fraction_near_identity_ro=float(near_identity_count / n_points) if n_points else 0.0,
        near_duplicate_coordinate_tolerance_rad=near_duplicate_coordinate_tolerance_rad,
        n_near_duplicate_coordinate_points=int(duplicate_stats["point_count"]),
        largest_near_duplicate_coordinate_cluster=int(duplicate_stats["largest_cluster"]),
        warnings=warnings,
    )


def _near_identity_ro_count(
    coordinates_analysis: np.ndarray,
    *,
    tolerance_rad: float,
) -> int:
    if tolerance_rad < 0:
        return 0
    norms = np.linalg.norm(np.asarray(coordinates_analysis, dtype=float), axis=1)
    return int((norms <= tolerance_rad).sum())


def _near_duplicate_coordinate_stats(
    coordinates_analysis: np.ndarray,
    *,
    tolerance_rad: float,
) -> dict[str, int]:
    coords = np.asarray(coordinates_analysis, dtype=float)
    if len(coords) == 0 or tolerance_rad <= 0:
        return {"point_count": 0, "largest_cluster": 0}
    quantized = np.round(coords / tolerance_rad).astype(np.int64)
    _, counts = np.unique(quantized, axis=0, return_counts=True)
    duplicate_counts = counts[counts > 1]
    if duplicate_counts.size == 0:
        return {"point_count": 0, "largest_cluster": 0}
    return {
        "point_count": int(duplicate_counts.sum()),
        "largest_cluster": int(duplicate_counts.max()),
    }


def _density_warnings(
    *,
    n_points: int,
    near_identity_count: int,
    duplicate_point_count: int,
    largest_duplicate_cluster: int,
) -> tuple[str, ...]:
    if n_points <= 0:
        return ()
    warnings = []
    if near_identity_count >= 10 and near_identity_count / n_points >= 0.01:
        warnings.append(
            "many_near_identity_relative_orientations: "
            f"{near_identity_count}/{n_points} particles have near-identity RO"
        )
    if duplicate_point_count >= 10 and duplicate_point_count / n_points >= 0.01:
        warnings.append(
            "many_near_duplicate_ro_coordinates: "
            f"{duplicate_point_count}/{n_points} particles are in near-duplicate "
            f"RO coordinate clusters; largest_cluster={largest_duplicate_cluster}"
        )
    return tuple(warnings)


def compute_landscape_density(
    ro_result: ROResult,
    *,
    policy: DensityPolicy | None = None,
) -> Landscape:
    """Compute a density-bearing landscape from RO-derived analysis coordinates."""

    _require_ro_result(ro_result)
    policy = policy or DensityPolicy()
    _validate_density_policy(policy)

    coords = _analysis_coordinates_from_ro(
        ro_result,
        coordinate_source=policy.coordinate_source,
    )
    sld = compute_sld_values(
        coords,
        k_neighbors=policy.k_neighbors,
        distance_floor_fraction=policy.distance_floor_fraction,
    )
    display = compute_sld_display_values(
        sld["sld_raw"],
        mode=policy.display_normalization_mode,
        outlier_mode=policy.display_outlier_mode,
        tail_search_fraction=policy.tail_search_fraction,
        tail_jump_factor=policy.tail_jump_factor,
        max_display_outlier_fraction=policy.max_display_outlier_fraction,
    )

    data = pd.DataFrame(
        {
            "particle_key": ro_result.data["particle_key"].tolist(),
            "coordinates_analysis": list(coords),
            "coordinates_display": list(coords.copy()),
            "sld_unfloored": sld["sld_unfloored"],
            "sld_raw": sld["sld_raw"],
            "sld_display": display["sld_display"],
            "sld_display_is_outlier": display["sld_display_is_outlier"],
            "sld_was_floored": sld["sld_was_floored"],
            "sld_local_k_mean": sld["sld_local_k_mean"],
            "sld_effective_local_k_mean": sld["sld_effective_local_k_mean"],
            "sld_distance_floor": sld["sld_distance_floor"],
        }
    )
    density_report = _build_density_report(
        data=data,
        coordinates_analysis=coords,
        requested_k_neighbors=policy.k_neighbors,
        effective_k_neighbors=min(policy.k_neighbors, len(data) - 1),
        global_local_k_mean=float(sld["global_local_k_mean"]),
        distance_floor=float(sld["sld_distance_floor"]),
        distance_floor_fraction=policy.distance_floor_fraction,
        display_diagnostics=display,
        near_identity_ro_tolerance_rad=policy.near_identity_ro_tolerance_rad,
        near_duplicate_coordinate_tolerance_rad=policy.near_duplicate_coordinate_tolerance_rad,
    )
    return Landscape(
        data=data,
        canonical_transform=None,
        active_policies={"density_policy": policy},
        density_report=density_report,
    )
