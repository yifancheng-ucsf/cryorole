"""Shared array-native radius evaluator for CLI and interactive selection."""

from __future__ import annotations

import math

import numpy as np
from scipy.spatial.transform import Rotation


def center_to_rotvec(
    center: np.ndarray | tuple[float, float, float] | list[float],
    *,
    representation: str,
    scipy_euler_sequence: str = "zyx",
    euler_degrees: bool = True,
) -> np.ndarray:
    values = np.asarray(center, dtype=float)
    if values.shape != (3,) or not np.isfinite(values).all():
        raise ValueError(f"{representation} center_input must be a finite length-3 vector")
    if representation == "rotvec":
        return values
    if representation == "euler":
        return Rotation.from_euler(
            scipy_euler_sequence, values, degrees=euler_degrees
        ).as_rotvec()
    raise ValueError(f"Unsupported center representation: {representation}")


def evaluate_radius_mask(
    coordinates: np.ndarray,
    center_rotvec: np.ndarray,
    *,
    radius: float,
    radius_unit: str | None,
    metric: str = "so3_geodesic",
) -> tuple[np.ndarray, np.ndarray]:
    """Return exact full-data mask and distances; never applies display filters."""

    points = np.asarray(coordinates, dtype=float)
    center = np.asarray(center_rotvec, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("selection coordinates must have shape (n, 3) and be finite")
    if center.shape != (3,) or not np.isfinite(center).all():
        raise ValueError("center_rotvec must be a finite length-3 vector")
    threshold = radius_to_radians(radius, radius_unit)
    if metric == "so3_geodesic":
        relative = Rotation.from_rotvec(center).inv() * Rotation.from_rotvec(points)
        distances = relative.magnitude()
    elif metric in {"rotvec_euclidean", "euclidean"}:
        distances = np.linalg.norm(points - center, axis=1)
    else:
        raise ValueError(f"Unsupported radius metric: {metric}")
    return distances <= threshold, distances


def radius_to_radians(radius: float, unit: str | None) -> float:
    value = float(radius)
    if not np.isfinite(value) or value <= 0:
        raise ValueError("selection radius must be a positive finite value")
    if unit == "degrees":
        return math.radians(value)
    if unit in {None, "radians"}:
        return value
    raise ValueError("radius unit must be degrees or radians")
