"""Deterministic sign policies for canonical axes."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class AxisSignResult:
    """Signed axes plus diagnostics for PCA sign disambiguation."""

    axes: np.ndarray
    axis_weighted_skewness: tuple[float, ...]
    sign_flips: tuple[bool, ...]
    ambiguous_sign_axes: tuple[str, ...]
    warnings: tuple[str, ...]


def apply_axis_sign_rule(axes: np.ndarray, *, sign_rule: str) -> np.ndarray:
    """Return axes with deterministic signs according to policy."""

    axes = np.asarray(axes, dtype=float).copy()
    if axes.shape != (3, 3):
        raise ValueError(f"axes must be shape (3, 3), got {axes.shape}")
    if sign_rule != "largest_component_positive":
        raise ValueError(f"Unsupported canonicalization sign_rule: {sign_rule}")

    for axis_index in range(axes.shape[0]):
        largest_component = int(np.argmax(np.abs(axes[axis_index])))
        if axes[axis_index, largest_component] < 0:
            axes[axis_index] *= -1.0
    return axes


def resolve_axis_signs(
    axes: np.ndarray,
    *,
    sign_rule: str,
    centered_coordinates: np.ndarray | None = None,
    density_weights: np.ndarray | None = None,
    positive_side: str = "high_density_skew",
    axis_names: tuple[str, ...] = ("alpha", "beta", "gamma"),
    ambiguity_threshold: float = 1e-8,
) -> AxisSignResult:
    """Resolve PCA axis signs under an explicit policy and return diagnostics."""

    original_axes = _validate_axes(axes)
    if sign_rule == "largest_component_positive":
        signed_axes = apply_axis_sign_rule(
            original_axes,
            sign_rule="largest_component_positive",
        )
        return AxisSignResult(
            axes=signed_axes,
            axis_weighted_skewness=(),
            sign_flips=_sign_flips(original_axes, signed_axes),
            ambiguous_sign_axes=(),
            warnings=(),
        )
    if sign_rule != "density_weighted_skewness":
        raise ValueError(f"Unsupported canonicalization sign_rule: {sign_rule}")
    if positive_side not in {"high_density_skew", "low_density_skew"}:
        raise ValueError(
            "positive_side must be 'high_density_skew' or 'low_density_skew'"
        )
    if ambiguity_threshold < 0:
        raise ValueError("sign_ambiguity_threshold must be non-negative")

    coordinates = _validate_centered_coordinates(centered_coordinates)
    weights = _validate_density_weights(density_weights, n_points=coordinates.shape[0])
    if tuple(axis_names) != ("alpha", "beta", "gamma"):
        raise ValueError("canonical axis names must be alpha, beta, gamma")

    signed_axes = apply_axis_sign_rule(
        original_axes,
        sign_rule="largest_component_positive",
    )
    skewness = list(
        weighted_skewness_for_axes(
            coordinates,
            signed_axes,
            weights,
            ambiguity_threshold=ambiguity_threshold,
        )
    )
    ambiguous_axes: list[str] = []
    for axis_index, value in enumerate(skewness):
        axis_name = axis_names[axis_index]
        if abs(value) <= ambiguity_threshold:
            ambiguous_axes.append(axis_name)
            continue
        if positive_side == "high_density_skew" and value < 0.0:
            signed_axes[axis_index] *= -1.0
            skewness[axis_index] *= -1.0
        elif positive_side == "low_density_skew" and value > 0.0:
            signed_axes[axis_index] *= -1.0
            skewness[axis_index] *= -1.0

    warnings = ()
    if ambiguous_axes:
        warnings = (
            "axis_sign_ambiguous: density-weighted skewness is near zero for axes "
            + ",".join(ambiguous_axes),
        )
    return AxisSignResult(
        axes=signed_axes,
        axis_weighted_skewness=tuple(float(value) for value in skewness),
        sign_flips=_sign_flips(original_axes, signed_axes),
        ambiguous_sign_axes=tuple(ambiguous_axes),
        warnings=warnings,
    )


def weighted_skewness_for_axes(
    centered_coordinates: np.ndarray,
    axes: np.ndarray,
    weights: np.ndarray,
    *,
    ambiguity_threshold: float = 1e-8,
) -> tuple[float, ...]:
    """Compute density-weighted skewness along row-wise canonical axes."""

    coordinates = _validate_centered_coordinates(centered_coordinates)
    signed_axes = _validate_axes(axes)
    density_weights = _validate_density_weights(weights, n_points=coordinates.shape[0])
    projections = coordinates @ signed_axes.T
    total_weight = float(density_weights.sum())
    values: list[float] = []
    for axis_index in range(projections.shape[1]):
        projected = projections[:, axis_index]
        weighted_mean = float(np.sum(density_weights * projected) / total_weight)
        centered_projection = projected - weighted_mean
        variance = float(
            np.sum(density_weights * centered_projection**2) / total_weight
        )
        if variance <= ambiguity_threshold:
            values.append(0.0)
            continue
        third_moment = float(
            np.sum(density_weights * centered_projection**3) / total_weight
        )
        values.append(third_moment / (variance ** 1.5))
    return tuple(values)


def _validate_axes(axes: np.ndarray) -> np.ndarray:
    axes = np.asarray(axes, dtype=float).copy()
    if axes.shape != (3, 3):
        raise ValueError(f"axes must be shape (3, 3), got {axes.shape}")
    if not np.isfinite(axes).all():
        raise ValueError("axes contain non-finite values")
    return axes


def _validate_centered_coordinates(value: np.ndarray | None) -> np.ndarray:
    if value is None:
        raise ValueError(
            "density_weighted_skewness sign_rule requires centered_coordinates"
        )
    coordinates = np.asarray(value, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("centered_coordinates must have shape (n, 3)")
    if coordinates.shape[0] < 2:
        raise ValueError("density_weighted_skewness requires at least two points")
    if not np.isfinite(coordinates).all():
        raise ValueError("centered_coordinates contains non-finite values")
    return coordinates


def _validate_density_weights(value: np.ndarray | None, *, n_points: int) -> np.ndarray:
    if value is None:
        raise ValueError("density_weighted_skewness sign_rule requires density_weights")
    weights = np.asarray(value, dtype=float)
    if weights.shape != (n_points,):
        raise ValueError(
            f"density_weights must have shape ({n_points},), got {weights.shape}"
        )
    if not np.isfinite(weights).all():
        raise ValueError("density_weights contains non-finite values")
    if (weights < 0.0).any():
        raise ValueError("density_weights must be non-negative")
    if float(weights.sum()) <= 0.0:
        raise ValueError("density_weights must have positive total weight")
    return weights


def _sign_flips(original_axes: np.ndarray, signed_axes: np.ndarray) -> tuple[bool, ...]:
    return tuple(
        bool(np.dot(original_axes[index], signed_axes[index]) < 0.0)
        for index in range(original_axes.shape[0])
    )
