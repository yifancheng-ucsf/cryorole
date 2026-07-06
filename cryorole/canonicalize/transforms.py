"""Canonical transform validation and application."""

from __future__ import annotations

import numpy as np


def transform_from_axes(axes: np.ndarray) -> np.ndarray:
    """Return a right-multiplication transform from row-wise RV output axes.

    Rows are ordered as canonical rotation-vector x/y/z axes. For the default
    ZYX-aware canonical assignment, those rows correspond to gamma/beta/alpha
    display axes, so PC1 appears in Euler alpha via the RV z column.
    """

    axes = np.asarray(axes, dtype=float)
    if axes.shape != (3, 3):
        raise ValueError(f"axes must be shape (3, 3), got {axes.shape}")
    transform = axes.T
    validate_canonical_transform(transform)
    return transform


def validate_canonical_transform(transform: np.ndarray) -> None:
    """Validate that a canonical transform is a proper 3D rotation."""

    transform = np.asarray(transform, dtype=float)
    if transform.shape != (3, 3):
        raise ValueError(f"canonical_transform must be shape (3, 3), got {transform.shape}")
    if not np.isfinite(transform).all():
        raise ValueError("canonical_transform contains non-finite values")
    if not np.allclose(transform.T @ transform, np.eye(3), atol=1e-8):
        raise ValueError("canonical_transform must be orthonormal")
    if not np.isclose(np.linalg.det(transform), 1.0, atol=1e-8):
        raise ValueError("canonical_transform must be right-handed with determinant +1")


def apply_canonical_transform(
    coordinates_analysis: np.ndarray,
    canonical_transform: np.ndarray,
) -> np.ndarray:
    """Apply a canonical transform to analysis coordinates without translating them."""

    coordinates = np.asarray(coordinates_analysis, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(
            "coordinates_analysis for canonicalization must be an array with shape (n, 3)"
        )
    validate_canonical_transform(canonical_transform)
    return coordinates @ canonical_transform
