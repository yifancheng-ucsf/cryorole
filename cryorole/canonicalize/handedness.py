"""Handedness policies for canonical transforms."""

from __future__ import annotations

import numpy as np


def enforce_handedness(axes: np.ndarray, *, handedness_rule: str) -> np.ndarray:
    """Return axes satisfying the requested handedness policy."""

    axes = np.asarray(axes, dtype=float).copy()
    if axes.shape != (3, 3):
        raise ValueError(f"axes must be shape (3, 3), got {axes.shape}")
    if handedness_rule != "right_handed":
        raise ValueError(f"Unsupported canonicalization handedness_rule: {handedness_rule}")

    if np.linalg.det(axes) < 0.0:
        axes[2] *= -1.0
    return axes
