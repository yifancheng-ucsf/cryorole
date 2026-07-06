"""Derived rotation representations.

This module consumes internal active rotation matrices only. It does not parse
source Euler conventions and must not contain passive/active bridge logic.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial.transform import Rotation

from cryorole.core.euler_conventions import resolve_euler_convention


def derive_rotation_representations(
    matrix_active: np.ndarray,
    *,
    euler_sequence: str | None = None,
    degrees: bool = True,
) -> dict[str, np.ndarray | float]:
    """Derive quaternion, rotvec, Euler, and angle from an active matrix."""

    matrix = np.asarray(matrix_active, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError(f"matrix_active must be 3x3, got {matrix.shape}")
    rotation = Rotation.from_matrix(matrix)
    resolved_euler = resolve_euler_convention(scipy_euler_sequence=euler_sequence)
    angle = rotation.magnitude()
    if degrees:
        angle = float(np.degrees(angle))
    return {
        "quat_ro": rotation.as_quat(),
        "rotvec_ro": rotation.as_rotvec(),
        "euler_ro": rotation.as_euler(resolved_euler.scipy_euler_sequence, degrees=degrees),
        "angle_ro": angle,
    }
