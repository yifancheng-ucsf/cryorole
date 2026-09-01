"""Compact array-native pose and relative-orientation models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PoseArrays:
    """Normalized poses for one source domain without per-row metadata objects."""

    particle_key: np.ndarray
    rotation_matrix_active: np.ndarray
    source_row_id: np.ndarray
    source_type: str
    domain_name: str

    def __post_init__(self) -> None:
        keys = np.asarray(self.particle_key).astype(str)
        matrices = np.asarray(self.rotation_matrix_active, dtype=float)
        rows = np.asarray(self.source_row_id, dtype=np.int64)
        n = len(keys)
        if matrices.shape != (n, 3, 3):
            raise ValueError(f"rotation_matrix_active must have shape ({n}, 3, 3)")
        if rows.shape != (n,):
            raise ValueError(f"source_row_id must have shape ({n},)")
        if not np.isfinite(matrices).all():
            raise ValueError("rotation_matrix_active contains non-finite values")
        object.__setattr__(self, "particle_key", keys)
        object.__setattr__(self, "rotation_matrix_active", matrices)
        object.__setattr__(self, "source_row_id", rows)

    @property
    def n_points(self) -> int:
        return int(len(self.particle_key))


@dataclass(frozen=True)
class MatchedPoseArrays:
    """Two normalized pose arrays already placed in matched particle order."""

    particle_key: np.ndarray
    ref_rotation_matrix_active: np.ndarray
    mov_rotation_matrix_active: np.ndarray
    ref_source_row_id: np.ndarray
    mov_source_row_id: np.ndarray

    def __post_init__(self) -> None:
        ref = PoseArrays(
            self.particle_key,
            self.ref_rotation_matrix_active,
            self.ref_source_row_id,
            "ref",
            "ref",
        )
        mov = PoseArrays(
            self.particle_key,
            self.mov_rotation_matrix_active,
            self.mov_source_row_id,
            "mov",
            "mov",
        )
        if ref.n_points == 0:
            raise ValueError("MatchedPoseArrays must contain at least one particle")
        object.__setattr__(self, "particle_key", ref.particle_key)
        object.__setattr__(self, "ref_rotation_matrix_active", ref.rotation_matrix_active)
        object.__setattr__(self, "mov_rotation_matrix_active", mov.rotation_matrix_active)
        object.__setattr__(self, "ref_source_row_id", ref.source_row_id)
        object.__setattr__(self, "mov_source_row_id", mov.source_row_id)

    @property
    def n_points(self) -> int:
        return int(len(self.particle_key))


@dataclass(frozen=True)
class ROArrays:
    """Vectorized relative-orientation truth and derived representations."""

    particle_key: np.ndarray
    rotation_matrix: np.ndarray
    quaternion_xyzw: np.ndarray
    rotation_vector: np.ndarray
    euler_zyx: np.ndarray
    angle_rad: np.ndarray
    ref_source_row_id: np.ndarray
    mov_source_row_id: np.ndarray

    def __post_init__(self) -> None:
        keys = np.asarray(self.particle_key).astype(str)
        n = len(keys)
        shapes = {
            "rotation_matrix": (n, 3, 3),
            "quaternion_xyzw": (n, 4),
            "rotation_vector": (n, 3),
            "euler_zyx": (n, 3),
            "angle_rad": (n,),
            "ref_source_row_id": (n,),
            "mov_source_row_id": (n,),
        }
        for name, shape in shapes.items():
            dtype = np.int64 if name.endswith("source_row_id") else float
            value = np.asarray(getattr(self, name), dtype=dtype)
            if value.shape != shape:
                raise ValueError(f"{name} must have shape {shape}")
            if dtype is float and not np.isfinite(value).all():
                raise ValueError(f"{name} contains non-finite values")
            object.__setattr__(self, name, value)
        object.__setattr__(self, "particle_key", keys)

    @property
    def n_points(self) -> int:
        return int(len(self.particle_key))
