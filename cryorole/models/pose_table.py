"""Normalized pose table model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


REQUIRED_POSE_COLUMNS = (
    "particle_key",
    "domain_name",
    "rotation_matrix_active",
    "shift_xy",
    "source_type",
    "source_row_id",
    "raw_metadata",
)


@dataclass
class PoseTable:
    """A normalized per-domain pose table using active matrices as truth."""

    data: pd.DataFrame

    def __post_init__(self) -> None:
        missing = [col for col in REQUIRED_POSE_COLUMNS if col not in self.data.columns]
        if missing:
            raise ValueError(f"PoseTable missing required columns: {missing}")

        for idx, matrix in self.data["rotation_matrix_active"].items():
            arr = np.asarray(matrix, dtype=float)
            if arr.shape != (3, 3):
                raise ValueError(
                    f"rotation_matrix_active at row {idx} must be 3x3, got {arr.shape}"
                )
            if not np.isfinite(arr).all():
                raise ValueError(f"rotation_matrix_active at row {idx} has non-finite values")

    def with_particle_keys(self, keys: pd.Series) -> "PoseTable":
        """Return a copy with resolved particle keys attached."""

        if len(keys) != len(self.data):
            raise ValueError("particle key count must match PoseTable row count")
        updated = self.data.copy()
        updated["particle_key"] = list(keys)
        return PoseTable(updated)

