"""Relative-orientation result model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


REQUIRED_RO_COLUMNS = (
    "particle_key",
    "R_ro",
    "quat_ro",
    "rotvec_ro",
    "euler_ro",
    "angle_ro",
    "reference_source_row_id",
    "moving_source_row_id",
)


@dataclass
class ROResult:
    """Per-particle relative orientation results derived from active matrices."""

    data: pd.DataFrame

    def __post_init__(self) -> None:
        missing = [col for col in REQUIRED_RO_COLUMNS if col not in self.data.columns]
        if missing:
            raise ValueError(f"ROResult missing required columns: {missing}")
        for idx, matrix in self.data["R_ro"].items():
            arr = np.asarray(matrix, dtype=float)
            if arr.shape != (3, 3):
                raise ValueError(f"R_ro at row {idx} must be 3x3, got {arr.shape}")
            if not np.isfinite(arr).all():
                raise ValueError(f"R_ro at row {idx} has non-finite values")
