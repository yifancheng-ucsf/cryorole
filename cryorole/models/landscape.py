"""Orientation landscape model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from cryorole.models.canonicalization_report import CanonicalizationReport
from cryorole.models.density_report import DensityReport


REQUIRED_LANDSCAPE_COLUMNS = (
    "particle_key",
    "coordinates_analysis",
    "coordinates_display",
    "sld_unfloored",
    "sld_raw",
    "sld_display",
    "sld_was_floored",
    "sld_local_k_mean",
    "sld_effective_local_k_mean",
    "sld_distance_floor",
)
OPTIONAL_LANDSCAPE_COLUMNS_WITH_DEFAULTS = {
    "sld_display_is_outlier": False,
}


@dataclass
class Landscape:
    """Orientation landscape with distinct analysis and display density values."""

    data: pd.DataFrame
    canonical_transform: np.ndarray | None = None
    active_policies: dict[str, object] | None = None
    density_report: DensityReport | None = None
    canonicalization_report: CanonicalizationReport | None = None

    def __post_init__(self) -> None:
        if (
            "sld_display_is_outlier" not in self.data.columns
            and "sld_display_was_clipped" in self.data.columns
        ):
            self.data["sld_display_is_outlier"] = self.data["sld_display_was_clipped"]
        for column, default in OPTIONAL_LANDSCAPE_COLUMNS_WITH_DEFAULTS.items():
            if column not in self.data.columns:
                self.data[column] = [default] * len(self.data)
        missing = [col for col in REQUIRED_LANDSCAPE_COLUMNS if col not in self.data.columns]
        if missing:
            raise ValueError(f"Landscape missing required columns: {missing}")
        if self.data["particle_key"].isna().any():
            raise ValueError("Landscape particle_key must not contain missing values")
        coordinate_columns = ["coordinates_analysis", "coordinates_display"]
        if "coordinates_canonical" in self.data.columns:
            coordinate_columns.append("coordinates_canonical")
        for column in coordinate_columns:
            for idx, coords in self.data[column].items():
                arr = np.asarray(coords, dtype=float)
                if arr.shape != (3,):
                    raise ValueError(
                        f"{column} at row {idx} must be length-3, got shape {arr.shape}"
                    )
        for column in ("sld_unfloored",):
            values = np.asarray(self.data[column], dtype=float)
            if np.isnan(values).any():
                raise ValueError(f"{column} contains NaN values")
        for column in (
            "sld_raw",
            "sld_display",
            "sld_local_k_mean",
            "sld_effective_local_k_mean",
            "sld_distance_floor",
        ):
            values = np.asarray(self.data[column], dtype=float)
            if not np.isfinite(values).all():
                raise ValueError(f"{column} contains non-finite values")
        if not self.data["sld_was_floored"].map(lambda value: isinstance(value, (bool, np.bool_))).all():
            raise ValueError("sld_was_floored must contain boolean values")
        if not self.data["sld_display_is_outlier"].map(
            lambda value: isinstance(value, (bool, np.bool_))
        ).all():
            raise ValueError("sld_display_is_outlier must contain boolean values")
        if self.density_report is not None and self.density_report.n_points != len(self.data):
            raise ValueError("density_report.n_points must equal Landscape row count")
        if self.canonical_transform is not None:
            transform = np.asarray(self.canonical_transform, dtype=float)
            if transform.shape != (3, 3):
                raise ValueError(
                    f"canonical_transform must be shape (3, 3), got {transform.shape}"
                )
            if not np.isfinite(transform).all():
                raise ValueError("canonical_transform contains non-finite values")
            if not np.allclose(transform.T @ transform, np.eye(3), atol=1e-8):
                raise ValueError("canonical_transform must be orthonormal")
            if not np.isclose(np.linalg.det(transform), 1.0, atol=1e-8):
                raise ValueError("canonical_transform determinant must be +1")
