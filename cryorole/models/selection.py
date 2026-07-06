"""First-class particle selection result model."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

import numpy as np

from cryorole.models.policies import SelectionPolicy


@dataclass(frozen=True)
class Selection:
    """Auditable, non-destructive particle selection from a Landscape."""

    selection_id: str
    parent_landscape_id: str | None
    parent_landscape_metadata: Mapping[str, object] = field(default_factory=dict)
    selection_mode: str = "top_fraction_by_density"
    selection_basis: str = "density:sld_raw"
    metric: str = "so3_geodesic"
    center_input: tuple[float, ...] | None = None
    center_input_representation: str | None = None
    center_input_space: str | None = None
    center_euler_sequence: str | None = None
    center_euler_convention: str | None = None
    center_scipy_euler_sequence: str | None = None
    center_degrees: bool | None = None
    center_evaluated: tuple[float, float, float] | None = None
    evaluation_space: str | None = None
    resolved_evaluation_space: str | None = None
    radius: float | None = None
    radius_unit: str | None = None
    threshold: float | None = None
    sld_min: float | None = None
    sld_max: float | None = None
    threshold_operator: str | None = None
    density_support_field: str | None = "sld_raw"
    density_artifact_policy: str = "include_all"
    density_artifact_flag_field: str | None = None
    density_artifact_candidate_count_before: int | None = None
    density_artifact_candidate_count_after: int | None = None
    density_artifact_excluded_count: int = 0
    top_fraction: float | None = 0.40
    random_fraction: float | None = None
    random_seed: int | None = None
    random_candidate_count: int | None = None
    metadata_domain: str | None = None
    metadata_source_file: str | None = None
    metadata_column: str | None = None
    metadata_values: tuple[str, ...] = ()
    metadata_source_row_id_field: str | None = None
    metadata_candidate_count: int | None = None
    metadata_missing_count: int | None = None
    metadata_invalid_count: int | None = None
    tie_break_rule: str | None = None
    range_coordinate_source: str | None = None
    resolved_range_coordinate_source: str | None = None
    range_representation: str | None = None
    range_euler_sequence: str | None = None
    range_euler_convention: str | None = None
    range_scipy_euler_sequence: str | None = None
    range_degrees: bool | None = None
    range_bounds: Mapping[str, tuple[float | None, float | None] | None] | None = None
    selected_particle_keys: tuple[object, ...] = ()
    selected_count: int = 0
    total_count: int = 0
    active_policy: SelectionPolicy | None = None

    def __post_init__(self) -> None:
        if not self.selection_id:
            raise ValueError("Selection selection_id must be non-empty")
        if self.evaluation_space is not None and self.evaluation_space not in {
            "analysis",
            "canonical",
            "canonical_if_available",
        }:
            raise ValueError(f"Unsupported selection evaluation_space: {self.evaluation_space}")
        if self.resolved_evaluation_space is not None and self.resolved_evaluation_space not in {
            "analysis",
            "canonical",
        }:
            raise ValueError("resolved_evaluation_space must be analysis or canonical")
        if self.center_input_space is not None and self.center_input_space not in {
            "analysis_display",
            "canonical_display",
            "evaluation",
        }:
            raise ValueError(f"Unsupported selection center_input_space: {self.center_input_space}")
        if self.metric not in {
            "so3_geodesic",
            "rotvec_euclidean",
            "euclidean",
            "coordinate_box",
            "random_without_replacement",
            "metadata_value_match",
        }:
            raise ValueError(f"Unsupported selection metric: {self.metric}")
        if self.density_support_field == "sld_display":
            raise ValueError("sld_display is display-only and cannot drive selection")
        if self.density_artifact_policy not in {"include_all", "exclude_display_outliers"}:
            raise ValueError(
                "density_artifact_policy must be include_all or exclude_display_outliers"
            )
        for field_name in (
            "density_artifact_candidate_count_before",
            "density_artifact_candidate_count_after",
            "density_artifact_excluded_count",
        ):
            value = getattr(self, field_name)
            if value is not None and value < 0:
                raise ValueError(f"{field_name} must be non-negative")
        if self.center_evaluated is not None:
            center = np.asarray(self.center_evaluated, dtype=float)
            if center.shape != (3,) or not np.isfinite(center).all():
                raise ValueError("Selection center_evaluated must be a finite length-3 vector")
        if self.radius is not None and (not np.isfinite(self.radius) or self.radius <= 0.0):
            raise ValueError("Selection radius must be a positive finite value")
        if self.radius_unit is not None and self.radius_unit not in {"radians", "degrees"}:
            raise ValueError("Selection radius_unit must be radians or degrees")
        if self.top_fraction is not None and not 0.0 < self.top_fraction <= 1.0:
            raise ValueError("Selection top_fraction must be > 0 and <= 1")
        if self.random_fraction is not None and not 0.0 < self.random_fraction <= 1.0:
            raise ValueError("Selection random_fraction must be > 0 and <= 1")
        if self.random_candidate_count is not None and self.random_candidate_count < 0:
            raise ValueError("Selection random_candidate_count must be non-negative")
        if self.metadata_domain is not None and self.metadata_domain not in {"ref", "mov"}:
            raise ValueError("Selection metadata_domain must be ref or mov")
        if self.metadata_candidate_count is not None and self.metadata_candidate_count < 0:
            raise ValueError("Selection metadata_candidate_count must be non-negative")
        if self.metadata_missing_count is not None and self.metadata_missing_count < 0:
            raise ValueError("Selection metadata_missing_count must be non-negative")
        if self.metadata_invalid_count is not None and self.metadata_invalid_count < 0:
            raise ValueError("Selection metadata_invalid_count must be non-negative")
        if self.threshold is not None and (not np.isfinite(self.threshold) or self.threshold < 0.0):
            raise ValueError("Selection threshold must be a finite non-negative value")
        if self.sld_min is not None and not np.isfinite(self.sld_min):
            raise ValueError("Selection sld_min must be finite")
        if self.sld_max is not None and not np.isfinite(self.sld_max):
            raise ValueError("Selection sld_max must be finite")
        if self.sld_min is not None and self.sld_max is not None and self.sld_min > self.sld_max:
            raise ValueError("Selection sld_min must be <= sld_max")
        if self.threshold_operator is not None and self.threshold_operator != ">=":
            raise ValueError(f"Unsupported selection threshold_operator: {self.threshold_operator}")
        if self.range_coordinate_source is not None and self.range_coordinate_source not in {
            "analysis",
            "canonical",
            "canonical_if_available",
        }:
            raise ValueError(
                f"Unsupported selection range_coordinate_source: {self.range_coordinate_source}"
            )
        if self.resolved_range_coordinate_source is not None and self.resolved_range_coordinate_source not in {
            "analysis",
            "canonical",
        }:
            raise ValueError(
                "resolved_range_coordinate_source must be analysis or canonical"
            )
        if self.range_representation is not None and self.range_representation not in {
            "euler",
            "rotvec",
        }:
            raise ValueError(f"Unsupported selection range_representation: {self.range_representation}")
        if self.selected_count != len(self.selected_particle_keys):
            raise ValueError("Selection selected_count must equal selected_particle_keys length")
        if self.total_count < self.selected_count:
            raise ValueError("Selection total_count must be >= selected_count")
