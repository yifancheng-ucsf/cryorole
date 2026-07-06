"""Diagnostics for policy-driven landscape canonicalization."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CanonicalizationReport:
    """Report PCA/SVD fit diagnostics for a canonicalization run."""

    fit_subset: str
    n_fit_points: int
    density_support_field: str
    singular_values: tuple[float, ...]
    explained_variance_ratios: tuple[float, ...]
    axis_degeneracy_detected: bool
    fit_top_fraction: float | None = None
    axis_assignment: str = "pc123_to_alpha_beta_gamma"
    pca_axis_order: tuple[str, ...] = (
        "PC1_largest_variance",
        "PC2_second_largest_variance",
        "PC3_third_variance",
    )
    assigned_coordinate_names: tuple[str, ...] = ("alpha", "beta", "gamma")
    assigned_rotvec_columns: tuple[str, ...] = ("z", "y", "x")
    sign_rule: str = "density_weighted_skewness"
    positive_side: str = "high_density_skew"
    sign_weight_field: str = "sld_raw"
    axis_weighted_skewness: tuple[float, ...] = ()
    flipped_axes: tuple[str, ...] = ()
    ambiguous_sign_axes: tuple[str, ...] = ()
    handedness_rule: str = "right_handed"
    handedness_adjusted_axes: tuple[str, ...] = ()
    transform_direction: str = "canonical_rv = raw_rv @ canonical_transform"
    warnings: tuple[str, ...] = ()
