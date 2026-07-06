"""Policy-driven canonicalization of density-bearing landscapes."""

from __future__ import annotations

import math

import numpy as np

from cryorole.canonicalize.axis_fit import fit_axes_pca_svd
from cryorole.canonicalize.transforms import apply_canonical_transform, transform_from_axes
from cryorole.models.canonicalization_report import CanonicalizationReport
from cryorole.models.landscape import Landscape
from cryorole.models.landscape_arrays import LandscapeArrays
from cryorole.models.policies import CanonicalizationPolicy


def canonicalize_landscape(
    landscape: Landscape,
    *,
    policy: CanonicalizationPolicy | None = None,
) -> Landscape:
    """Return a new Landscape with canonical coordinates and transform added.

    PCA centering is used only to estimate canonical axes. With the current
    origin_policy="preserve" behavior, the canonical rotation is applied to
    coordinates_analysis without translating or recentering coordinates.
    """

    policy = policy or CanonicalizationPolicy()
    _validate_policy(policy)
    fit_coordinates, sign_weights = _select_fit_data(landscape, policy)
    axis_fit = fit_axes_pca_svd(
        fit_coordinates,
        policy=policy,
        density_weights=sign_weights,
    )
    canonical_transform = transform_from_axes(axis_fit.axes)

    output = landscape.data.copy(deep=True)
    analysis_coordinates = _stack_coordinate_column(output["coordinates_analysis"])
    canonical_coordinates = apply_canonical_transform(
        analysis_coordinates,
        canonical_transform,
    )
    output["coordinates_canonical"] = [row.copy() for row in canonical_coordinates]

    active_policies = dict(landscape.active_policies or {})
    active_policies["canonicalization_policy"] = policy
    canonicalization_report = CanonicalizationReport(
        fit_subset=policy.fit_subset,
        n_fit_points=len(fit_coordinates),
        density_support_field=policy.density_support_field,
        singular_values=axis_fit.singular_values,
        explained_variance_ratios=axis_fit.explained_variance_ratios,
        axis_degeneracy_detected=axis_fit.axis_degeneracy_detected,
        fit_top_fraction=policy.fit_top_fraction,
        axis_assignment=axis_fit.axis_assignment,
        pca_axis_order=axis_fit.pca_axis_order,
        assigned_coordinate_names=axis_fit.assigned_coordinate_names,
        assigned_rotvec_columns=axis_fit.assigned_rotvec_columns,
        sign_rule=axis_fit.sign_rule,
        positive_side=axis_fit.positive_side,
        sign_weight_field=axis_fit.sign_weight_field,
        axis_weighted_skewness=axis_fit.axis_weighted_skewness,
        flipped_axes=axis_fit.flipped_axes,
        ambiguous_sign_axes=axis_fit.ambiguous_sign_axes,
        handedness_rule=axis_fit.handedness_rule,
        handedness_adjusted_axes=axis_fit.handedness_adjusted_axes,
        transform_direction=axis_fit.transform_direction,
        warnings=axis_fit.warnings,
    )

    return Landscape(
        data=output,
        canonical_transform=canonical_transform,
        active_policies=active_policies,
        density_report=landscape.density_report,
        canonicalization_report=canonicalization_report,
    )


def canonicalize_landscape_arrays(
    arrays: LandscapeArrays,
    *,
    policy: CanonicalizationPolicy | None = None,
) -> tuple[LandscapeArrays, CanonicalizationReport]:
    """Canonicalize a production NPZ landscape without DataFrame object columns."""

    policy = policy or CanonicalizationPolicy()
    _validate_policy(policy)
    fit_coordinates, sign_weights = _select_fit_data_arrays(arrays, policy)
    axis_fit = fit_axes_pca_svd(
        fit_coordinates,
        policy=policy,
        density_weights=sign_weights,
    )
    canonical_transform = transform_from_axes(axis_fit.axes)
    canonical_coordinates = apply_canonical_transform(
        arrays.coordinates_analysis,
        canonical_transform,
    )
    report = CanonicalizationReport(
        fit_subset=policy.fit_subset,
        n_fit_points=len(fit_coordinates),
        density_support_field=policy.density_support_field,
        singular_values=axis_fit.singular_values,
        explained_variance_ratios=axis_fit.explained_variance_ratios,
        axis_degeneracy_detected=axis_fit.axis_degeneracy_detected,
        fit_top_fraction=policy.fit_top_fraction,
        axis_assignment=axis_fit.axis_assignment,
        pca_axis_order=axis_fit.pca_axis_order,
        assigned_coordinate_names=axis_fit.assigned_coordinate_names,
        assigned_rotvec_columns=axis_fit.assigned_rotvec_columns,
        sign_rule=axis_fit.sign_rule,
        positive_side=axis_fit.positive_side,
        sign_weight_field=axis_fit.sign_weight_field,
        axis_weighted_skewness=axis_fit.axis_weighted_skewness,
        flipped_axes=axis_fit.flipped_axes,
        ambiguous_sign_axes=axis_fit.ambiguous_sign_axes,
        handedness_rule=axis_fit.handedness_rule,
        handedness_adjusted_axes=axis_fit.handedness_adjusted_axes,
        transform_direction=axis_fit.transform_direction,
        warnings=axis_fit.warnings,
    )
    return (
        LandscapeArrays(
            particle_key=arrays.particle_key,
            coordinates_analysis=arrays.coordinates_analysis,
            coordinates_display=arrays.coordinates_display,
            sld_unfloored=arrays.sld_unfloored,
            sld_raw=arrays.sld_raw,
            sld_display=arrays.sld_display,
            sld_was_floored=arrays.sld_was_floored,
            sld_local_k_mean=arrays.sld_local_k_mean,
            sld_effective_local_k_mean=arrays.sld_effective_local_k_mean,
            sld_distance_floor=arrays.sld_distance_floor,
            ref_source_row_id=arrays.ref_source_row_id,
            mov_source_row_id=arrays.mov_source_row_id,
            coordinates_canonical=canonical_coordinates,
            canonical_transform=canonical_transform,
        ),
        report,
    )


def _validate_policy(policy: CanonicalizationPolicy) -> None:
    if policy.method != "motion_axis_pca":
        raise ValueError(f"Unsupported canonicalization method: {policy.method}")
    if policy.axis_assignment != "pc123_to_alpha_beta_gamma":
        raise ValueError(
            f"Unsupported canonicalization axis_assignment: {policy.axis_assignment}"
        )
    if policy.fit_subset not in {"all", "top_fraction"}:
        raise ValueError(f"Unsupported canonicalization fit_subset: {policy.fit_subset}")
    if not 0.0 < policy.fit_top_fraction <= 1.0:
        raise ValueError("fit_top_fraction must be > 0 and <= 1")
    if policy.density_support_field == "sld_display":
        raise ValueError("sld_display is display-only and cannot drive canonicalization fitting")
    if policy.sign_weight_field == "sld_display":
        raise ValueError(
            "sld_display is display-only and cannot drive canonicalization sign disambiguation"
        )
    if policy.sign_rule not in {"density_weighted_skewness", "largest_component_positive"}:
        raise ValueError(f"Unsupported canonicalization sign_rule: {policy.sign_rule}")
    if policy.positive_side not in {"high_density_skew", "low_density_skew"}:
        raise ValueError(
            "positive_side must be 'high_density_skew' or 'low_density_skew'"
        )
    if policy.sign_ambiguity_threshold < 0:
        raise ValueError("sign_ambiguity_threshold must be non-negative")
    if policy.origin_policy != "preserve":
        raise ValueError(f"Unsupported canonicalization origin_policy: {policy.origin_policy}")


def _select_fit_data(
    landscape: Landscape,
    policy: CanonicalizationPolicy,
) -> tuple[np.ndarray, np.ndarray]:
    data = landscape.data
    if policy.density_support_field not in data.columns:
        raise ValueError(
            f"Landscape missing density support field: {policy.density_support_field}"
        )
    if policy.sign_weight_field not in data.columns:
        raise ValueError(f"Landscape missing sign weight field: {policy.sign_weight_field}")

    coordinates = _stack_coordinate_column(data["coordinates_analysis"])
    sign_weights = np.asarray(data[policy.sign_weight_field], dtype=float)
    if not np.isfinite(sign_weights).all():
        raise ValueError(f"{policy.sign_weight_field} contains non-finite values")
    if policy.fit_subset == "all":
        return coordinates, sign_weights

    support_values = np.asarray(data[policy.density_support_field], dtype=float)
    if not np.isfinite(support_values).all():
        raise ValueError(f"{policy.density_support_field} contains non-finite values")

    n_fit = max(2, int(math.ceil(len(data) * policy.fit_top_fraction)))
    n_fit = min(n_fit, len(data))
    selected = np.argsort(support_values, kind="mergesort")[-n_fit:]
    return coordinates[selected], sign_weights[selected]


def _select_fit_data_arrays(
    arrays: LandscapeArrays,
    policy: CanonicalizationPolicy,
) -> tuple[np.ndarray, np.ndarray]:
    if not hasattr(arrays, policy.density_support_field):
        raise ValueError(
            f"Landscape missing density support field: {policy.density_support_field}"
        )
    if not hasattr(arrays, policy.sign_weight_field):
        raise ValueError(f"Landscape missing sign weight field: {policy.sign_weight_field}")
    coordinates = arrays.coordinates_analysis
    sign_weights = np.asarray(getattr(arrays, policy.sign_weight_field), dtype=float)
    if not np.isfinite(sign_weights).all():
        raise ValueError(f"{policy.sign_weight_field} contains non-finite values")
    if policy.fit_subset == "all":
        return coordinates, sign_weights

    support_values = np.asarray(getattr(arrays, policy.density_support_field), dtype=float)
    if not np.isfinite(support_values).all():
        raise ValueError(f"{policy.density_support_field} contains non-finite values")

    n_fit = max(2, int(math.ceil(arrays.n_points * policy.fit_top_fraction)))
    n_fit = min(n_fit, arrays.n_points)
    selected = np.argsort(support_values, kind="mergesort")[-n_fit:]
    return coordinates[selected], sign_weights[selected]


def _stack_coordinate_column(column) -> np.ndarray:
    return np.vstack([np.asarray(coords, dtype=float) for coords in column])
