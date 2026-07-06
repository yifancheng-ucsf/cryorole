"""Axis fitting in RO analysis coordinate space."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cryorole.canonicalize.sign_rules import resolve_axis_signs
from cryorole.canonicalize.sign_rules import weighted_skewness_for_axes
from cryorole.models.policies import CanonicalizationPolicy


AXIS_NAMES = ("alpha", "beta", "gamma")
ASSIGNED_ROTVEC_COLUMNS = ("z", "y", "x")
ROTVEC_ROW_TO_AXIS_NAME = ("gamma", "beta", "alpha")
PCA_AXIS_ORDER = (
    "PC1_largest_variance",
    "PC2_second_largest_variance",
    "PC3_third_variance",
)


@dataclass(frozen=True)
class AxisFitResult:
    """PCA/SVD axes and diagnostics from canonicalization fitting."""

    axes: np.ndarray
    singular_values: tuple[float, ...]
    explained_variance_ratios: tuple[float, ...]
    axis_degeneracy_detected: bool
    axis_assignment: str
    pca_axis_order: tuple[str, ...]
    assigned_coordinate_names: tuple[str, ...]
    assigned_rotvec_columns: tuple[str, ...]
    sign_rule: str
    positive_side: str
    sign_weight_field: str
    axis_weighted_skewness: tuple[float, ...]
    flipped_axes: tuple[str, ...]
    ambiguous_sign_axes: tuple[str, ...]
    handedness_rule: str
    handedness_adjusted_axes: tuple[str, ...]
    transform_direction: str
    warnings: tuple[str, ...]


def fit_axes_pca_svd(
    coordinates_analysis: np.ndarray,
    *,
    policy: CanonicalizationPolicy,
    density_weights: np.ndarray | None = None,
) -> AxisFitResult:
    """Fit canonical axes from analysis coordinates using policy-controlled SVD.

    The fit coordinates are centered only to estimate PCA axes. The returned
    axes do not encode a translation; origin handling is applied downstream by
    the canonicalization origin policy.
    """

    if policy.axis_method != "pca_svd":
        raise ValueError(f"Unsupported canonicalization axis_method: {policy.axis_method}")
    if policy.axis_assignment != "pc123_to_alpha_beta_gamma":
        raise ValueError(
            f"Unsupported canonicalization axis_assignment: {policy.axis_assignment}"
        )

    coordinates = np.asarray(coordinates_analysis, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(
            "coordinates_analysis for canonicalization must be an array with shape (n, 3)"
        )
    if coordinates.shape[0] < 2:
        raise ValueError("Canonicalization requires at least two fit coordinates")
    if not np.isfinite(coordinates).all():
        raise ValueError("coordinates_analysis contains non-finite values")

    centered = coordinates - coordinates.mean(axis=0)
    if np.allclose(centered, 0.0):
        raise ValueError("Canonicalization axis fitting is undefined for collapsed coordinates")

    _, singular_values, vh = np.linalg.svd(centered, full_matrices=False)
    vh = _complete_3d_axes(vh)
    sign_result = resolve_axis_signs(
        vh,
        sign_rule=policy.sign_rule,
        centered_coordinates=centered,
        density_weights=density_weights,
        positive_side=policy.positive_side,
        axis_names=AXIS_NAMES,
        ambiguity_threshold=policy.sign_ambiguity_threshold,
    )
    signed_axes = sign_result.axes
    assigned_axes = _assign_pc_axes_to_rotvec_rows(signed_axes, policy=policy)
    axes, handedness_adjusted_axes = _enforce_assigned_handedness(
        assigned_axes,
        policy=policy,
    )
    display_order_axes = _rotvec_rows_to_display_order_axes(axes)
    if policy.sign_rule == "density_weighted_skewness":
        axis_weighted_skewness = weighted_skewness_for_axes(
            centered,
            display_order_axes,
            np.asarray(density_weights, dtype=float),
            ambiguity_threshold=policy.sign_ambiguity_threshold,
        )
    else:
        axis_weighted_skewness = sign_result.axis_weighted_skewness
    flipped_axes = tuple(
        AXIS_NAMES[index]
        for index in range(vh.shape[0])
        if np.dot(vh[index], display_order_axes[index]) < 0.0
    )
    variance_like = singular_values**2
    variance_total = float(variance_like.sum())
    if variance_total == 0.0:
        explained_ratios = np.zeros_like(variance_like)
    else:
        explained_ratios = variance_like / variance_total
    degeneracy_detected = _detect_near_degenerate_axes(singular_values)
    warnings = list(sign_result.warnings)
    if degeneracy_detected:
        warnings.append(
            "near_degenerate_pca_axes: adjacent singular values are too close "
            "for PCA axes to be uniquely meaningful",
        )
    if handedness_adjusted_axes:
        warnings.append(
            "handedness_enforcement_adjusted_axis: "
            f"{','.join(handedness_adjusted_axes)} axis sign adjusted to satisfy "
            f"{policy.handedness_rule} handedness"
        )
    return AxisFitResult(
        axes=axes,
        singular_values=tuple(float(value) for value in singular_values),
        explained_variance_ratios=tuple(float(value) for value in explained_ratios),
        axis_degeneracy_detected=degeneracy_detected,
        axis_assignment=policy.axis_assignment,
        pca_axis_order=PCA_AXIS_ORDER,
        assigned_coordinate_names=AXIS_NAMES,
        assigned_rotvec_columns=ASSIGNED_ROTVEC_COLUMNS,
        sign_rule=policy.sign_rule,
        positive_side=policy.positive_side,
        sign_weight_field=policy.sign_weight_field,
        axis_weighted_skewness=tuple(float(value) for value in axis_weighted_skewness),
        flipped_axes=flipped_axes,
        ambiguous_sign_axes=sign_result.ambiguous_sign_axes,
        handedness_rule=policy.handedness_rule,
        handedness_adjusted_axes=handedness_adjusted_axes,
        transform_direction="canonical_rv = raw_rv @ canonical_transform",
        warnings=tuple(warnings),
    )


def _detect_near_degenerate_axes(singular_values: np.ndarray) -> bool:
    nonzero = np.asarray(singular_values, dtype=float)
    nonzero = nonzero[nonzero > 1e-12]
    if nonzero.size < 2:
        return False
    return bool(np.any(np.isclose(nonzero[:-1], nonzero[1:], rtol=1e-3, atol=1e-12)))


def _complete_3d_axes(vh: np.ndarray) -> np.ndarray:
    """Complete thin-SVD right singular vectors to a 3x3 orthonormal basis."""

    axes = np.asarray(vh, dtype=float)
    if axes.shape == (3, 3):
        return axes
    if axes.shape == (2, 3):
        third = np.cross(axes[0], axes[1])
        norm = np.linalg.norm(third)
        if norm == 0.0:
            raise ValueError("Canonicalization axis fitting produced degenerate axes")
        return np.vstack([axes, third / norm])
    raise ValueError(f"SVD axes must have shape (2, 3) or (3, 3), got {axes.shape}")


def _assign_pc_axes_to_rotvec_rows(
    axes: np.ndarray,
    *,
    policy: CanonicalizationPolicy,
) -> np.ndarray:
    """Map PC axes to RV rows so ZYX Euler alpha receives PC1.

    The returned row order is RV x/y/z. For ZYX display, alpha is primarily
    rotation around z, beta around y, and gamma around x near the identity.
    Therefore PC1 -> z, PC2 -> y, and PC3 -> x.
    """

    if policy.axis_assignment != "pc123_to_alpha_beta_gamma":
        raise ValueError(
            f"Unsupported canonicalization axis_assignment: {policy.axis_assignment}"
        )
    return np.vstack([axes[2], axes[1], axes[0]])


def _rotvec_rows_to_display_order_axes(axes: np.ndarray) -> np.ndarray:
    """Return axes in alpha/beta/gamma display order from RV x/y/z row order."""

    return np.vstack([axes[2], axes[1], axes[0]])


def _enforce_assigned_handedness(
    axes: np.ndarray,
    *,
    policy: CanonicalizationPolicy,
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Enforce handedness after ZYX-aware axis assignment.

    Reversing PC order into RV x/y/z rows is an odd permutation. When the
    assigned basis is left-handed, adjust the gamma/RV-x row so PC1->alpha and
    PC2->beta signs remain governed by the skewness rule.
    """

    if policy.handedness_rule != "right_handed":
        raise ValueError(f"Unsupported canonicalization handedness_rule: {policy.handedness_rule}")
    output = np.asarray(axes, dtype=float).copy()
    adjusted_axes: tuple[str, ...] = ()
    if np.linalg.det(output) < 0.0:
        output[0] *= -1.0
        adjusted_axes = (ROTVEC_ROW_TO_AXIS_NAME[0],)
    if not np.isclose(np.linalg.det(output), 1.0, atol=1e-8):
        raise ValueError("Assigned canonical axes must be right-handed with determinant +1")
    return output, adjusted_axes
