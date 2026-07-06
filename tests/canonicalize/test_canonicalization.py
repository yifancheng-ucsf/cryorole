from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.transform import Rotation

from cryorole.canonicalize import canonicalize_landscape, canonicalize_landscape_arrays
from cryorole.canonicalize.sign_rules import apply_axis_sign_rule
from cryorole.models.landscape import Landscape
from cryorole.models.landscape_arrays import LandscapeArrays
from cryorole.models.policies import CanonicalizationPolicy


def _make_landscape(
    coordinates_analysis,
    *,
    coordinates_display=None,
    sld_raw=None,
    sld_display=None,
) -> Landscape:
    n_points = len(coordinates_analysis)
    if coordinates_display is None:
        coordinates_display = coordinates_analysis
    if sld_raw is None:
        sld_raw = np.ones(n_points)
    if sld_display is None:
        sld_display = sld_raw
    data = pd.DataFrame(
        {
            "particle_key": [f"p{i}" for i in range(n_points)],
            "coordinates_analysis": [np.asarray(row, dtype=float) for row in coordinates_analysis],
            "coordinates_display": [np.asarray(row, dtype=float) for row in coordinates_display],
            "sld_unfloored": np.asarray(sld_raw, dtype=float) + 0.5,
            "sld_raw": np.asarray(sld_raw, dtype=float),
            "sld_display": np.asarray(sld_display, dtype=float),
            "sld_was_floored": [False] * n_points,
            "sld_local_k_mean": np.ones(n_points),
            "sld_effective_local_k_mean": np.ones(n_points),
            "sld_distance_floor": np.full(n_points, 1e-4),
        }
    )
    return Landscape(data=data)


def _make_landscape_arrays(coordinates_analysis, *, sld_raw=None) -> LandscapeArrays:
    n_points = len(coordinates_analysis)
    if sld_raw is None:
        sld_raw = np.ones(n_points)
    coordinates = np.asarray(coordinates_analysis, dtype=float)
    return LandscapeArrays(
        particle_key=np.asarray([f"p{i}" for i in range(n_points)]),
        coordinates_analysis=coordinates,
        coordinates_display=coordinates,
        sld_unfloored=np.asarray(sld_raw, dtype=float) + 0.5,
        sld_raw=np.asarray(sld_raw, dtype=float),
        sld_display=np.asarray(sld_raw, dtype=float),
        sld_was_floored=np.zeros(n_points, dtype=bool),
        sld_local_k_mean=np.ones(n_points),
        sld_effective_local_k_mean=np.ones(n_points),
        sld_distance_floor=np.full(n_points, 1e-4),
    )


def test_canonicalization_consumes_analysis_not_display_coordinates():
    analysis = np.array(
        [
            [-2.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ]
    )
    display = np.array(
        [
            [0.0, -2.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 2.0, 0.0],
        ]
    )
    result = canonicalize_landscape(
        _make_landscape(analysis, coordinates_display=display),
        policy=CanonicalizationPolicy(fit_subset="all"),
    )

    assert np.allclose(np.abs(result.canonical_transform[:, 2]), [1.0, 0.0, 0.0])


def test_array_native_canonicalization_matches_dataframe_transform():
    coordinates = np.array(
        [
            [-2.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 0.5, 0.0],
        ]
    )
    policy = CanonicalizationPolicy(fit_subset="all")

    dataframe_result = canonicalize_landscape(_make_landscape(coordinates), policy=policy)
    arrays_result, report = canonicalize_landscape_arrays(
        _make_landscape_arrays(coordinates),
        policy=policy,
    )

    np.testing.assert_allclose(
        arrays_result.canonical_transform,
        dataframe_result.canonical_transform,
    )
    np.testing.assert_allclose(
        arrays_result.coordinates_canonical,
        np.vstack(dataframe_result.data["coordinates_canonical"]),
    )
    assert report.n_fit_points == len(coordinates)


def test_array_native_canonicalization_large_smoke_keeps_ndarray_coordinates():
    n_points = 100_000
    x = np.linspace(-1.0, 1.0, n_points)
    coordinates = np.column_stack((x, 0.1 * x**2, 0.01 * x))

    arrays_result, report = canonicalize_landscape_arrays(
        _make_landscape_arrays(coordinates),
        policy=CanonicalizationPolicy(fit_subset="all"),
    )

    assert arrays_result.coordinates_analysis.shape == (n_points, 3)
    assert arrays_result.coordinates_canonical.shape == (n_points, 3)
    assert isinstance(arrays_result.coordinates_canonical, np.ndarray)
    assert report.n_fit_points == n_points


def test_default_fit_subset_uses_sld_raw_not_sld_display():
    coordinates = np.array(
        [
            [-2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, -5.0, 0.0],
            [0.0, 5.0, 0.0],
            [0.0, 0.0, 4.0],
        ]
    )
    sld_raw = np.array([10.0, 10.0, 10.0, 1.0, 1.0, 1.0])
    sld_display = np.array([1.0, 1.0, 1.0, 100.0, 100.0, 100.0])

    result = canonicalize_landscape(
        _make_landscape(coordinates, sld_raw=sld_raw, sld_display=sld_display)
    )

    assert np.allclose(np.abs(result.canonical_transform[:, 2]), [1.0, 0.0, 0.0])


def test_sld_display_cannot_drive_scientific_canonicalization():
    landscape = _make_landscape(
        np.array(
            [
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ]
        )
    )

    with pytest.raises(ValueError, match="sld_display is display-only"):
        canonicalize_landscape(
            landscape,
            policy=CanonicalizationPolicy(density_support_field="sld_display"),
        )


def test_canonicalization_does_not_overwrite_raw_coordinates_or_sld_values():
    coordinates = np.array(
        [
            [-2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
        ]
    )
    landscape = _make_landscape(coordinates, sld_raw=np.array([3.0, 4.0, 5.0]))
    original_analysis = [row.copy() for row in landscape.data["coordinates_analysis"]]
    original_sld = landscape.data[["sld_unfloored", "sld_raw", "sld_display"]].copy()

    result = canonicalize_landscape(
        landscape,
        policy=CanonicalizationPolicy(fit_subset="all"),
    )

    assert "coordinates_canonical" in result.data.columns
    for before, after in zip(original_analysis, result.data["coordinates_analysis"]):
        assert np.allclose(before, after)
    pd.testing.assert_frame_equal(
        original_sld,
        result.data[["sld_unfloored", "sld_raw", "sld_display"]],
    )
    assert len(result.data) == len(landscape.data)


def test_pca_svd_axis_fit_is_deterministic_on_synthetic_data():
    coordinates = np.array(
        [
            [-3.0, 0.1, 0.0],
            [-1.0, 0.0, 0.0],
            [1.0, -0.1, 0.0],
            [3.0, 0.0, 0.0],
        ]
    )
    landscape = _make_landscape(coordinates)
    policy = CanonicalizationPolicy(fit_subset="all")

    first = canonicalize_landscape(landscape, policy=policy)
    second = canonicalize_landscape(landscape, policy=policy)

    assert np.allclose(first.canonical_transform, second.canonical_transform)
    assert np.allclose(
        np.vstack(first.data["coordinates_canonical"]),
        np.vstack(second.data["coordinates_canonical"]),
    )


def test_largest_component_sign_rule_is_deterministic():
    axes = np.array(
        [
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
        ]
    )

    signed = apply_axis_sign_rule(axes, sign_rule="largest_component_positive")

    assert np.allclose(signed, np.eye(3))


def test_default_canonical_policy_uses_density_weighted_skewness():
    policy = CanonicalizationPolicy()

    assert policy.axis_assignment == "pc123_to_alpha_beta_gamma"
    assert policy.sign_rule == "density_weighted_skewness"
    assert policy.positive_side == "high_density_skew"
    assert policy.sign_weight_field == "sld_raw"


def test_axis_assignment_maps_pc123_to_zyx_euler_alpha_beta_gamma_spread_order():
    coordinates = np.array(
        [
            [-0.30, 0.0, 0.0],
            [-0.10, 0.0, 0.0],
            [0.0, -0.20, 0.0],
            [0.0, 0.20, 0.0],
            [0.10, 0.0, 0.0],
            [0.30, 0.0, 0.0],
            [0.0, 0.0, -0.05],
            [0.0, 0.0, 0.05],
        ]
    )

    result = canonicalize_landscape(
        _make_landscape(coordinates),
        policy=CanonicalizationPolicy(fit_subset="all"),
    )
    canonical = np.vstack(result.data["coordinates_canonical"])
    rv_spread = np.var(canonical, axis=0)
    euler = Rotation.from_rotvec(canonical).as_euler("zyx", degrees=True)
    euler_spread = np.var(euler, axis=0)
    report = result.canonicalization_report

    assert rv_spread[2] > rv_spread[1] > rv_spread[0]
    assert euler_spread[0] > euler_spread[1] > euler_spread[2]
    assert report.axis_assignment == "pc123_to_alpha_beta_gamma"
    assert report.pca_axis_order == (
        "PC1_largest_variance",
        "PC2_second_largest_variance",
        "PC3_third_variance",
    )
    assert report.assigned_coordinate_names == ("alpha", "beta", "gamma")
    assert report.assigned_rotvec_columns == ("z", "y", "x")
    assert report.transform_direction == "canonical_rv = raw_rv @ canonical_transform"


def test_density_weighted_skewness_positive_side_is_policy_controlled():
    coordinates = np.array(
        [
            [-3.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, -2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.0, 1.0],
        ]
    )
    sld_raw = np.array([20.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    high = canonicalize_landscape(
        _make_landscape(coordinates, sld_raw=sld_raw),
        policy=CanonicalizationPolicy(fit_subset="all", positive_side="high_density_skew"),
    )
    low = canonicalize_landscape(
        _make_landscape(coordinates, sld_raw=sld_raw),
        policy=CanonicalizationPolicy(fit_subset="all", positive_side="low_density_skew"),
    )

    high_report = high.canonicalization_report
    low_report = low.canonicalization_report
    assert high_report.axis_weighted_skewness[0] > 0.0
    assert low_report.axis_weighted_skewness[0] < 0.0
    assert high_report.positive_side == "high_density_skew"
    assert low_report.positive_side == "low_density_skew"
    assert high_report.sign_weight_field == "sld_raw"
    assert not np.allclose(high.canonical_transform, low.canonical_transform)


def test_density_weighted_skewness_uses_sld_raw_not_sld_display():
    coordinates = np.array(
        [
            [-3.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, -2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, -1.0],
            [0.0, 0.0, 1.0],
        ]
    )
    sld_raw = np.array([20.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    sld_display = np.array([1.0, 1.0, 1.0, 1.0, 20.0, 1.0, 1.0, 1.0, 1.0])

    result = canonicalize_landscape(
        _make_landscape(coordinates, sld_raw=sld_raw, sld_display=sld_display),
        policy=CanonicalizationPolicy(fit_subset="all"),
    )

    assert result.canonicalization_report.sign_weight_field == "sld_raw"
    assert result.canonicalization_report.axis_weighted_skewness[0] > 0.0


def test_canonical_transform_is_right_handed():
    landscape = _make_landscape(
        np.array(
            [
                [-2.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
    )

    result = canonicalize_landscape(
        landscape,
        policy=CanonicalizationPolicy(fit_subset="all"),
    )

    assert np.isclose(np.linalg.det(result.canonical_transform), 1.0)


def test_origin_policy_preserves_origin_without_recentering():
    coordinates = np.array(
        [
            [10.0, -2.0, 0.0],
            [11.0, -2.0, 0.0],
            [12.0, -2.0, 0.0],
            [13.0, -2.0, 0.0],
        ]
    )

    result = canonicalize_landscape(
        _make_landscape(coordinates),
        policy=CanonicalizationPolicy(fit_subset="all", origin_policy="preserve"),
    )

    canonical = np.vstack(result.data["coordinates_canonical"])
    assert np.allclose(canonical, coordinates @ result.canonical_transform)
    assert np.allclose(canonical.mean(axis=0), coordinates.mean(axis=0) @ result.canonical_transform)
    assert not np.allclose(canonical.mean(axis=0), np.zeros(3))


def test_canonicalization_policy_is_recorded():
    landscape = _make_landscape(
        np.array(
            [
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
            ]
        )
    )
    policy = CanonicalizationPolicy(fit_subset="all")

    result = canonicalize_landscape(landscape, policy=policy)

    assert result.active_policies["canonicalization_policy"] == policy


def test_landscape_rejects_non_orthonormal_canonical_transform():
    landscape = _make_landscape(
        np.array(
            [
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
            ]
        )
    )

    with pytest.raises(ValueError, match="canonical_transform must be orthonormal"):
        Landscape(
            data=landscape.data,
            canonical_transform=np.diag([1.0, 2.0, 1.0]),
        )


def test_landscape_rejects_canonical_transform_with_determinant_not_plus_one():
    landscape = _make_landscape(
        np.array(
            [
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
            ]
        )
    )

    with pytest.raises(ValueError, match="canonical_transform determinant must be \\+1"):
        Landscape(
            data=landscape.data,
            canonical_transform=np.diag([-1.0, 1.0, 1.0]),
        )


def test_canonicalization_diagnostics_include_singular_values_and_n_fit_points():
    coordinates = np.array(
        [
            [-2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, -5.0, 0.0],
            [0.0, 5.0, 0.0],
            [0.0, 0.0, 4.0],
        ]
    )
    landscape = _make_landscape(
        coordinates,
        sld_raw=np.array([10.0, 10.0, 10.0, 1.0, 1.0, 1.0]),
    )

    result = canonicalize_landscape(landscape)
    report = result.canonicalization_report

    assert report is not None
    assert len(report.singular_values) == 3
    assert len(report.explained_variance_ratios) == 3
    assert report.n_fit_points == 3
    assert report.fit_subset == "top_fraction"
    assert report.density_support_field == "sld_raw"


def test_near_degenerate_pca_axes_are_reported_not_failed():
    coordinates = np.array(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ]
    )

    result = canonicalize_landscape(
        _make_landscape(coordinates),
        policy=CanonicalizationPolicy(fit_subset="all"),
    )
    report = result.canonicalization_report

    assert report.axis_degeneracy_detected is True
    assert any(
        warning.startswith("near_degenerate_pca_axes")
        for warning in report.warnings
    )


def test_ambiguous_density_weighted_signs_are_reported():
    coordinates = np.array(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ]
    )

    result = canonicalize_landscape(
        _make_landscape(coordinates),
        policy=CanonicalizationPolicy(fit_subset="all"),
    )
    report = result.canonicalization_report

    assert report.ambiguous_sign_axes
    assert any("axis_sign_ambiguous" in warning for warning in report.warnings)
