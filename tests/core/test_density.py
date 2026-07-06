import numpy as np
import pandas as pd
import pytest

from cryorole.core import (
    compute_landscape_density,
    compute_sld_display_values,
    compute_sld_values,
)
from cryorole.models.density_report import DensityReport
from cryorole.models.landscape import Landscape
from cryorole.models.policies import DensityPolicy
from cryorole.models.ro_result import ROResult


def _ro_result_from_rotvecs(rotvecs: list[np.ndarray]) -> ROResult:
    return ROResult(
        pd.DataFrame(
            {
                "particle_key": [f"p{i}" for i in range(len(rotvecs))],
                "R_ro": [np.eye(3) for _ in rotvecs],
                "quat_ro": [np.array([0.0, 0.0, 0.0, 1.0]) for _ in rotvecs],
                "rotvec_ro": rotvecs,
                "euler_ro": [np.zeros(3) for _ in rotvecs],
                "angle_ro": [0.0 for _ in rotvecs],
                "reference_source_row_id": list(range(len(rotvecs))),
                "moving_source_row_id": list(range(len(rotvecs))),
            }
        )
    )


def test_sld_formula_on_small_synthetic_point_set() -> None:
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ]
    )

    result = compute_sld_values(
        coords,
        k_neighbors=1,
        distance_floor_fraction=1e-4,
    )

    expected_local = np.array([1.0, 1.0, 2.0])
    expected_global = np.mean(expected_local)
    np.testing.assert_allclose(result["sld_local_k_mean"], expected_local)
    assert result["global_local_k_mean"] == expected_global
    np.testing.assert_allclose(result["sld_unfloored"], expected_global / expected_local)
    np.testing.assert_allclose(result["sld_raw"], expected_global / expected_local)


def test_relative_floor_behavior_and_unfloored_retention() -> None:
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
        ]
    )

    result = compute_sld_values(
        coords,
        k_neighbors=1,
        distance_floor_fraction=1e-4,
    )

    assert np.isinf(result["sld_unfloored"][0])
    assert np.isinf(result["sld_unfloored"][1])
    assert result["sld_was_floored"][0]
    assert result["sld_was_floored"][1]
    assert result["sld_raw"][0] == 1 / 1e-4
    assert result["sld_raw"][1] == 1 / 1e-4
    assert result["sld_effective_local_k_mean"][0] == result["sld_distance_floor"]


def test_floored_points_are_reported_and_not_dropped() -> None:
    ro_result = _ro_result_from_rotvecs(
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, 0.0, 0.0]),
            np.array([10.0, 0.0, 0.0]),
        ]
    )

    landscape = compute_landscape_density(
        ro_result,
        policy=DensityPolicy(k_neighbors=1),
    )

    assert len(landscape.data) == 3
    assert landscape.density_report is not None
    assert landscape.density_report.n_floored_points == 2
    assert landscape.density_report.fraction_floored_points == 2 / 3
    assert landscape.density_report.n_inf_sld_unfloored == 2
    assert landscape.density_report.n_inf_sld_raw == 0
    assert landscape.density_report.floored_particle_keys == ("p0", "p1")
    assert {row["reason"] for row in landscape.density_report.floored_rows} == {
        "local_k_mean_below_distance_floor"
    }


def test_density_output_is_deterministic() -> None:
    ro_result = _ro_result_from_rotvecs(
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.3, 0.0, 0.0]),
            np.array([0.0, 0.4, 0.0]),
        ]
    )
    policy = DensityPolicy(k_neighbors=2, display_normalization_mode="minmax")

    first = compute_landscape_density(ro_result, policy=policy)
    second = compute_landscape_density(ro_result, policy=policy)

    np.testing.assert_allclose(first.data["sld_raw"], second.data["sld_raw"])
    np.testing.assert_allclose(first.data["sld_display"], second.data["sld_display"])
    assert first.density_report == second.density_report


def test_raw_and_display_sld_are_separate_fields() -> None:
    ro_result = _ro_result_from_rotvecs(
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 2.0, 0.0]),
        ]
    )

    landscape = compute_landscape_density(
        ro_result,
        policy=DensityPolicy(k_neighbors=2, display_normalization_mode="minmax"),
    )

    assert "sld_raw" in landscape.data.columns
    assert "sld_display" in landscape.data.columns
    assert "sld_display_is_outlier" in landscape.data.columns
    assert "sld_unfloored" in landscape.data.columns
    assert not np.allclose(landscape.data["sld_raw"], landscape.data["sld_display"])


def test_display_normalization_does_not_change_raw_sld_truth() -> None:
    ro_result = _ro_result_from_rotvecs(
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.8, 0.0, 0.0]),
            np.array([0.0, 0.5, 0.0]),
            np.array([0.0, 0.0, 0.6]),
        ]
    )
    raw_landscape = compute_landscape_density(
        ro_result,
        policy=DensityPolicy(k_neighbors=2, display_normalization_mode="identity"),
    )
    display_landscape = compute_landscape_density(
        ro_result,
        policy=DensityPolicy(k_neighbors=2, display_normalization_mode="minmax"),
    )

    np.testing.assert_allclose(raw_landscape.data["sld_raw"], display_landscape.data["sld_raw"])
    np.testing.assert_allclose(
        raw_landscape.data["sld_unfloored"],
        display_landscape.data["sld_unfloored"],
    )
    assert not np.allclose(raw_landscape.data["sld_display"], display_landscape.data["sld_display"])


def test_sld_display_defaults_to_identity_without_changing_raw() -> None:
    ro_result = _ro_result_from_rotvecs(
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, 0.0, 0.0]),
            np.array([10.0, 0.0, 0.0]),
            np.array([20.0, 0.0, 0.0]),
        ]
    )

    landscape = compute_landscape_density(ro_result, policy=DensityPolicy(k_neighbors=1))

    raw = landscape.data["sld_raw"].to_numpy(dtype=float)
    np.testing.assert_allclose(landscape.data["sld_display"], raw)
    np.testing.assert_allclose(landscape.data["sld_raw"], raw)
    assert landscape.density_report is not None
    assert landscape.density_report.sld_display_mode == "identity"
    assert landscape.density_report.sld_display_outlier_mode == "tail_jump"


def test_sld_display_identity_mode_preserves_legacy_display_values() -> None:
    ro_result = _ro_result_from_rotvecs(
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.8, 0.0, 0.0]),
        ]
    )

    landscape = compute_landscape_density(
        ro_result,
        policy=DensityPolicy(k_neighbors=1, display_normalization_mode="identity"),
    )

    np.testing.assert_allclose(landscape.data["sld_display"], landscape.data["sld_raw"])
    assert not landscape.data["sld_display_is_outlier"].any()


def test_tail_jump_display_outlier_detection_marks_qualifying_tail_jump() -> None:
    sld_raw = np.array([1.0, 2.0, 3.0, 4.0, 25.0])

    result = compute_sld_display_values(
        sld_raw,
        mode="identity",
        outlier_mode="tail_jump",
        tail_search_fraction=1.0,
        tail_jump_factor=5.0,
        max_display_outlier_fraction=0.25,
    )

    np.testing.assert_allclose(result["sld_display"], sld_raw)
    assert result["sld_display_is_outlier"].tolist() == [False, False, False, False, True]
    assert result["sld_display_outlier_threshold"] == 4.0
    assert result["sld_display_color_vmax"] == 4.0
    assert result["n_sld_display_outliers"] == 1


def test_tail_jump_display_outlier_detection_ignores_nonqualifying_cases() -> None:
    no_jump = compute_sld_display_values(
        np.array([1.0, 2.0, 3.0, 4.0, 18.0]),
        mode="identity",
        outlier_mode="tail_jump",
        tail_search_fraction=1.0,
        tail_jump_factor=5.0,
        max_display_outlier_fraction=0.25,
    )
    outside_tail = compute_sld_display_values(
        np.array([1.0, 2.0, 20.0, 21.0, 22.0]),
        mode="identity",
        outlier_mode="tail_jump",
        tail_search_fraction=0.4,
        tail_jump_factor=5.0,
        max_display_outlier_fraction=0.5,
    )
    too_many_after_jump = compute_sld_display_values(
        np.array([1.0, 2.0, 3.0, 20.0, 21.0]),
        mode="identity",
        outlier_mode="tail_jump",
        tail_search_fraction=1.0,
        tail_jump_factor=5.0,
        max_display_outlier_fraction=0.25,
    )

    assert not no_jump["sld_display_is_outlier"].any()
    assert not outside_tail["sld_display_is_outlier"].any()
    assert not too_many_after_jump["sld_display_is_outlier"].any()


def test_density_report_records_near_identity_ro_diagnostics() -> None:
    ro_result = _ro_result_from_rotvecs(
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
        ]
    )

    landscape = compute_landscape_density(ro_result, policy=DensityPolicy(k_neighbors=1))

    assert landscape.density_report is not None
    assert landscape.density_report.near_identity_ro_tolerance_rad == 1e-8
    assert landscape.density_report.n_near_identity_ro == 1
    assert landscape.density_report.fraction_near_identity_ro == 1 / 3


def test_complete_collapse_coordinates_raise_clear_error() -> None:
    coords = np.zeros((3, 3))

    with pytest.raises(ValueError, match="SLD is undefined.*particle identity matching.*pose metadata"):
        compute_sld_values(coords, k_neighbors=1, distance_floor_fraction=1e-4)


def test_percentile_diagnostics_ignore_inf_when_finite_values_exist() -> None:
    ro_result = _ro_result_from_rotvecs(
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, 0.0, 0.0]),
            np.array([10.0, 0.0, 0.0]),
            np.array([20.0, 0.0, 0.0]),
        ]
    )

    landscape = compute_landscape_density(ro_result, policy=DensityPolicy(k_neighbors=1))

    assert landscape.density_report is not None
    assert landscape.density_report.n_inf_sld_unfloored == 2
    assert np.isfinite(landscape.density_report.p99_sld_unfloored)
    assert landscape.density_report.p99_sld_unfloored < np.inf


def test_density_report_distinguishes_requested_and_effective_k() -> None:
    ro_result = _ro_result_from_rotvecs(
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
        ]
    )

    landscape = compute_landscape_density(ro_result, policy=DensityPolicy(k_neighbors=50))

    assert landscape.density_report is not None
    assert landscape.density_report.requested_k_neighbors == 50
    assert landscape.density_report.effective_k_neighbors == 2


def _valid_landscape_data() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "particle_key": ["p0"],
            "coordinates_analysis": [np.zeros(3)],
            "coordinates_display": [np.zeros(3)],
            "sld_unfloored": [1.0],
            "sld_raw": [1.0],
            "sld_display": [1.0],
            "sld_display_is_outlier": [False],
            "sld_was_floored": [False],
            "sld_local_k_mean": [1.0],
            "sld_effective_local_k_mean": [1.0],
            "sld_distance_floor": [1e-4],
        }
    )


def _density_report(n_points: int) -> DensityReport:
    return DensityReport(
        n_points=n_points,
        requested_k_neighbors=1,
        effective_k_neighbors=1,
        global_local_k_mean=1.0,
        distance_floor=1e-4,
        distance_floor_fraction=1e-4,
        n_floored_points=0,
        fraction_floored_points=0.0,
        n_inf_sld_unfloored=0,
        n_inf_sld_raw=0,
        max_sld_unfloored=1.0,
        max_sld_raw=1.0,
        p99_sld_unfloored=1.0,
        p99_sld_raw=1.0,
        max_over_p99_unfloored=1.0,
        max_over_p99_raw=1.0,
        floored_particle_keys=(),
        floored_rows=(),
    )


def test_landscape_rejects_missing_particle_key() -> None:
    data = _valid_landscape_data()
    data.loc[0, "particle_key"] = None

    with pytest.raises(ValueError, match="particle_key"):
        Landscape(data=data)


def test_landscape_checks_density_report_n_points_consistency() -> None:
    with pytest.raises(ValueError, match="density_report.n_points"):
        Landscape(data=_valid_landscape_data(), density_report=_density_report(n_points=2))
