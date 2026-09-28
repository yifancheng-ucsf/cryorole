import numpy as np
import pandas as pd
import pytest
from scipy.spatial.transform import Rotation

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


def _brute_force_so3_local_k_mean(rotvecs: np.ndarray, k_neighbors: int) -> np.ndarray:
    quaternions = Rotation.from_rotvec(rotvecs).as_quat()
    dots = np.clip(np.abs(quaternions @ quaternions.T), 0.0, 1.0)
    distances = 2.0 * np.arccos(dots)
    np.fill_diagonal(distances, np.inf)
    k = min(k_neighbors, len(rotvecs) - 1)
    return np.mean(np.sort(distances, axis=1)[:, :k], axis=1)


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


def test_default_sld_metric_matches_explicit_rotvec_euclidean() -> None:
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.2, 0.1, 0.0],
            [-0.3, 0.0, 0.1],
            [0.0, 0.4, -0.2],
        ]
    )

    implicit = compute_sld_values(
        coords,
        k_neighbors=2,
        distance_floor_fraction=1e-4,
    )
    explicit = compute_sld_values(
        coords,
        k_neighbors=2,
        distance_floor_fraction=1e-4,
        sld_metric="rotvec_euclidean",
    )

    for field in (
        "sld_local_k_mean",
        "global_local_k_mean",
        "sld_distance_floor",
        "sld_effective_local_k_mean",
        "sld_unfloored",
        "sld_raw",
        "sld_was_floored",
    ):
        np.testing.assert_array_equal(implicit[field], explicit[field])


def test_so3_sld_matches_brute_force_reference_with_small_batches() -> None:
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.3, -0.1, 0.2],
            [-0.4, 0.2, 0.1],
            [0.1, 0.5, -0.2],
            [-0.2, -0.1, 0.6],
        ]
    )

    result = compute_sld_values(
        coords,
        k_neighbors=3,
        distance_floor_fraction=1e-4,
        sld_metric="so3_geodesic",
        query_batch_size=2,
    )
    expected_local = _brute_force_so3_local_k_mean(coords, 3)
    expected_global = float(np.mean(expected_local))

    np.testing.assert_allclose(result["sld_local_k_mean"], expected_local, atol=1e-14)
    assert result["global_local_k_mean"] == pytest.approx(expected_global)
    np.testing.assert_allclose(result["sld_raw"], expected_global / expected_local)


def test_so3_sld_handles_quaternion_double_cover_and_rotvec_pi_boundary() -> None:
    epsilon = 1e-6
    exact_double_cover = np.array(
        [
            [np.pi, 0.0, 0.0],
            [-np.pi, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )
    coords = np.array(
        [
            [np.pi - epsilon, 0.0, 0.0],
            [-np.pi + epsilon, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )

    exact_result = compute_sld_values(
        exact_double_cover,
        k_neighbors=1,
        distance_floor_fraction=1e-4,
        sld_metric="so3_geodesic",
        query_batch_size=1,
    )
    result = compute_sld_values(
        coords,
        k_neighbors=1,
        distance_floor_fraction=1e-4,
        sld_metric="so3_geodesic",
        query_batch_size=1,
    )

    np.testing.assert_allclose(exact_result["sld_local_k_mean"][:2], 0.0, atol=1e-15)
    assert result["sld_local_k_mean"][0] == pytest.approx(2.0 * epsilon, abs=1e-10)
    assert result["sld_local_k_mean"][1] == pytest.approx(2.0 * epsilon, abs=1e-10)


def test_so3_sld_keeps_distinct_particles_with_identical_rotations() -> None:
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.4, 0.0, 0.0],
        ]
    )

    result = compute_sld_values(
        coords,
        k_neighbors=1,
        distance_floor_fraction=1e-4,
        sld_metric="so3_geodesic",
        query_batch_size=1,
    )

    np.testing.assert_allclose(result["sld_local_k_mean"][:2], 0.0, atol=1e-15)
    assert np.isinf(result["sld_unfloored"][0])
    assert np.isinf(result["sld_unfloored"][1])
    assert np.isfinite(result["sld_raw"]).all()


def test_so3_sld_supports_k_equal_to_n_minus_one() -> None:
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.2, 0.0, 0.0],
            [0.0, 0.3, 0.0],
            [0.0, 0.0, 0.4],
        ]
    )

    result = compute_sld_values(
        coords,
        k_neighbors=99,
        distance_floor_fraction=1e-4,
        sld_metric="so3_geodesic",
        query_batch_size=2,
    )

    np.testing.assert_allclose(
        result["sld_local_k_mean"],
        _brute_force_so3_local_k_mean(coords, len(coords) - 1),
    )


def test_compute_sld_values_rejects_unknown_metric_and_nonfinite_coordinates() -> None:
    coords = np.array([[0.0, 0.0, 0.0], [0.2, 0.0, 0.0]])

    with pytest.raises(ValueError, match="Unsupported SLD metric"):
        compute_sld_values(
            coords,
            k_neighbors=1,
            distance_floor_fraction=1e-4,
            sld_metric="unknown",
        )
    with pytest.raises(ValueError, match="non-finite"):
        compute_sld_values(
            np.array([[0.0, 0.0, 0.0], [np.nan, 0.0, 0.0]]),
            k_neighbors=1,
            distance_floor_fraction=1e-4,
            sld_metric="so3_geodesic",
        )


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


def test_density_report_distinguishes_high_sld_from_discontinuous_tail() -> None:
    rotvecs = [np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0])]
    rotvecs.extend(np.array([index / 1000.0, 0.0, 0.0]) for index in range(2, 1000))

    landscape = compute_landscape_density(
        _ro_result_from_rotvecs(rotvecs),
        policy=DensityPolicy(k_neighbors=1),
    )

    report = landscape.density_report
    assert report is not None
    assert report.high_sld_threshold == 100.0
    assert report.n_high_sld_points == 2
    assert report.fraction_high_sld_points == pytest.approx(0.002)
    assert "HIGH_SLD_PRESENT" in report.warning_codes
    assert "DISCONTINUOUS_SLD_TAIL" in report.warning_codes
    assert "DISTANCE_FLOOR_APPLIED" in report.warning_codes
    assert report.sld_display_outlier_threshold is not None
    np.testing.assert_allclose(
        landscape.data["sld_display"],
        landscape.data["sld_raw"],
    )
    assert len(landscape.data) == 1000


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


def test_density_report_records_requested_and_resolved_so3_metric() -> None:
    ro_result = _ro_result_from_rotvecs(
        [
            np.array([0.0, 0.0, 0.0]),
            np.array([0.3, 0.0, 0.0]),
            np.array([0.0, 0.4, 0.0]),
        ]
    )

    landscape = compute_landscape_density(
        ro_result,
        policy=DensityPolicy(k_neighbors=1, sld_metric="so3_geodesic"),
        query_batch_size=1,
    )

    assert landscape.density_report is not None
    assert landscape.density_report.requested_sld_metric == "so3_geodesic"
    assert landscape.density_report.resolved_sld_metric == "so3_geodesic"
    assert landscape.active_policies["density_policy"].sld_metric == "so3_geodesic"


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


def test_rotvec_euclidean_query_batch_size_bounds_each_tree_query(monkeypatch):
    from cryorole.core import density as density_module

    query_sizes: list[int] = []

    class RecordingTree:
        def __init__(self, points):
            self.points = np.asarray(points, dtype=float)

        def query(self, queries, k):
            values = np.asarray(queries, dtype=float)
            query_sizes.append(len(values))
            distances = np.linalg.norm(values[:, None, :] - self.points[None, :, :], axis=2)
            order = np.argsort(distances, axis=1)[:, :k]
            return np.take_along_axis(distances, order, axis=1), order

    monkeypatch.setattr(density_module, "cKDTree", RecordingTree)
    coords = np.arange(30, dtype=float).reshape(10, 3) / 100.0

    compute_sld_values(
        coords,
        k_neighbors=2,
        distance_floor_fraction=1e-4,
        query_batch_size=3,
    )

    assert query_sizes == [3, 3, 3, 1]


@pytest.mark.parametrize("sld_metric", ["rotvec_euclidean", "so3_geodesic"])
@pytest.mark.parametrize("distance_floor_fraction", [1e-4, 0.0])
def test_array_native_density_accepts_inf_sld_for_more_than_k_identical_ros(
    sld_metric, distance_floor_fraction
) -> None:
    """P0-3: >k identical ROs give +inf unfloored SLD, which is reported, not fatal."""

    from cryorole.core.density import compute_landscape_density_arrays

    rng = np.random.default_rng(0)
    coords = np.vstack([np.tile([0.1, 0.2, 0.3], (60, 1)), rng.normal(0.0, 0.3, (100, 3))])
    n = len(coords)
    policy = DensityPolicy(
        k_neighbors=50, sld_metric=sld_metric, distance_floor_fraction=distance_floor_fraction
    )

    arrays, report = compute_landscape_density_arrays(
        np.arange(n).astype(str), coords, np.arange(n), np.arange(n), policy=policy
    )

    assert report.n_inf_sld_unfloored == 60
    assert np.isposinf(arrays.sld_unfloored[:60]).all()
    assert not np.isnan(arrays.sld_raw).any()
    assert report.n_inf_sld_raw == (60 if distance_floor_fraction == 0.0 else 0)


@pytest.mark.parametrize("sld_metric", ["rotvec_euclidean", "so3_geodesic"])
@pytest.mark.parametrize("distance_floor_fraction", [1e-4, 0.0])
def test_identical_ros_with_rounding_noise_still_give_inf_unfloored_sld(
    sld_metric, distance_floor_fraction
) -> None:
    """Identical poses give ROs that differ by rounding only (~1e-17 rad on macOS wheels).

    They must be classified exactly like bitwise-identical ROs, on every platform.
    """

    from cryorole.core.density import compute_landscape_density_arrays

    rng = np.random.default_rng(0)
    near_identity = rng.normal(0.0, 2e-17, (30, 3))
    base = Rotation.from_rotvec([0.1, 0.2, 0.3])
    near_base = (base * Rotation.from_rotvec(rng.normal(0.0, 1e-16, (30, 3)))).as_rotvec()
    coords = np.vstack([near_identity, near_base, rng.normal(0.0, 0.3, (100, 3))])
    assert np.unique(coords[:30], axis=0).shape[0] == 30  # really not bitwise identical
    n = len(coords)
    policy = DensityPolicy(
        k_neighbors=20, sld_metric=sld_metric, distance_floor_fraction=distance_floor_fraction
    )

    arrays, report = compute_landscape_density_arrays(
        np.arange(n).astype(str), coords, np.arange(n), np.arange(n), policy=policy
    )

    assert report.n_inf_sld_unfloored == 60
    assert np.isposinf(arrays.sld_unfloored[:60]).all()
    assert np.isfinite(arrays.sld_unfloored[60:]).all()
    assert report.n_inf_sld_raw == (60 if distance_floor_fraction == 0.0 else 0)


@pytest.mark.parametrize("sld_metric", ["rotvec_euclidean", "so3_geodesic"])
def test_tiny_but_real_distances_are_resolved_and_not_treated_as_coincident(sld_metric) -> None:
    """The coincidence tolerance (1e-12 rad) must not swallow small physical distances."""

    spacing = 1e-9
    coords = np.array([[0.2, 0.0, 0.0], [0.2 + spacing, 0.0, 0.0], [0.5, 0.1, 0.0], [-0.3, 0.2, 0.1]])

    result = compute_sld_values(
        coords, k_neighbors=1, distance_floor_fraction=0.0, sld_metric=sld_metric, query_batch_size=2
    )

    # arccos(|q.p|) could not resolve this (it returns 0 or ~2e-8 rad); the chord form does.
    assert result["sld_local_k_mean"][0] == pytest.approx(spacing, rel=1e-6)
    assert result["sld_local_k_mean"][1] == pytest.approx(spacing, rel=1e-6)
    assert np.isfinite(result["sld_unfloored"]).all()


def test_landscape_arrays_accept_posinf_sld_but_reject_nan_and_neginf() -> None:
    from cryorole.models.landscape_arrays import LandscapeArrays

    def build(sld_unfloored):
        n = 3
        return LandscapeArrays(
            particle_key=np.array(["a", "b", "c"]),
            coordinates_analysis=np.zeros((n, 3)),
            coordinates_display=None,
            sld_unfloored=sld_unfloored,
            sld_raw=np.ones(n),
            sld_display=np.ones(n),
            sld_was_floored=np.zeros(n, dtype=bool),
            sld_local_k_mean=np.ones(n),
            sld_effective_local_k_mean=np.ones(n),
            sld_distance_floor=np.zeros(n),
        )

    assert np.isposinf(build(np.array([1.0, np.inf, 2.0])).sld_unfloored).sum() == 1
    with pytest.raises(ValueError, match="NaN or -inf"):
        build(np.array([1.0, np.nan, 2.0]))
    with pytest.raises(ValueError, match="NaN or -inf"):
        build(np.array([1.0, -np.inf, 2.0]))
