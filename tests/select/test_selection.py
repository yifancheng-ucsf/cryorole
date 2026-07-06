from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.transform import Rotation

from cryorole.models.landscape import Landscape
from cryorole.models.policies import SelectionPolicy
from cryorole.select import select_particles


def _make_landscape(
    coordinates_analysis,
    *,
    coordinates_canonical=None,
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
    if coordinates_canonical is not None:
        data["coordinates_canonical"] = [
            np.asarray(row, dtype=float) for row in coordinates_canonical
        ]
    return Landscape(data=data)


def test_top_fraction_selection_uses_sld_raw_by_default():
    landscape = _make_landscape(
        np.zeros((5, 3)),
        sld_raw=np.array([1.0, 10.0, 2.0, 9.0, 3.0]),
        sld_display=np.array([100.0, 1.0, 90.0, 1.0, 80.0]),
    )

    selection = select_particles(landscape)

    assert selection.density_support_field == "sld_raw"
    assert selection.selected_particle_keys == ("p1", "p3")
    assert selection.selected_count == 2
    assert selection.total_count == 5
    assert selection.density_artifact_policy == "include_all"
    assert selection.density_artifact_excluded_count == 0


def test_top_fraction_selection_includes_display_outlier_rows_by_default():
    landscape = _make_landscape(
        np.zeros((4, 3)),
        sld_raw=np.array([100.0, 9.0, 8.0, 1.0]),
    )
    landscape.data["sld_display_is_outlier"] = [True, False, False, False]

    selection = select_particles(landscape, policy=SelectionPolicy(top_fraction=0.25))

    assert selection.selected_particle_keys == ("p0",)
    assert selection.density_artifact_policy == "include_all"
    assert selection.density_artifact_candidate_count_before == 4
    assert selection.density_artifact_candidate_count_after == 4
    assert selection.density_artifact_excluded_count == 0


def test_top_fraction_selection_can_exclude_display_outlier_rows_explicitly():
    landscape = _make_landscape(
        np.zeros((4, 3)),
        sld_raw=np.array([100.0, 9.0, 8.0, 1.0]),
    )
    landscape.data["sld_display_is_outlier"] = [True, False, False, False]

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            top_fraction=0.25,
            density_artifact_policy="exclude_display_outliers",
        ),
    )

    assert selection.selected_particle_keys == ("p1",)
    assert selection.density_artifact_policy == "exclude_display_outliers"
    assert selection.density_artifact_flag_field == "sld_display_is_outlier"
    assert selection.density_artifact_candidate_count_before == 4
    assert selection.density_artifact_candidate_count_after == 3
    assert selection.density_artifact_excluded_count == 1


def test_top_fraction_selection_uses_top_fraction_and_records_tie_break_rule():
    landscape = _make_landscape(
        np.zeros((4, 3)),
        sld_raw=np.array([2.0, 4.0, 4.0, 1.0]),
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(top_fraction=0.50),
    )

    assert selection.selected_particle_keys == ("p1", "p2")
    assert selection.top_fraction == 0.50
    assert selection.tie_break_rule == "stable_input_order"
    assert selection.threshold is None


def test_random_selection_is_reproducible_and_records_policy():
    landscape = _make_landscape(np.zeros((10, 3)))
    policy = SelectionPolicy(
        selection_mode="random",
        random_fraction=0.30,
        random_seed=123,
        selection_id="random-30",
    )

    first = select_particles(landscape, policy=policy)
    second = select_particles(landscape, policy=policy)

    assert first.selected_particle_keys == second.selected_particle_keys
    assert first.selection_mode == "random"
    assert first.selection_basis == "random_fraction"
    assert first.metric == "random_without_replacement"
    assert first.random_fraction == 0.30
    assert first.random_seed == 123
    assert first.random_candidate_count == 10
    assert first.selected_count == 3
    assert first.total_count == 10


def test_random_selection_different_seed_can_change_keys():
    landscape = _make_landscape(np.zeros((20, 3)))

    first = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="random",
            random_fraction=0.25,
            random_seed=1,
        ),
    )
    second = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="random",
            random_fraction=0.25,
            random_seed=2,
        ),
    )

    assert first.selected_count == second.selected_count == 5
    assert first.selected_particle_keys != second.selected_particle_keys


@pytest.mark.parametrize("random_fraction", [None, 0.0, -0.1, 1.1])
def test_random_selection_requires_valid_fraction(random_fraction):
    landscape = _make_landscape(np.zeros((3, 3)))
    message = "requires random_fraction" if random_fraction is None else "random_fraction"

    with pytest.raises(ValueError, match=message):
        select_particles(
            landscape,
            policy=SelectionPolicy(
                selection_mode="random",
                random_fraction=random_fraction,
                random_seed=0,
            ),
        )


def test_metadata_value_selection_uses_ref_source_row_id():
    landscape = _make_landscape(np.zeros((4, 3)))
    landscape.data["ref_source_row_id"] = [10, 11, 12, 13]
    source_metadata = pd.DataFrame(
        {"_rlnClassNumber": ["1", "2", "1", "3"]},
        index=[10, 11, 12, 13],
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="metadata_value",
            metadata_domain="ref",
            metadata_column="rlnClassNumber",
            metadata_values=("1",),
            metadata_source_file="ref.star",
            metadata_source_row_id_field="ref_source_row_id",
        ),
        source_metadata=source_metadata,
    )

    assert selection.selected_particle_keys == ("p0", "p2")
    assert selection.selection_basis == "source_metadata:ref:rlnClassNumber"
    assert selection.metric == "metadata_value_match"
    assert selection.metadata_domain == "ref"
    assert selection.metadata_source_file == "ref.star"
    assert selection.metadata_column == "rlnClassNumber"
    assert selection.metadata_values == ("1",)
    assert selection.metadata_source_row_id_field == "ref_source_row_id"
    assert selection.metadata_candidate_count == 4
    assert selection.metadata_missing_count == 0


def test_metadata_value_selection_uses_mov_source_row_id():
    landscape = _make_landscape(np.zeros((4, 3)))
    landscape.data["mov_source_row_id"] = [3, 2, 1, 0]
    source_metadata = pd.DataFrame(
        {"_rlnOpticsGroup": ["A", "B", "A", "C"]},
        index=[0, 1, 2, 3],
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="metadata_value",
            metadata_domain="mov",
            metadata_column="rlnOpticsGroup",
            metadata_values=("A",),
            metadata_source_row_id_field="mov_source_row_id",
        ),
        source_metadata=source_metadata,
    )

    assert selection.selected_particle_keys == ("p1", "p3")
    assert selection.metadata_source_row_id_field == "mov_source_row_id"


def test_metadata_multi_value_selection_selects_union():
    landscape = _make_landscape(np.zeros((5, 3)))
    landscape.data["ref_source_row_id"] = [0, 1, 2, 3, 4]
    source_metadata = pd.DataFrame({"_rlnClassNumber": [1, 2, 3, 2, 1]})

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="metadata_value",
            metadata_domain="ref",
            metadata_column="rlnClassNumber",
            metadata_values=("1", "3"),
        ),
        source_metadata=source_metadata,
    )

    assert selection.selected_particle_keys == ("p0", "p2", "p4")


def test_metadata_selection_requires_source_row_field():
    landscape = _make_landscape(np.zeros((2, 3)))
    source_metadata = pd.DataFrame({"_rlnClassNumber": ["1", "2"]})

    with pytest.raises(ValueError, match="ref_source_row_id"):
        select_particles(
            landscape,
            policy=SelectionPolicy(
                selection_mode="metadata_value",
                metadata_domain="ref",
                metadata_column="rlnClassNumber",
                metadata_values=("1",),
            ),
            source_metadata=source_metadata,
        )


def test_top_fraction_selection_does_not_use_threshold():
    landscape = _make_landscape(np.zeros((3, 3)))

    with pytest.raises(ValueError, match="must not use threshold"):
        select_particles(
            landscape,
            policy=SelectionPolicy(threshold=1.0),
        )


def test_threshold_by_density_selects_expected_keys_by_sld_raw():
    landscape = _make_landscape(
        np.zeros((4, 3)),
        sld_raw=np.array([0.5, 4.0, 1.0, 3.0]),
        sld_display=np.array([9.0, 0.1, 9.0, 0.1]),
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="threshold_by_density",
            threshold=3.0,
        ),
    )

    assert selection.selected_particle_keys == ("p1", "p3")
    assert selection.threshold == 3.0
    assert selection.threshold_operator == ">="
    assert selection.top_fraction is None
    assert selection.selected_count == 2
    assert selection.total_count == 4


def test_threshold_by_density_rejects_sld_display():
    landscape = _make_landscape(np.zeros((3, 3)))

    with pytest.raises(ValueError, match="sld_display is display-only"):
        select_particles(
            landscape,
            policy=SelectionPolicy(
                selection_mode="threshold_by_density",
                density_support_field="sld_display",
                threshold=1.0,
            ),
        )


@pytest.mark.parametrize("threshold", [None, -1.0, float("nan"), float("inf")])
def test_threshold_by_density_requires_finite_non_negative_threshold(threshold):
    landscape = _make_landscape(np.zeros((3, 3)))
    message = "requires threshold" if threshold is None else "finite non-negative"

    with pytest.raises(ValueError, match=message):
        select_particles(
            landscape,
            policy=SelectionPolicy(
                selection_mode="threshold_by_density",
                threshold=threshold,
            ),
        )


def test_threshold_by_density_allows_empty_selection_and_records_counts():
    landscape = _make_landscape(
        np.zeros((3, 3)),
        sld_raw=np.array([1.0, 2.0, 3.0]),
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="threshold_by_density",
            threshold=10.0,
            top_fraction=0.0,
        ),
    )

    assert selection.selected_particle_keys == ()
    assert selection.selected_count == 0
    assert selection.total_count == 3
    assert selection.top_fraction is None


def test_radius_selection_works_in_analysis_coordinates_with_rotvec_center():
    landscape = _make_landscape(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ]
        )
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center",
            evaluation_space="analysis",
            center_input=(0.0, 0.0, 0.0),
            center_input_representation="rotvec",
            radius=1.01,
        ),
    )

    assert selection.selected_particle_keys == ("p0", "p1")
    assert selection.center_evaluated == (0.0, 0.0, 0.0)
    assert selection.selection_basis == "coordinates_analysis"
    assert selection.metric == "so3_geodesic"


def test_radius_selection_works_in_canonical_coordinates_when_present():
    landscape = _make_landscape(
        np.array(
            [
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
            ]
        ),
        coordinates_canonical=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.5, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ]
        ),
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center",
            evaluation_space="canonical",
            center_input=(0.0, 0.0, 0.0),
            radius=0.6,
        ),
    )

    assert selection.selected_particle_keys == ("p0", "p1")
    assert selection.selection_basis == "coordinates_canonical"


def test_default_radius_selection_resolves_to_canonical_when_available():
    landscape = _make_landscape(
        np.array(
            [
                [10.0, 0.0, 0.0],
                [11.0, 0.0, 0.0],
                [12.0, 0.0, 0.0],
            ]
        ),
        coordinates_canonical=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.5, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ]
        ),
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center",
            center_input=(0.0, 0.0, 0.0),
            radius=0.6,
        ),
    )

    assert selection.evaluation_space == "canonical_if_available"
    assert selection.resolved_evaluation_space == "canonical"
    assert selection.selection_basis == "coordinates_canonical"
    assert selection.selected_particle_keys == ("p0", "p1")


def test_default_radius_selection_resolves_to_analysis_without_canonical():
    landscape = _make_landscape(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [0.5, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ]
        )
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center",
            center_input=(0.0, 0.0, 0.0),
            radius=0.6,
        ),
    )

    assert selection.evaluation_space == "canonical_if_available"
    assert selection.resolved_evaluation_space == "analysis"
    assert selection.selection_basis == "coordinates_analysis"
    assert selection.selected_particle_keys == ("p0", "p1")


def test_radius_selection_uses_so3_geodesic_distance_by_default():
    near_pi = np.array([np.pi, 0.0, 0.0])
    same_rotation_near_branch = np.array([-np.pi + 0.05, 0.0, 0.0])
    landscape = _make_landscape(
        np.array(
            [
                near_pi,
                same_rotation_near_branch,
                [0.0, 0.0, 0.0],
            ]
        )
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center",
            evaluation_space="analysis",
            center_input=tuple(near_pi),
            radius=0.1,
        ),
    )

    assert selection.metric == "so3_geodesic"
    assert selection.selected_particle_keys == ("p0", "p1")


def test_radius_selection_can_use_explicit_rotvec_euclidean_metric():
    near_pi = np.array([np.pi, 0.0, 0.0])
    same_rotation_near_branch = np.array([-np.pi + 0.05, 0.0, 0.0])
    landscape = _make_landscape(
        np.array(
            [
                near_pi,
                same_rotation_near_branch,
                [0.0, 0.0, 0.0],
            ]
        )
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center",
            evaluation_space="analysis",
            center_input=tuple(near_pi),
            radius=0.1,
            metric="rotvec_euclidean",
        ),
    )

    assert selection.metric == "rotvec_euclidean"
    assert selection.selected_particle_keys == ("p0",)


def test_radius_selection_requires_canonical_coordinates_when_requested():
    landscape = _make_landscape(np.zeros((3, 3)))

    with pytest.raises(ValueError, match="requires coordinates_canonical"):
        select_particles(
            landscape,
            policy=SelectionPolicy(
                selection_mode="radius_around_center",
                evaluation_space="canonical",
                center_input=(0.0, 0.0, 0.0),
                radius=1.0,
            ),
        )


def test_euler_center_input_converts_to_rotvec_and_uses_analysis_coordinates():
    center = Rotation.from_euler("ZYX", (90.0, 0.0, 0.0), degrees=True).as_rotvec()
    landscape = _make_landscape(
        np.array(
            [
                center,
                center + np.array([0.1, 0.0, 0.0]),
                np.array([3.0, 3.0, 3.0]),
            ]
        )
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center",
            evaluation_space="analysis",
            center_input=(90.0, 0.0, 0.0),
            center_input_representation="euler",
            center_euler_sequence="ZYX",
            center_degrees=True,
            radius=0.2,
        ),
    )

    assert np.allclose(selection.center_evaluated, center)
    assert selection.center_input == (90.0, 0.0, 0.0)
    assert selection.center_input_representation == "euler"
    assert selection.center_euler_convention == "intrinsic_zyx"
    assert selection.center_scipy_euler_sequence == "ZYX"
    assert selection.selected_particle_keys == ("p0", "p1")


def test_default_euler_center_uses_extrinsic_fixed_axis_zyx():
    center = Rotation.from_euler("zyx", (90.0, 0.0, 0.0), degrees=True).as_rotvec()
    landscape = _make_landscape(np.array([center, np.array([3.0, 3.0, 3.0])]))

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center",
            evaluation_space="analysis",
            center_input=(90.0, 0.0, 0.0),
            center_input_representation="euler",
            center_degrees=True,
            radius=0.2,
        ),
    )

    assert np.allclose(selection.center_evaluated, center)
    assert selection.center_euler_convention == "extrinsic_zyx"
    assert selection.center_scipy_euler_sequence == "zyx"
    assert selection.selected_particle_keys == ("p0",)


def test_euler_center_uses_canonical_when_available_with_canonical_if_available():
    center = Rotation.from_euler("ZYX", (90.0, 0.0, 0.0), degrees=True).as_rotvec()
    landscape = _make_landscape(
        np.array(
            [
                [20.0, 0.0, 0.0],
                [21.0, 0.0, 0.0],
                [22.0, 0.0, 0.0],
            ]
        ),
        coordinates_canonical=np.array(
            [
                center,
                center + np.array([0.1, 0.0, 0.0]),
                np.array([3.0, 3.0, 3.0]),
            ]
        ),
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center",
            evaluation_space="canonical_if_available",
            center_input=(90.0, 0.0, 0.0),
            center_input_representation="euler",
            center_euler_sequence="ZYX",
            center_degrees=True,
            radius=0.2,
        ),
    )

    assert selection.selection_basis == "coordinates_canonical"
    assert selection.selected_particle_keys == ("p0", "p1")


def test_coordinates_display_is_not_used_for_scientific_radius_selection():
    landscape = _make_landscape(
        np.array(
            [
                [10.0, 0.0, 0.0],
                [11.0, 0.0, 0.0],
                [12.0, 0.0, 0.0],
            ]
        ),
        coordinates_display=np.array(
            [
                [0.0, 0.0, 0.0],
                [0.1, 0.0, 0.0],
                [0.2, 0.0, 0.0],
            ]
        ),
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="radius_around_center",
            evaluation_space="analysis",
            center_input=(0.0, 0.0, 0.0),
            radius=0.5,
        ),
    )

    assert selection.selected_particle_keys == ()


def test_alpha_range_selects_expected_particles_from_canonical_euler_display():
    canonical = Rotation.from_euler(
        "ZYX",
        [
            [45.0, 0.0, 0.0],
            [75.0, 0.0, 0.0],
            [110.0, 0.0, 0.0],
            [140.0, 0.0, 0.0],
        ],
        degrees=True,
    ).as_rotvec()
    landscape = _make_landscape(
        np.zeros((4, 3)),
        coordinates_canonical=canonical,
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_coordinate_source="canonical",
            range_representation="euler",
            range_euler_sequence="ZYX",
            range_degrees=True,
            range_bounds={"alpha": (60.0, 120.0)},
        ),
    )

    assert selection.selected_particle_keys == ("p1", "p2")
    assert selection.selection_basis == "display_coordinates:canonical:euler"
    assert selection.range_coordinate_source == "canonical"
    assert selection.resolved_range_coordinate_source == "canonical"
    assert selection.range_representation == "euler"
    assert selection.range_euler_convention == "intrinsic_zyx"
    assert selection.range_scipy_euler_sequence == "ZYX"
    assert selection.range_bounds == {"alpha": (60.0, 120.0)}


def test_beta_and_gamma_euler_bounds_can_be_combined():
    coordinates = Rotation.from_euler(
        "zyx",
        [
            [0.0, 20.0, 30.0],
            [0.0, 20.0, 70.0],
            [0.0, 45.0, 30.0],
            [0.0, 45.0, 70.0],
        ],
        degrees=True,
    ).as_rotvec()
    landscape = _make_landscape(coordinates)

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_representation="euler",
            range_bounds={"beta": (10.0, 30.0), "gamma": (20.0, 40.0)},
        ),
    )

    assert selection.selected_particle_keys == ("p0",)


def test_unconstrained_range_axes_are_allowed_with_one_constrained_axis():
    coordinates = Rotation.from_euler(
        "zyx",
        [
            [10.0, -20.0, 0.0],
            [70.0, 40.0, 90.0],
            [130.0, -40.0, -90.0],
        ],
        degrees=True,
    ).as_rotvec()
    landscape = _make_landscape(coordinates)

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_representation="euler",
            range_bounds={"alpha": None, "beta": (None, None), "gamma": (-100.0, 100.0)},
        ),
    )

    assert selection.selected_particle_keys == ("p0", "p1", "p2")


def test_wraparound_euler_bounds_work():
    coordinates = Rotation.from_euler(
        "zyx",
        [
            [175.0, 0.0, 0.0],
            [-175.0, 0.0, 0.0],
            [160.0, 0.0, 0.0],
            [-160.0, 0.0, 0.0],
        ],
        degrees=True,
    ).as_rotvec()
    landscape = _make_landscape(coordinates)

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_representation="euler",
            range_bounds={"alpha": (170.0, -170.0)},
        ),
    )

    assert selection.selected_particle_keys == ("p0", "p1")


def test_range_by_coordinates_with_rotvec_works():
    landscape = _make_landscape(
        np.array(
            [
                [0.1, 0.0, 0.0],
                [0.5, 0.0, 0.0],
                [0.9, 0.0, 0.0],
            ]
        )
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_representation="rotvec",
            range_bounds={"x": (0.4, 1.0), "y": (-0.1, 0.1)},
        ),
    )

    assert selection.selected_particle_keys == ("p1", "p2")
    assert selection.selection_basis == "coordinates_analysis:rotvec"
    assert selection.range_euler_sequence is None
    assert selection.range_degrees is None


def test_range_canonical_if_available_resolves_to_canonical_when_present():
    analysis = Rotation.from_euler(
        "zyx",
        [[10.0, 0.0, 0.0], [20.0, 0.0, 0.0], [30.0, 0.0, 0.0]],
        degrees=True,
    ).as_rotvec()
    canonical = Rotation.from_euler(
        "zyx",
        [[80.0, 0.0, 0.0], [90.0, 0.0, 0.0], [130.0, 0.0, 0.0]],
        degrees=True,
    ).as_rotvec()
    landscape = _make_landscape(analysis, coordinates_canonical=canonical)

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_coordinate_source="canonical_if_available",
            range_representation="euler",
            range_bounds={"alpha": (70.0, 100.0)},
        ),
    )

    assert selection.selected_particle_keys == ("p0", "p1")
    assert selection.resolved_range_coordinate_source == "canonical"
    assert selection.selection_basis == "display_coordinates:canonical:euler"


def test_default_range_selection_resolves_to_canonical_when_available():
    analysis = Rotation.from_euler(
        "zyx",
        [[10.0, 0.0, 0.0], [20.0, 0.0, 0.0], [30.0, 0.0, 0.0]],
        degrees=True,
    ).as_rotvec()
    canonical = Rotation.from_euler(
        "zyx",
        [[80.0, 0.0, 0.0], [90.0, 0.0, 0.0], [130.0, 0.0, 0.0]],
        degrees=True,
    ).as_rotvec()
    landscape = _make_landscape(analysis, coordinates_canonical=canonical)

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_representation="euler",
            range_bounds={"alpha": (70.0, 100.0)},
        ),
    )

    assert selection.range_coordinate_source == "canonical_if_available"
    assert selection.resolved_range_coordinate_source == "canonical"
    assert selection.selection_basis == "display_coordinates:canonical:euler"
    assert selection.selected_particle_keys == ("p0", "p1")


def test_default_range_selection_resolves_to_analysis_without_canonical():
    analysis = Rotation.from_euler(
        "zyx",
        [[80.0, 0.0, 0.0], [90.0, 0.0, 0.0], [130.0, 0.0, 0.0]],
        degrees=True,
    ).as_rotvec()
    landscape = _make_landscape(analysis)

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_representation="euler",
            range_bounds={"alpha": (70.0, 100.0)},
        ),
    )

    assert selection.range_coordinate_source == "canonical_if_available"
    assert selection.resolved_range_coordinate_source == "analysis"
    assert selection.selection_basis == "display_coordinates:analysis:euler"
    assert selection.selected_particle_keys == ("p0", "p1")


def test_range_canonical_source_fails_when_canonical_coordinates_are_missing():
    landscape = _make_landscape(np.zeros((3, 3)))

    with pytest.raises(ValueError, match="requires coordinates_canonical"):
        select_particles(
            landscape,
            policy=SelectionPolicy(
                selection_mode="range_by_coordinates",
                range_coordinate_source="canonical",
                range_representation="euler",
                range_bounds={"alpha": (0.0, 10.0)},
            ),
        )


def test_range_selection_preserves_selected_keys_in_input_row_order():
    coordinates = Rotation.from_euler(
        "zyx",
        [
            [90.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [70.0, 0.0, 0.0],
            [110.0, 0.0, 0.0],
        ],
        degrees=True,
    ).as_rotvec()
    landscape = _make_landscape(coordinates)

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_representation="euler",
            range_bounds={"alpha": (60.0, 120.0)},
        ),
    )

    assert selection.selected_particle_keys == ("p0", "p2", "p3")


def test_range_by_coordinates_with_no_constrained_axes_fails_clearly():
    landscape = _make_landscape(np.zeros((3, 3)))

    with pytest.raises(ValueError, match="requires at least one constrained axis"):
        select_particles(
            landscape,
            policy=SelectionPolicy(
                selection_mode="range_by_coordinates",
                range_bounds={"alpha": None, "beta": (None, None)},
            ),
        )


def test_selection_is_non_destructive_and_does_not_modify_landscape():
    landscape = _make_landscape(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ]
        ),
        sld_raw=np.array([3.0, 1.0, 2.0]),
    )
    columns_before = tuple(landscape.data.columns)
    analysis_before = [row.copy() for row in landscape.data["coordinates_analysis"]]
    sld_before = landscape.data[["sld_unfloored", "sld_raw", "sld_display"]].copy()

    select_particles(landscape, policy=SelectionPolicy(top_fraction=0.34))

    assert tuple(landscape.data.columns) == columns_before
    for before, after in zip(analysis_before, landscape.data["coordinates_analysis"]):
        assert np.allclose(before, after)
    pd.testing.assert_frame_equal(
        sld_before,
        landscape.data[["sld_unfloored", "sld_raw", "sld_display"]],
    )


def test_selection_records_active_policy_and_selected_particle_keys():
    landscape = _make_landscape(
        np.zeros((3, 3)),
        sld_raw=np.array([1.0, 5.0, 2.0]),
    )
    policy = SelectionPolicy(
        top_fraction=0.33,
        parent_landscape_id="landscape-1",
        parent_landscape_metadata={"source": "synthetic"},
        selection_id="selection-1",
    )

    selection = select_particles(landscape, policy=policy)

    assert selection.selection_id == "selection-1"
    assert selection.parent_landscape_id == "landscape-1"
    assert selection.parent_landscape_metadata == {"source": "synthetic"}
    assert selection.active_policy == policy
    assert selection.selected_particle_keys == ("p1",)
    assert selection.selected_count == 1
    assert selection.total_count == 3


def test_top_fraction_tie_at_boundary_prefers_earlier_input_rows():
    landscape = _make_landscape(
        np.zeros((4, 3)),
        sld_raw=np.array([10.0, 5.0, 5.0, 1.0]),
    )

    selection = select_particles(
        landscape,
        policy=SelectionPolicy(top_fraction=0.50),
    )

    assert selection.selected_particle_keys == ("p0", "p1")
    assert selection.tie_break_rule == "stable_input_order"


@pytest.mark.parametrize(
    ("policy", "message"),
    [
        (
            SelectionPolicy(
                selection_mode="radius_around_center",
                center_input=(0.0, 0.0, 0.0),
                radius=0.0,
            ),
            "radius must be a positive finite value",
        ),
        (
            SelectionPolicy(
                selection_mode="radius_around_center",
                center_input=(0.0, 0.0),
                radius=1.0,
            ),
            "rotvec center_input must be a finite length-3 vector",
        ),
        (
            SelectionPolicy(evaluation_space="display"),
            "Unsupported selection evaluation_space",
        ),
        (
            SelectionPolicy(center_input_space="analysis"),
            "Unsupported selection center_input_space",
        ),
    ],
)
def test_invalid_radius_center_or_evaluation_space_fails_clearly(policy, message):
    landscape = _make_landscape(np.zeros((3, 3)))

    with pytest.raises(ValueError, match=message):
        select_particles(landscape, policy=policy)
