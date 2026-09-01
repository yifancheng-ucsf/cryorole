from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cryorole.export import landscape_from_arrays
from cryorole.models.landscape_arrays import LandscapeArrays
from cryorole.models.policies import SelectionPolicy
from cryorole.select.selectors import select_particles
from cryorole.select.service import _selection_landscape_from_arrays


def _arrays(rows: int = 20) -> LandscapeArrays:
    coordinates = np.column_stack(
        (np.linspace(-0.2, 0.2, rows), np.zeros(rows), np.zeros(rows))
    )
    density = np.linspace(0.0, 10.0, rows)
    return LandscapeArrays(
        particle_key=np.asarray([f"p{index}" for index in range(rows)]),
        coordinates_analysis=coordinates,
        coordinates_display=coordinates,
        coordinates_canonical=coordinates,
        canonical_transform=np.eye(3),
        sld_unfloored=density,
        sld_raw=density,
        sld_display=density,
        sld_display_is_outlier=np.zeros(rows, dtype=bool),
        sld_was_floored=np.zeros(rows, dtype=bool),
        sld_local_k_mean=np.ones(rows),
        sld_effective_local_k_mean=np.ones(rows),
        sld_distance_floor=np.full(rows, 1e-4),
        ref_source_row_id=np.arange(rows),
        mov_source_row_id=np.arange(rows),
    )


@pytest.mark.parametrize(
    "policy",
    [
        SelectionPolicy(selection_mode="top_fraction_by_density", top_fraction=0.25),
        SelectionPolicy(selection_mode="threshold_by_density", threshold=6.0),
        SelectionPolicy(
            selection_mode="radius_around_center",
            center_input=(0.0, 0.0, 0.0),
            center_input_representation="rotvec",
            evaluation_space="analysis",
            radius=0.08,
            radius_unit="radians",
        ),
        SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_coordinate_source="canonical",
            range_representation="rotvec",
            range_bounds={"x": (-0.05, 0.1)},
        ),
        SelectionPolicy(selection_mode="random", random_fraction=0.3, random_seed=17),
    ],
)
def test_array_adapter_matches_dataframe_selection_membership(policy: SelectionPolicy) -> None:
    arrays = _arrays()
    expected = select_particles(landscape_from_arrays(arrays), policy=policy)
    actual = select_particles(_selection_landscape_from_arrays(arrays, policy), policy=policy)

    assert actual.selected_particle_keys == expected.selected_particle_keys
    assert actual.selected_count == expected.selected_count
    assert actual.total_count == expected.total_count
    assert actual.selection_basis == expected.selection_basis
    assert actual.metric == expected.metric


def test_array_adapter_matches_metadata_selection() -> None:
    arrays = _arrays()
    policy = SelectionPolicy(
        selection_mode="metadata_value",
        metadata_domain="ref",
        metadata_column="class",
        metadata_values=("1",),
        metadata_source_row_id_field="ref_source_row_id",
    )
    metadata = pd.DataFrame({"class": np.arange(arrays.n_points) % 2})

    expected = select_particles(
        landscape_from_arrays(arrays),
        policy=policy,
        source_metadata=metadata,
    )
    actual = select_particles(
        _selection_landscape_from_arrays(arrays, policy),
        policy=policy,
        source_metadata=metadata,
    )

    assert actual.selected_particle_keys == expected.selected_particle_keys
    assert actual.metadata_candidate_count == arrays.n_points
    assert actual.metadata_missing_count == 0
