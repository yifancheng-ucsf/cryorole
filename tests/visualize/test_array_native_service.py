from __future__ import annotations

import numpy as np

from cryorole.export import write_landscape_npz_arrays
from cryorole.models.landscape_arrays import LandscapeArrays
from cryorole.visualize.service import (
    VisualizationRequest,
    _prepare_array_native_visualization,
)


def _arrays(rows: int = 12) -> LandscapeArrays:
    coordinates = np.column_stack(
        (np.linspace(-0.3, 0.3, rows), np.zeros(rows), np.zeros(rows))
    )
    density = np.arange(rows, dtype=float)
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


def test_npz_visualization_filters_full_arrays_before_dataframe_materialization(tmp_path) -> None:
    source = write_landscape_npz_arrays(_arrays(), tmp_path / "landscape.npz")
    request = VisualizationRequest(
        run_dir=str(tmp_path),
        representation="rotvec",
        threshold=5.0,
    )

    landscape, metadata, stats, _scale = _prepare_array_native_visualization(
        request,
        source,
        euler_metadata={"scipy_euler_sequence": "zyx"},
        range_bounds={},
        selected_landscape_report=None,
    )

    assert metadata is None
    assert stats["n_points_full_parent"] == 12
    assert stats["n_points_input"] == 12
    assert stats["n_points_after_display_filter"] == {"analysis": 7}
    assert stats["n_points_2d"] == {"analysis": 7}
    assert stats["n_points_3d"] == {"analysis": 0}
    assert stats["materialized_plotting_rows"] == 7
    assert landscape.data["particle_key"].tolist() == [f"p{index}" for index in range(5, 12)]


def test_one_d_materializes_every_filtered_row_even_with_small_point_cap(tmp_path) -> None:
    source = write_landscape_npz_arrays(_arrays(rows=120), tmp_path / "landscape.npz")
    request = VisualizationRequest(
        run_dir=str(tmp_path),
        representation="rotvec",
        view="1d",
        threshold=1.0,
        max_points=5,
    )

    landscape, _metadata, stats, _scale = _prepare_array_native_visualization(
        request,
        source,
        euler_metadata={"scipy_euler_sequence": "zyx"},
        range_bounds={},
        selected_landscape_report=None,
    )

    assert stats["n_points_after_display_filter"] == {"analysis": 119}
    assert stats["n_points_1d"] == {"analysis": 119}
    assert stats["materialized_plotting_rows"] == 119
    assert len(landscape.data) == 119
