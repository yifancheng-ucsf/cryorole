from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from cryorole.export.visualization import write_landscape_visualizations
from cryorole.models.landscape import Landscape


def _make_landscape(*, canonical: bool = False, n_points: int = 3) -> Landscape:
    base = np.arange(n_points, dtype=float)
    data = pd.DataFrame(
        {
            "particle_key": [f"p{i}" for i in range(n_points)],
            "coordinates_analysis": [
                np.array([value, value / 2.0, value / 3.0]) for value in base
            ],
            "coordinates_display": [
                np.array([value, value / 2.0, value / 3.0]) for value in base
            ],
            "sld_unfloored": base + 1.0,
            "sld_raw": base + 1.0,
            "sld_display": (base + 1.0) / 10.0,
            "sld_was_floored": [False] * n_points,
            "sld_local_k_mean": [1.0] * n_points,
            "sld_effective_local_k_mean": [1.0] * n_points,
            "sld_distance_floor": [0.0001] * n_points,
        }
    )
    if canonical:
        data["coordinates_canonical"] = [
            np.array([value / 3.0, value / 2.0, value]) for value in base
        ]
    return Landscape(data=data)


def test_visualization_writes_display_table_figures_and_report(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(),
        tmp_path,
        coordinate_source="analysis",
        representation="both",
        formats=("png",),
    )

    assert report["status"] == "ok"
    assert report["backend"] == "matplotlib_static"
    assert report["color_field"] == "sld_display"
    assert report["color_field_usage"] == "display_only"
    assert report["display_coordinates_contains"] == "display_filtered_rows_only"
    assert report["figure_point_rasterized"] is True
    assert report["euler_axis_limits"] == [-180.0, 180.0]
    assert report["euler_convention"] == "extrinsic_zyx"
    assert report["scipy_euler_sequence"] == "zyx"
    assert report["euler_convention_source"] == "cli_default"
    assert report["coordinate_source_resolved"] == ["analysis"]
    assert report["representations"] == ["euler", "rotvec"]
    assert (tmp_path / "display_coordinates.csv").exists()
    assert (tmp_path / "analysis_euler_projection_2d.png").exists()
    assert (tmp_path / "analysis_euler_cloud_3d.png").exists()
    assert (tmp_path / "analysis_rotvec_projection_2d.png").exists()
    assert (tmp_path / "analysis_rotvec_cloud_3d.png").exists()
    payload = json.loads((tmp_path / "visualization_report.json").read_text(encoding="utf-8"))
    assert payload["generated_files"]["display_coordinates_csv"] == str(
        tmp_path / "display_coordinates.csv"
    )


def test_visualization_records_explicit_intrinsic_euler_convention(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(),
        tmp_path,
        coordinate_source="analysis",
        representation="euler",
        euler_convention="intrinsic_zyx",
        euler_convention_source="cli_override",
        formats=("png",),
    )

    assert report["euler_convention"] == "intrinsic_zyx"
    assert report["scipy_euler_sequence"] == "ZYX"
    assert report["euler_angle_columns"] == [
        "raw_ea_zyx_alpha_deg",
        "raw_ea_zyx_beta_deg",
        "raw_ea_zyx_gamma_deg",
    ]


def test_visualization_report_includes_selection_metadata(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(n_points=2),
        tmp_path,
        coordinate_source="analysis",
        representation="rotvec",
        formats=("png",),
        selection_metadata={
            "selection_id": "subset",
            "selected_count": 2,
            "total_count": 3,
            "selection_filter_applied": True,
            "landscape_count_before_selection_filter": 3,
            "landscape_count_after_selection_filter": 2,
            "missing_selected_key_count": 0,
        },
    )

    assert report["selection_id"] == "subset"
    payload = json.loads((tmp_path / "visualization_report.json").read_text(encoding="utf-8"))
    assert payload["selected_count"] == 2
    assert payload["total_count"] == 3
    assert payload["landscape_count_after_selection_filter"] == 2


def test_visualization_debug_projection_csvs_are_optional(tmp_path) -> None:
    write_landscape_visualizations(
        _make_landscape(),
        tmp_path,
        coordinate_source="analysis",
        representation="rotvec",
        formats=("png",),
    )
    assert not (tmp_path / "analysis_projection_xy.csv").exists()

    debug_dir = tmp_path / "debug"
    write_landscape_visualizations(
        _make_landscape(),
        debug_dir,
        coordinate_source="analysis",
        representation="rotvec",
        formats=("png",),
        write_projection_csvs=True,
    )
    rows = (debug_dir / "analysis_projection_xy.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    assert rows[0] == "particle_key,x,y,sld_display"


def test_visualization_coordinate_source_canonical_requires_canonical_coordinates(tmp_path) -> None:
    with pytest.raises(ValueError, match="coordinates_canonical"):
        write_landscape_visualizations(
            _make_landscape(),
            tmp_path,
            coordinate_source="canonical",
            formats=("png",),
        )


def test_visualization_coordinate_source_canonical_if_available_prefers_canonical(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(canonical=True),
        tmp_path,
        coordinate_source="canonical_if_available",
        representation="euler",
        formats=("png",),
    )

    assert report["coordinate_source_resolved"] == ["canonical"]
    assert not (tmp_path / "analysis_euler_projection_2d.png").exists()
    assert (tmp_path / "canonical_euler_projection_2d.png").exists()


def test_visualization_does_not_mutate_landscape(tmp_path) -> None:
    landscape = _make_landscape()
    before = landscape.data.copy(deep=True)

    write_landscape_visualizations(landscape, tmp_path, formats=("png",))

    pd.testing.assert_frame_equal(landscape.data, before)


def test_visualization_display_filter_and_range_are_display_only(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(n_points=5),
        tmp_path,
        coordinate_source="analysis",
        representation="rotvec",
        display_top_fraction=0.40,
        range_bounds={"x": (3.0, None)},
        formats=("png",),
    )

    display = pd.read_csv(tmp_path / "display_coordinates.csv")
    assert len(display) == 2
    assert set(display["particle_key"]) == {"p3", "p4"}
    assert report["n_points_input"] == 5
    assert report["n_points_after_display_filter"]["analysis"] == 2
    assert report["display_filter_mode"] == "top_fraction"
    assert report["display_range_filter_applied"] is True


def test_visualization_3d_downsampling_is_deterministic(tmp_path) -> None:
    kwargs = dict(
        landscape=_make_landscape(n_points=20),
        coordinate_source="analysis",
        representation="rotvec",
        formats=("png",),
        max_points_3d=5,
        random_seed=7,
    )
    report_a = write_landscape_visualizations(output_dir=tmp_path / "a", **kwargs)
    report_b = write_landscape_visualizations(output_dir=tmp_path / "b", **kwargs)

    display_a = pd.read_csv(tmp_path / "a" / "display_coordinates.csv")
    display_b = pd.read_csv(tmp_path / "b" / "display_coordinates.csv")
    assert display_a["particle_key"].tolist() == display_b["particle_key"].tolist()
    assert report_a["n_points_3d"]["analysis"] == 5
    assert report_a["downsampled"]["analysis"] is True
    assert report_a["random_seed"] == 7


def test_visualization_2d_downsampling_is_reported_deterministically(tmp_path) -> None:
    kwargs = dict(
        landscape=_make_landscape(n_points=20),
        coordinate_source="analysis",
        representation="rotvec",
        formats=("png",),
        max_points_2d=4,
        random_seed=11,
    )
    report_a = write_landscape_visualizations(output_dir=tmp_path / "a", **kwargs)
    report_b = write_landscape_visualizations(output_dir=tmp_path / "b", **kwargs)

    assert report_a["n_points_2d"] == report_b["n_points_2d"]
    assert report_a["n_points_2d"]["analysis"] == 4
    assert report_a["random_seed"] == 11


def test_visualization_rejects_unsupported_formats(tmp_path) -> None:
    with pytest.raises(ValueError, match="Unsupported figure format"):
        write_landscape_visualizations(
            _make_landscape(),
            tmp_path,
            formats=("jpg",),
        )


def test_legacy_visualization_style_records_legacy_compatible_parameters(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(),
        tmp_path,
        coordinate_source="analysis",
        representation="both",
        visual_style="legacy",
        formats=("png",),
    )

    assert report["visual_style"] == "legacy"
    assert report["color_map"] == "rainbow_r"
    assert report["point_size"] == 1.0
    assert report["colorbar_position"] == "bottom"
    assert report["aspect"] == "equal"
    assert report["sort_points_by_color"] == "ascending"
    assert report["figure_size"]["projection_2d"] == [20.0, 9.0]
    assert report["axis_limits"]["alpha"] == [-180.0, 180.0]
    assert report["axis_limits"]["x"] == [-3.14, 3.14]


def test_legacy_visualization_reports_color_scale_and_max_divisor_filter(tmp_path) -> None:
    landscape = _make_landscape(n_points=5)

    report = write_landscape_visualizations(
        landscape,
        tmp_path,
        coordinate_source="analysis",
        representation="rotvec",
        visual_style="legacy",
        display_filter_mode="legacy_max_divisor",
        display_max_divisor=3.0,
        formats=("png",),
    )

    display = pd.read_csv(tmp_path / "display_coordinates.csv")
    assert report["color_vmin"] == 0.2
    assert report["color_vmin_source"] == "legacy_display_threshold"
    assert report["color_vmax"] == 0.5
    assert report["color_vmax_source"] == "legacy_color_max_ceiling"
    assert report["display_filter_mode"] == "legacy_max_divisor"
    assert report["display_filter_threshold"] == 0.2
    assert len(display) == 4
    assert len(landscape.data) == 5


def test_explicit_color_bounds_record_explicit_sources(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(n_points=5),
        tmp_path,
        coordinate_source="analysis",
        representation="rotvec",
        visual_style="legacy",
        display_filter_mode="legacy_max_divisor",
        color_vmin=0.1,
        color_vmax=1.0,
        formats=("png",),
    )

    assert report["color_vmin"] == 0.1
    assert report["color_vmax"] == 1.0
    assert report["color_vmin_source"] == "explicit"
    assert report["color_vmax_source"] == "explicit"


def test_visualization_auto_color_vmax_excludes_display_outliers(tmp_path) -> None:
    landscape = _make_landscape(n_points=5)
    landscape.data["sld_display"] = [1.0, 2.0, 3.0, 4.0, 100.0]
    landscape.data["sld_display_is_outlier"] = [False, False, False, False, True]

    report = write_landscape_visualizations(
        landscape,
        tmp_path,
        coordinate_source="analysis",
        representation="rotvec",
        formats=("png",),
    )

    assert report["color_vmax"] == 4.0
    assert report["color_vmax_source"] == "tail_jump_display_outliers"
    assert report["color_scale_mode"] == "tail_jump_exclude_outliers"
    assert report["display_color_vmax"] == 4.0
    assert report["n_sld_display_outliers"] == 1
    assert report["fraction_sld_display_outliers"] == 0.2


def test_axis_limits_do_not_filter_display_rows(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(n_points=5),
        tmp_path,
        coordinate_source="analysis",
        representation="rotvec",
        axis_limits={"x": (-1.0, 1.0)},
        formats=("png",),
    )

    display = pd.read_csv(tmp_path / "display_coordinates.csv")
    assert len(display) == 5
    assert report["axis_limits"]["x"] == [-1.0, 1.0]
    assert report["display_filter_mode"] == "none"
    assert report["display_range_filter_applied"] is False


def test_unused_axis_limits_are_reported_for_unselected_representation(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(n_points=5),
        tmp_path,
        coordinate_source="analysis",
        representation="rotvec",
        axis_limits={"alpha": (-90.0, 90.0), "x": (-1.0, 1.0)},
        formats=("png",),
    )

    assert report["axis_limits"]["alpha"] == [-90.0, 90.0]
    assert report["axis_limits"]["x"] == [-1.0, 1.0]
    assert report["unused_axis_limits"]["alpha"] == [-90.0, 90.0]
    assert "x" not in report["unused_axis_limits"]


def test_sort_points_by_color_is_reported_and_deterministic(tmp_path) -> None:
    kwargs = dict(
        landscape=_make_landscape(n_points=6),
        coordinate_source="analysis",
        representation="rotvec",
        sort_points_by_color="ascending",
        formats=("png",),
    )

    report_a = write_landscape_visualizations(output_dir=tmp_path / "a", **kwargs)
    report_b = write_landscape_visualizations(output_dir=tmp_path / "b", **kwargs)

    assert report_a["sort_points_by_color"] == "ascending"
    assert report_b["sort_points_by_color"] == "ascending"
    display_a = pd.read_csv(tmp_path / "a" / "display_coordinates.csv")
    display_b = pd.read_csv(tmp_path / "b" / "display_coordinates.csv")
    assert display_a["particle_key"].tolist() == display_b["particle_key"].tolist()


def test_optional_histograms_and_axis_direction_map_are_generated(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(n_points=5),
        tmp_path,
        coordinate_source="analysis",
        representation="rotvec",
        formats=("png",),
        generate_histograms=True,
        generate_axis_direction_map=True,
    )

    assert (tmp_path / "analysis_sld_histogram.png").exists()
    assert (tmp_path / "analysis_rotation_angle_histogram.png").exists()
    assert (tmp_path / "analysis_euler_histograms.png").exists()
    assert (tmp_path / "analysis_axis_direction_azimuth_elevation.png").exists()
    assert "analysis_sld_histogram_png" in report["generated_histograms"]
    assert "analysis_axis_direction_azimuth_elevation_png" in report[
        "generated_axis_direction_map"
    ]


def test_visualization_report_does_not_reintroduce_rnd_public_keys(tmp_path) -> None:
    report = write_landscape_visualizations(
        _make_landscape(),
        tmp_path,
        coordinate_source="analysis",
        representation="rotvec",
        visual_style="legacy",
        display_filter_mode="legacy_max_divisor",
        generate_histograms=True,
        formats=("png",),
    )

    assert "rnd" not in json.dumps(report).lower()
