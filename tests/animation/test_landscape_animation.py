from __future__ import annotations

import json

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy.spatial.transform import Rotation

matplotlib.use("Agg")

from cryorole.animation.artifacts import (
    load_animation_landscape,
    load_canonical_frame,
    resolve_recorded_euler_convention,
)
from cryorole.animation.projection import ProjectionRenderer, resolve_axis_limits
from cryorole.animation.schemas import TrajectoryFrame
from cryorole.core.display_policy import EULER_PROJECTIONS
from cryorole.export.visualization import _display_indices, _resolve_color_scale


def _write_landscape(run_dir, *, canonical=False):
    target = (
        run_dir / "canonical" / "default" / "canonical_landscape.npz"
        if canonical
        else run_dir / "data" / "raw_landscape.npz"
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    raw = np.asarray([[0.0, 0.0, 0.0], [0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
    payload = {
        "artifact_type": np.asarray("canonical_landscape" if canonical else "raw_landscape"),
        "schema_version": np.asarray("1"),
        "particle_key": np.asarray(["a", "b", "c"]),
        "coordinates_analysis": raw,
        "sld_display": np.asarray([1.0, 2.0, 10.0]),
        "sld_display_is_outlier": np.asarray([False, False, True]),
    }
    if canonical:
        payload["coordinates_canonical"] = raw @ Rotation.from_euler(
            "z", 20, degrees=True
        ).as_matrix()
        payload["canonical_transform"] = Rotation.from_euler(
            "z", 20, degrees=True
        ).as_matrix()
    np.savez(target, **payload)
    return target


def _write_provenance(run_dir, convention="extrinsic_zyx"):
    (run_dir / "run_summary.json").write_text(
        json.dumps({"euler_convention": convention}),
        encoding="utf-8",
    )


def test_raw_and_canonical_artifact_resolution_and_filtering(tmp_path):
    run_dir = tmp_path / "run"
    _write_landscape(run_dir)
    _write_landscape(run_dir, canonical=True)
    _write_provenance(run_dir)

    raw = load_animation_landscape(
        run_dir,
        coordinate_set="raw",
        canonical_id=None,
        scipy_euler_sequence="zyx",
        sld_threshold=1.5,
        top_fraction=None,
    )
    canonical = load_animation_landscape(
        run_dir,
        coordinate_set="canonical",
        canonical_id="default",
        scipy_euler_sequence="zyx",
        sld_threshold=None,
        top_fraction=1 / 3,
    )
    assert raw.total_rows == canonical.total_rows == 3
    assert raw.displayed_rows == 2
    assert canonical.displayed_rows == 1
    assert raw.sld_field == "sld_display"
    assert raw.color_scale.vmin == 1.5
    assert raw.color_scale.vmax == 2.0


def test_animation_and_visualize_threshold_and_top_fraction_row_parity(tmp_path):
    run_dir = tmp_path / "run"
    source = _write_landscape(run_dir)
    with np.load(source, allow_pickle=False) as payload:
        values = {key: payload[key] for key in payload.files}
    values["sld_display"] = np.asarray([2.0, 2.0, 3.0])
    np.savez(source, **values)

    threshold = load_animation_landscape(
        run_dir,
        coordinate_set="raw",
        canonical_id=None,
        scipy_euler_sequence="zyx",
        sld_threshold=2.0,
        top_fraction=None,
    )
    top = load_animation_landscape(
        run_dir,
        coordinate_set="raw",
        canonical_id=None,
        scipy_euler_sequence="zyx",
        sld_threshold=None,
        top_fraction=2 / 3,
    )
    all_euler = Rotation.from_rotvec(values["coordinates_analysis"]).as_euler(
        "zyx", degrees=True
    )
    table = __import__("pandas").DataFrame({"sld_display": values["sld_display"]})
    visualize_threshold = _display_indices(
        table,
        values["coordinates_analysis"],
        all_euler,
        representations=["euler"],
        display_top_fraction=None,
        display_sld_threshold=2.0,
        display_density_field="sld_display",
        range_bounds={},
    )
    visualize_top = _display_indices(
        table,
        values["coordinates_analysis"],
        all_euler,
        representations=["euler"],
        display_top_fraction=2 / 3,
        display_sld_threshold=None,
        display_density_field="sld_display",
        range_bounds={},
    )
    assert threshold.source_row_indices.tolist() == visualize_threshold.tolist()
    assert top.source_row_indices.tolist() == visualize_top.tolist()


def test_strict_euler_inheritance_and_conflict(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    _write_provenance(run_dir, "intrinsic_zyx")
    resolved = resolve_recorded_euler_convention(run_dir, coordinate_set="raw")
    assert resolved.euler_convention == "intrinsic_zyx"
    with pytest.raises(ValueError, match="conflicts"):
        resolve_recorded_euler_convention(
            run_dir,
            coordinate_set="raw",
            requested="extrinsic_zyx",
        )
    (run_dir / "run_summary.json").unlink()
    with pytest.raises(ValueError, match="provenance"):
        resolve_recorded_euler_convention(run_dir, coordinate_set="raw")


def test_canonical_frame_contract_validation(tmp_path):
    run_dir = tmp_path / "run"
    frame_dir = run_dir / "canonical" / "c1"
    frame_dir.mkdir(parents=True)
    transform = Rotation.from_euler("x", 30, degrees=True).as_matrix()
    (frame_dir / "canonical_frame.json").write_text(
        json.dumps(
            {
                "artifact_type": "canonical_frame",
                "canonical_transform": transform.tolist(),
                "transform_direction": "canonical_rv = raw_rv @ canonical_transform",
                "coordinate_space": "rotvec_ro_radians",
                "source_canonical_id": "c1",
                "physical_change_of_basis": True,
            }
        ),
        encoding="utf-8",
    )
    loaded = load_canonical_frame(run_dir, "c1")
    assert np.allclose(loaded.transform, transform)
    assert loaded.raw_to_canonical_column_convention == "x_canonical = C.T @ x_raw"
    assert loaded.physical_change_of_basis is True


def test_canonical_frame_is_not_implicitly_a_physical_map_frame(tmp_path):
    run_dir = tmp_path / "run"
    frame_dir = run_dir / "canonical" / "c1"
    frame_dir.mkdir(parents=True)
    (frame_dir / "canonical_frame.json").write_text(
        json.dumps(
            {
                "canonical_transform": np.eye(3).tolist(),
                "transform_direction": "canonical_rv = raw_rv @ canonical_transform",
                "coordinate_space": "rotvec_ro_radians",
            }
        ),
        encoding="utf-8",
    )
    loaded = load_canonical_frame(run_dir, "c1")
    assert loaded.physical_change_of_basis is False


def _frames(annotation=""):
    rotations = Rotation.from_euler("zyx", [[0, 0, 0], [15, 25, 35]], degrees=True)
    return tuple(
        TrajectoryFrame(
            frame_index=index,
            time_sec=index / 10,
            segment_index=0,
            from_label="a",
            to_label="b",
            segment_fraction=float(index),
            rotation=rotation,
            is_waypoint=True,
            waypoint_label=("a", "b")[index],
            annotation=annotation if index == 1 else "",
        )
        for index, rotation in enumerate(rotations)
    )


def test_axis_expansion_explicit_failure_and_cached_rendering(tmp_path):
    landscape_euler = np.asarray([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
    trajectory_euler = np.asarray([[0.0, 0.0, 0.0], [15.0, 25.0, 35.0]])
    limits, expansions = resolve_axis_limits(
        landscape_euler,
        trajectory_euler,
        explicit=None,
        range_bounds={},
    )
    assert expansions == ()
    assert limits["alpha_beta"] == (-180.0, 180.0, -180.0, 180.0)
    with pytest.raises(ValueError, match="outside explicit"):
        resolve_axis_limits(
            landscape_euler,
            trajectory_euler,
            explicit={
                "alpha_beta": [-1, 1, -1, 1],
                "alpha_gamma": [-1, 1, -1, 1],
                "beta_gamma": [-1, 1, -1, 1],
            },
            range_bounds={},
        )
    with pytest.raises(ValueError, match="outside explicit display range"):
        resolve_axis_limits(
            landscape_euler,
            trajectory_euler,
            explicit=None,
            range_bounds={"alpha": (-5.0, 5.0)},
        )

    renderer = ProjectionRenderer(
        landscape_euler=landscape_euler,
        sld_display=np.asarray([1.0, 2.0]),
        trajectory_frames=_frames(),
        scipy_euler_sequence="zyx",
        axis_limits=limits,
        canvas_width=480,
        canvas_height=240,
        point_size=2,
        landscape_alpha=0.5,
        trail_frames=1,
        show_current_values=True,
    )
    outputs = renderer.render(tmp_path / "frames")
    assert renderer.static_scatter_count == 3
    assert renderer.static_background_draw_count == 1
    assert len(outputs) == len(_frames())
    assert [path.name for path in outputs] == ["frame_000000.png", "frame_000001.png"]
    assert len({path.read_bytes()[16:24] for path in outputs}) == 1


def test_projection_style_color_scale_and_static_selection_are_shared(tmp_path):
    landscape_euler = np.asarray([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [2.0, 3.0, 4.0]])
    sld = np.asarray([1.0, 2.0, 100.0])
    outliers = np.asarray([False, False, True])
    limits, _ = resolve_axis_limits(
        landscape_euler,
        np.asarray([[0.0, 0.0, 0.0], [15.0, 25.0, 35.0]]),
        explicit=None,
        range_bounds={},
    )
    renderer = ProjectionRenderer(
        landscape_euler=landscape_euler,
        sld_display=sld,
        sld_display_is_outlier=outliers,
        trajectory_frames=_frames(),
        scipy_euler_sequence="zyx",
        axis_limits=limits,
        canvas_width=480,
        canvas_height=240,
        point_size=1,
        landscape_alpha=1,
        colormap="rainbow_r",
        color_vmin=None,
        color_vmax=None,
        trail_frames=1,
        show_current_values=True,
    )
    visualize_scale = _resolve_color_scale(
        sld,
        display_outlier_mask=outliers,
        visual_style="legacy",
        color_vmin=None,
        color_vmax=None,
        display_threshold=None,
    )
    before = renderer.static_source_row_count
    renderer.render(tmp_path / "frames")
    assert tuple(item[0] for item in EULER_PROJECTIONS) == tuple(
        item[0] for item in renderer.panel_specs
    )
    assert renderer.resolved_vmin == visualize_scale[0]
    assert renderer.resolved_vmax == visualize_scale[1]
    assert renderer.static_source_row_count == before == 3
    assert renderer.static_background_draw_count == 1
    assert all(axis.get_aspect() == 1.0 for axis in renderer.axes)
    assert renderer.colorbar.orientation == "horizontal"
    assert renderer.colorbar.ax.get_xlabel() == "SLD"
    assert renderer.colorbar_label == "SLD"
    assert renderer.sld_field == "sld_display"


def test_projection_geometry_is_compact_and_uses_one_coordinate_only_annotation(tmp_path):
    canvas_width = 900
    landscape_euler = np.asarray([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
    frames = _frames(annotation="waypoint narrative")
    limits, _ = resolve_axis_limits(
        landscape_euler,
        np.asarray([[0.0, 0.0, 0.0], [15.0, 25.0, 35.0]]),
        explicit=None,
        range_bounds={},
    )
    renderer = ProjectionRenderer(
        landscape_euler=landscape_euler,
        sld_display=np.asarray([1.0, 2.0]),
        trajectory_frames=frames,
        scipy_euler_sequence="zyx",
        axis_limits=limits,
        canvas_width=canvas_width,
        canvas_height=300,
        point_size=1,
        landscape_alpha=1,
        trail_frames=1,
        show_current_values=True,
    )
    renderer.render(tmp_path / "shown")

    widths = [rectangle[2] for rectangle in renderer.panel_rectangles]
    gaps = [
        right[0] - (left[0] + left[2])
        for left, right in zip(
            renderer.panel_rectangles,
            renderer.panel_rectangles[1:],
        )
    ]
    assert max(widths) - min(widths) <= 1
    assert all(0 <= gap <= canvas_width * 0.0125 for gap in gaps)
    assert len(renderer.annotations) == 1
    assert renderer.annotations[0].get_text() == "α=15.00 β=25.00 γ=35.00"
    assert "waypoint narrative" not in renderer.annotations[0].get_text()
    assert len(renderer.figure.axes) == 4
    assert renderer.colorbar.orientation == "horizontal"
    assert renderer.shared_colorbar_policy == "single_horizontal_centered"


def test_projection_marker_and_font_scale_are_legible_and_point_size_independent():
    landscape_euler = np.asarray([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
    limits, _ = resolve_axis_limits(
        landscape_euler,
        np.asarray([[0.0, 0.0, 0.0], [15.0, 25.0, 35.0]]),
        explicit=None,
        range_bounds={},
    )

    def build(point_size):
        return ProjectionRenderer(
            landscape_euler=landscape_euler,
            sld_display=np.asarray([1.0, 2.0]),
            trajectory_frames=_frames(),
            scipy_euler_sequence="zyx",
            axis_limits=limits,
            canvas_width=1800,
            canvas_height=600,
            point_size=point_size,
            landscape_alpha=1,
            trail_frames=1,
            show_current_values=True,
            composite_scale=0.58,
        )

    small = build(1)
    large = build(12)
    assert small.markers[0].get_sizes()[0] == pytest.approx(
        large.markers[0].get_sizes()[0]
    )
    assert small.resolved_marker_diameter_pixels == pytest.approx(13.0)
    assert small.resolved_current_value_font_size_pixels == pytest.approx(17.0)
    assert (
        np.sqrt(small.marker_size_points2) * 100 / 72 * 0.58
    ) == pytest.approx(13.0)
    assert (
        small.current_value_font_size_points * 100 / 72 * 0.58
    ) == pytest.approx(17.0)
    assert small.annotation.figure is small.figure
    assert small.annotation.get_bbox_patch().get_alpha() < 1
    plt.close(small.figure)
    plt.close(large.figure)


def test_projection_annotation_is_empty_when_current_values_are_disabled(tmp_path):
    landscape_euler = np.asarray([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]])
    frames = _frames(annotation="waypoint narrative")
    limits, _ = resolve_axis_limits(
        landscape_euler,
        np.asarray([[0.0, 0.0, 0.0], [15.0, 25.0, 35.0]]),
        explicit=None,
        range_bounds={},
    )
    renderer = ProjectionRenderer(
        landscape_euler=landscape_euler,
        sld_display=np.asarray([1.0, 2.0]),
        trajectory_frames=frames,
        scipy_euler_sequence="zyx",
        axis_limits=limits,
        canvas_width=900,
        canvas_height=300,
        point_size=1,
        landscape_alpha=1,
        trail_frames=0,
        show_current_values=False,
    )
    renderer.render(tmp_path / "hidden")
    assert len(renderer.annotations) == 1
    assert renderer.annotations[0].get_text() == ""
