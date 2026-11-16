"""Display model behind the GUI View page: same display policy as `cryorole visualize`."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole import api
from cryorole.core.display_policy import visible_value_raster
from cryorole.errors import CryoroleError

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run"))
from test_core_productionization import _write_cs  # noqa: E402

N = 300


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    root = tmp_path_factory.mktemp("landscape_view")
    rng = np.random.default_rng(11)
    uids = list(range(7000, 7000 + N))
    ref = rng.normal(0.0, 0.6, (N, 3))
    mov = (Rotation.from_rotvec(ref) * Rotation.from_rotvec(rng.normal(0.0, 0.15, (N, 3)))).as_rotvec()
    _write_cs(root / "ref.cs", uids, ref)
    _write_cs(root / "mov.cs", uids, mov)
    api.run(root / "ref.cs", root / "mov.cs", output_dir=root / "run", k_neighbors=10, no_visualize=True)
    api.canonicalize(root / "run", no_visualize=True)
    return root / "run"


def test_open_raw_landscape_exposes_identity_and_frames(run_dir):
    view = api.open_landscape(run_dir)
    assert view.space == "raw" and view.canonical_id is None
    assert view.n_points == N
    assert view.canonical_ids == ["default"]
    assert view.run_id == view.session.run_id
    assert view.euler_convention == "extrinsic_zyx"


def test_default_frame_shows_everything_and_samples_3d_deterministically(run_dir):
    view = api.open_landscape(run_dir)
    frame = view.apply(max_points_3d=100)
    assert frame.filtered_count == N and frame.full_count == N
    assert frame.sample_count == 100 and frame.downsampled
    sld = view.sld_display[frame.sample_indices]
    assert np.all(np.diff(sld) >= 0)  # draw order: ascending SLD, as in the figures
    again = view.apply(max_points_3d=100)
    assert np.array_equal(np.sort(again.sample_indices), np.sort(frame.sample_indices))
    assert frame.vmax == pytest.approx(min(float(view.sld_display.max()), 100.0))
    assert frame.vmin == pytest.approx(float(view.sld_display.min()))


def test_display_policy_matches_cryorole_visualize(run_dir):
    report = api.visualize(run_dir, threshold=1.0, view="2d,1d", max_points=100, visual_id="parity").report
    frame = api.open_landscape(run_dir).apply(threshold=1.0, max_points_3d=100)
    assert frame.filtered_count == report["view_reports"]["1d"]["input_row_count"]
    assert frame.vmin == pytest.approx(report["color_vmin"])
    assert frame.vmax == pytest.approx(report["color_vmax"])
    view = api.open_landscape(run_dir)
    assert view.style["color_map"] == report["color_map"]
    assert view.style["sort_points_by_color"] == report["sort_points_by_color"]


def test_top_fraction_keeps_ties_and_color_vmax_is_display_only(run_dir):
    view = api.open_landscape(run_dir)
    frame = view.apply(top_fraction=0.25)
    cutoff = np.sort(view.sld_display)[-int(np.ceil(N * 0.25))]
    assert frame.filtered_count == int(np.count_nonzero(view.sld_display >= cutoff))
    capped = view.apply(top_fraction=0.25, color_vmax=1.5)
    assert capped.vmax == 1.5
    assert np.array_equal(capped.filtered_indices, frame.filtered_indices)


def test_unusable_display_settings_are_cryorole_errors(run_dir):
    view = api.open_landscape(run_dir)
    with pytest.raises(CryoroleError) as empty:
        view.apply(threshold=1e9)
    assert empty.value.code == "display_empty" and empty.value.remedy
    with pytest.raises(CryoroleError) as bad_vmax:
        view.apply(threshold=1.0, color_vmax=0.5)
    assert bad_vmax.value.code == "input_invalid"
    with pytest.raises(CryoroleError):
        view.apply(max_points_3d=0)


def test_coordinates_use_the_recorded_euler_convention(run_dir):
    view = api.open_landscape(run_dir)
    rotvec = view.coordinates("rotvec")
    assert np.array_equal(rotvec, view.session.arrays.coordinates_analysis)
    expected = Rotation.from_rotvec(rotvec).as_euler(view.session.euler_sequence, degrees=True)
    assert np.allclose(view.coordinates("euler"), expected)
    assert view.axis_names("euler") == ("alpha", "beta", "gamma")
    assert view.axis_limits("rotvec") == ((-3.14, 3.14),) * 3


def test_projection_images_and_histograms_use_every_filtered_row(run_dir):
    view = api.open_landscape(run_dir)
    frame = view.apply(threshold=1.0, max_points_3d=10)
    coords = view.coordinates("euler", frame.filtered_indices)
    image = view.projection_image(
        frame, "euler", (0, 1), x_range=(-180, 180), y_range=(-180, 180), width=90, height=60
    )
    expected = visible_value_raster(
        coords[:, 0], coords[:, 1], view.sld_display[frame.filtered_indices],
        x_range=(-180, 180), y_range=(-180, 180), width=90, height=60,
    )
    assert np.array_equal(image, expected, equal_nan=True)
    counts, edges = view.histogram(frame, "rotvec", 2)
    assert counts.sum() == pytest.approx(100.0) and len(edges) == len(counts) + 1
    sld_counts, _ = view.sld_histogram(frame)
    assert sld_counts.sum() == pytest.approx(100.0)


def test_canonical_landscape_uses_canonical_coordinates(run_dir):
    view = api.open_landscape(run_dir, space="canonical")  # only one frame: resolved automatically
    assert view.space == "canonical" and view.canonical_id == "default"
    assert np.array_equal(view.coordinates("rotvec"), view.session.arrays.coordinates_canonical)


def test_wrong_folders_get_a_clear_remedy(tmp_path):
    with pytest.raises(CryoroleError) as missing:
        api.open_landscape(tmp_path / "nothing_here")
    assert missing.value.code == "input_not_found" and "cryorole run" in missing.value.remedy
    with pytest.raises(CryoroleError) as not_a_bundle:
        api.open_landscape(tmp_path)
    assert "is not a cryoROLE run bundle" in not_a_bundle.value.message
    assert "run_manifest.json" in not_a_bundle.value.remedy


# --- radius selection (M2) ---------------------------------------------------------

from cryorole.export import read_selection_json  # noqa: E402


def _keys(run_dir, selection_id):
    selection = read_selection_json(Path(run_dir) / "selections" / selection_id / "selection.json")
    return selection, set(map(str, selection.selected_particle_keys))


@pytest.mark.parametrize("space", ["raw", "canonical"])
def test_gui_radius_selection_matches_cryorole_select(run_dir, space):
    view = api.open_landscape(run_dir, space=space)
    center_rv = tuple(view.coordinates("rotvec", [7])[0])
    draft = view.evaluate_radius(center_rv, representation="rotvec", radius_deg=12.0)
    assert draft.total_count == N and 0 < draft.selected_count < N
    assert draft.selected_mask.shape == (N,) and int(draft.selected_mask.sum()) == draft.selected_count
    view.confirm(f"gui_{space}", draft)
    api.select(run_dir, selection_id=f"cli_{space}", mode="radius", space=space,
               center=center_rv, center_representation="rotvec", radius=12.0)
    gui_selection, gui_keys = _keys(run_dir, f"gui_{space}")
    _, cli_keys = _keys(run_dir, f"cli_{space}")
    assert gui_keys == cli_keys and len(gui_keys) == draft.selected_count
    assert gui_selection.interaction_provenance["source"] == "cryorole_gui"
    assert gui_selection.metric == "so3_geodesic" and gui_selection.radius == 12.0
    assert f"--selection-id gui_{space}" in draft.command(f"gui_{space}")


def test_euler_centre_gives_the_same_selection_as_its_rotation_vector(run_dir):
    view = api.open_landscape(run_dir)
    center_rv = view.coordinates("rotvec", [3])[0]
    euler = view.represent(center_rv, "euler")
    assert np.allclose(euler, view.coordinates("euler", [3])[0])
    from_rv = view.evaluate_radius(center_rv, representation="rotvec", radius_deg=8.0)
    from_euler = view.evaluate_radius(euler, representation="euler", radius_deg=8.0)
    assert np.array_equal(from_rv.selected_mask, from_euler.selected_mask)


def test_confirm_refuses_bad_names_existing_names_and_stale_counts(run_dir):
    view = api.open_landscape(run_dir)
    first = view.evaluate_radius((0.0, 0.0, 0.0), representation="rotvec", radius_deg=20.0)
    with pytest.raises(CryoroleError) as bad:
        view.confirm("../escape", first)
    assert bad.value.code == "input_invalid" and "region_01" in bad.value.remedy
    newer = view.evaluate_radius((0.0, 0.0, 0.0), representation="rotvec", radius_deg=25.0)
    with pytest.raises(CryoroleError) as stale:
        view.confirm("stale_one", first)
    assert "changed after this count" in stale.value.message
    assert not (Path(run_dir) / "selections" / "stale_one").exists()
    view.confirm("kept_once", newer)
    again = view.evaluate_radius((0.0, 0.0, 0.0), representation="rotvec", radius_deg=5.0)
    with pytest.raises(CryoroleError) as exists:
        view.confirm("kept_once", again)
    assert exists.value.code == "output_exists" and "never overwrites" in exists.value.remedy
    _, keys = _keys(run_dir, "kept_once")
    assert len(keys) == newer.selected_count  # the first confirm is untouched
    assert "kept_once" in view.selection_ids()


def test_subset_images_and_histograms_share_bins_and_scale(run_dir):
    view = api.open_landscape(run_dir)
    frame = view.apply(threshold=None)
    draft = view.evaluate_radius((0.0, 0.0, 0.0), representation="rotvec", radius_deg=30.0)
    rows = frame.filtered_indices[draft.selected_mask[frame.filtered_indices]]
    full, edges = view.histogram(frame, "rotvec", 0)
    part, part_edges = view.histogram(frame, "rotvec", 0, rows=rows, edges=edges)
    assert np.array_equal(edges, part_edges)
    assert part.sum() == pytest.approx(100.0 * rows.size / frame.filtered_count)
    assert np.all(part <= full + 1e-9)
    image = view.projection_image(frame, "rotvec", (0, 1), x_range=(-3, 3), y_range=(-3, 3),
                                  width=40, height=40, rows=rows)
    assert np.isfinite(image).sum() > 0


def test_saved_selections_list_names_and_counts_from_summaries(run_dir):
    view = api.open_landscape(run_dir)
    draft = view.evaluate_radius((0.0, 0.0, 0.0), representation="rotvec", radius_deg=15.0)
    view.confirm("listed_gui", draft)
    api.select(run_dir, selection_id="listed_cli", mode="threshold", sld_min=1.0)
    saved = dict(view.saved_selections())
    assert saved["listed_gui"] == draft.selected_count
    assert saved["listed_cli"] == int(np.count_nonzero(view.session.arrays.sld_raw >= 1.0))
