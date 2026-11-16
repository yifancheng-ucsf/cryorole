"""M5 display artifacts: identity, read-only overlays, renderer parity and publication."""

from __future__ import annotations

import hashlib
import json
import shlex
import sys
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from cryorole import api
from cryorole.errors import CryoroleError

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run"))
from test_core_productionization import _write_cs  # noqa: E402


@pytest.fixture
def bundle(tmp_path):
    rng = np.random.default_rng(41)
    uids = list(range(240))
    _write_cs(tmp_path / "ref.cs", uids, rng.normal(0, 0.2, (240, 3)))
    _write_cs(tmp_path / "mov.cs", uids, rng.normal(0, 0.3, (240, 3)))
    run = tmp_path / "run"
    api.run(tmp_path / "ref.cs", tmp_path / "mov.cs", output_dir=run, k_neighbors=10, no_visualize=True)
    api.canonicalize(run, no_visualize=True)
    api.select(run, selection_id="sample", mode="random", fraction=0.25, seed=3)
    return run


def _hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob("*") if p.is_file()}


def test_overlay_maps_saved_keys_in_both_spaces_without_writing_or_creating_a_draft(bundle):
    before = _hashes(bundle)
    raw = api.open_landscape(bundle)
    canonical = api.open_landscape(bundle, space="canonical", canonical_id="default")
    a, b = (view.saved_selection_overlay("sample") for view in (raw, canonical))
    assert a.selected_count == 60 and a.total_count == 240
    assert np.array_equal(a.selected_mask, b.selected_mask)
    assert not a.selected_mask.flags.writeable
    assert raw.session.draft_selected_mask() is None
    assert _hashes(bundle) == before


@pytest.mark.parametrize("change", ["run", "unknown_key", "count"])
def test_overlay_rejects_another_run_or_damaged_keys(bundle, change):
    path = bundle / "selections/sample/selection.json"
    payload = json.loads(path.read_text())
    if change == "run":
        payload["parent_run_id"] = "another-run"
    elif change == "unknown_key":
        payload["selected_particle_keys"][0] = "not-in-this-run"
    else:
        payload["total_count"] += 1
    path.write_text(json.dumps(payload))
    with pytest.raises((ValueError, CryoroleError)):
        api.open_landscape(bundle).saved_selection_overlay("sample")


def test_overlay_rejects_path_traversal_and_changed_parent(bundle):
    view = api.open_landscape(bundle)
    with pytest.raises(CryoroleError):
        view.saved_selection_overlay("../sample")
    with (bundle / "data/raw_landscape.npz").open("ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError, match="landscape changed"):
        view.saved_selection_overlay("sample")


@pytest.mark.parametrize("space", ["raw", "canonical"])
def test_saved_figure_snapshot_matches_public_renderer_and_preserves_scientific_artifacts(bundle, space):
    view = api.open_landscape(bundle, space=space, canonical_id="default" if space == "canonical" else None)
    frame = view.apply(top_fraction=0.4, max_points_3d=5, color_vmax=4.0)
    before = _hashes(bundle)
    limits = ((-0.8, 0.8), (-0.6, 0.6), (-1.0, 1.0))
    result = view.save_figures("snapshot", frame=frame, representation="rotvec", axis_ranges=limits)
    report = json.loads((result.output_dir / "visualization_report.json").read_text())
    assert report["view_reports"]["1d"]["input_row_count"] == frame.filtered_count
    assert report["gui_display_snapshot"]["run_id"] == view.run_id
    assert report["gui_display_snapshot"]["highlight_overlays_included"] is False
    assert ".stage-" not in json.dumps(report)
    command = shlex.split(report["gui_display_snapshot"]["equivalent_cli_command"])
    assert command[command.index("--space") + 1] == space
    assert command[command.index("--representation") + 1] == "rotvec"
    standard = api.visualize(
        bundle, space=space, canonical_id=view.canonical_id, visual_id="parity", view="2d,1d",
        representation="rotvec", formats=("png",), top_fraction=0.4, vmin=frame.vmin, vmax=frame.vmax,
        max_points=view.n_points, axis_limit=tuple(zip(("x", "y", "z"), limits)),
    )
    pngs = list(result.output_dir.glob("*.png"))
    assert pngs
    for png in pngs:
        with Image.open(png) as image:
            assert min(image.size) > 100
        assert png.read_bytes() == (standard.output_dir / png.name).read_bytes()
    after = _hashes(bundle)
    assert all(after[name] == digest for name, digest in before.items())
    assert not any(p.name.startswith(".") for p in result.output_dir.parent.iterdir())
    with pytest.raises(FileExistsError):
        view.save_figures("snapshot", frame=frame, representation="rotvec", axis_ranges=limits)


def test_render_failure_cleans_only_its_stage_and_preserves_existing_outputs(bundle, monkeypatch):
    from cryorole.interactive import landscape_view as module

    view = api.open_landscape(bundle)
    frame = view.apply()
    before = _hashes(bundle)

    def fail(_request, *, _output_dir):
        (_output_dir / "partial.png").write_bytes(b"incomplete")
        raise RuntimeError("renderer stopped")

    monkeypatch.setattr(module, "visualize", fail)
    with pytest.raises(RuntimeError, match="renderer stopped"):
        view.save_figures("failed", frame=frame, representation="euler", axis_ranges=((-180, 180),) * 3)
    assert _hashes(bundle) == before
    parent = bundle / "visualizations/raw"
    assert not parent.exists() or not list(parent.iterdir())


def test_figures_refuse_a_parent_changed_during_rendering(bundle, monkeypatch):
    from cryorole.interactive import landscape_view as module

    view = api.open_landscape(bundle)
    frame = view.apply()
    actual = module.visualize

    def changed(request, *, _output_dir):
        result = actual(request, _output_dir=_output_dir)
        manifest = bundle / "run_manifest.json"
        payload = json.loads(manifest.read_text())
        payload["run_id"] = "replaced"
        manifest.write_text(json.dumps(payload))
        return result

    monkeypatch.setattr(module, "visualize", changed)
    with pytest.raises(ValueError, match="run_id changed"):
        view.save_figures("failed", frame=frame, representation="euler", axis_ranges=((-180, 180),) * 3)
    assert not (bundle / "visualizations/raw/failed").exists()
