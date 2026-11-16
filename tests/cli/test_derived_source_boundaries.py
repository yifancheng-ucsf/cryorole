"""Derived display/frame outputs preserve the input artifacts they consume."""
import json
import os

import numpy as np
import pytest

from cryorole import api
from cryorole.cli.main import main
from cryorole.workflows.rotate_landscape import main as rotate_main


@pytest.fixture
def bundle(tmp_path):
    data = np.zeros(30, dtype=[("uid", "u8"), ("alignments3D/pose", "f8", (3,))])
    data["uid"] = np.arange(100, 130)
    for seed, domain in enumerate(["ref", "mov"]):
        data["alignments3D/pose"] = np.random.default_rng(seed).normal(0, .2, (30, 3))
        with (tmp_path / f"{domain}.cs").open("wb") as stream:
            np.save(stream, data)
    run = tmp_path / "run"
    api.run(tmp_path / "ref.cs", tmp_path / "mov.cs", output_dir=run,
            no_visualize=True, no_raw_csv=True, k_neighbors=3)
    api.canonicalize(run, no_visualize=True, fit_top_fraction=1)
    (tmp_path / "scene.cxs").write_bytes(b"script-only fixture")
    (tmp_path / "path.csv").write_text("label,rv_x_rad,rv_y_rad,rv_z_rad\nstart,0,0,0\nend,.1,.2,.3\n")
    config = tmp_path / "config"
    config.mkdir()
    (config / "basis.json").write_text(json.dumps(np.eye(3).tolist()))
    (config / "limits.json").write_text(json.dumps({"alpha": [-180, 180], "beta": [-90, 90], "gamma": [-180, 180]}))
    return tmp_path, run, config


def _animate(root, run):
    return ["animate", "--run-dir", str(run), "--coordinate-set", "raw", "--path-csv", str(root / "path.csv"),
            "--path-space", "rv", "--chimerax-session", str(root / "scene.cxs"), "--reference-model-id", "#1",
            "--moving-model-id", "#2", "--pivot", "0", "0", "0", "--baseline-ro", "identity",
            "--map-frame", "raw", "--frames-per-segment", "2", "--canvas-width", "600", "--canvas-height", "300",
            "--structure-width", "120", "--structure-height", "80"]


def _views(root, run):
    return ["canonical-views", "--run-dir", str(run), "--chimerax-session", str(root / "scene.cxs"),
            "--map-frame", "raw", "--width", "320", "--height", "240"]


def _refuses_preserving(args, paths, entry=main):
    before = {path: path.read_bytes() for path in paths}
    try:
        status = entry([str(value) for value in args])
    except SystemExit as error:
        status = error.code
    assert all(path.exists() and path.read_bytes() == contents for path, contents in before.items())
    assert status == 2


@pytest.mark.parametrize("case", ["raw-data", "canonical-frame", "animate-basis", "views-basis", "axis-limits",
                                  "canonical-input", "canonical-use-frame"])
def test_normal_entries_refuse_outputs_containing_inputs(bundle, case):
    root, run, config = bundle
    output = config
    protected = config / "basis.json"
    if case in {"canonical-frame", "views-basis"}:
        args = _views(root, run)
        if case == "canonical-frame":
            output = run / "canonical/default"
            protected = output / "canonical_frame.json"
        else:
            args += ["--map-frame", "explicit", "--map-frame-transform", protected]
    elif case.startswith("canonical-"):
        output = run / "canonical/default"
        if case == "canonical-input":
            protected = output / "canonical_landscape.npz"
            args = ["canonicalize", "--landscape", protected, "--fit-top", ".6", "--no-visualize"]
        else:
            protected = output / "canonical_frame.json"
            args = ["canonicalize", "--run-dir", run, "--use-frame", protected, "--no-visualize"]
    else:
        args = _animate(root, run)
        if case == "raw-data":
            output = run / "data"
            protected = output / "raw_landscape.npz"
        elif case == "animate-basis":
            args += ["--map-frame", "explicit", "--map-frame-transform", protected]
        else:
            protected = config / "limits.json"
            args += ["--axis-limits", protected]
    _refuses_preserving([*args, "--output-dir", output, "--overwrite"], [protected])


@pytest.mark.parametrize("command", ["animate", "canonical-views", "canonicalize"])
def test_referenced_frame_npz_is_protected_too(bundle, command):
    root, run, config = bundle
    frame = run / "canonical/default/canonical_frame.json"
    payload = json.loads(frame.read_text())
    transform = payload.pop("canonical_transform")
    source = config / "canonical_frame.npz"
    np.savez(source, canonical_transform=transform)
    payload["canonical_frame_npz"] = str(source)
    frame.write_text(json.dumps(payload))
    if command == "animate":
        args = [*_animate(root, run), "--coordinate-set", "canonical"]
    elif command == "canonical-views":
        args = _views(root, run)
    else:
        args = ["canonicalize", "--run-dir", run, "--use-frame", frame, "--no-visualize"]
    _refuses_preserving([*args, "--output-dir", config, "--overwrite"], [source, frame])


@pytest.mark.parametrize("name,flag", [("path.csv", "--path-csv"), ("scene.cxs", "--chimerax-session")])
def test_animation_preserves_lexical_source_symlink(bundle, name, flag):
    root, run, config = bundle
    alias = config / name
    alias.symlink_to(root / name)
    _refuses_preserving([*_animate(root, run), flag, alias, "--output-dir", config, "--overwrite"],
                       [alias, root / name])
    assert alias.is_symlink()


def test_visualization_parent_symlink_cannot_redirect_overwrite_into_raw_data(bundle):
    root, run, _ = bundle
    parent = run / "visualizations"
    parent.mkdir(exist_ok=True)
    (parent / "raw").symlink_to(run, target_is_directory=True)
    _refuses_preserving(["visualize", "--run-dir", run, "--visual-id", "data", "--overwrite", "--all"],
                       [run / "data/raw_landscape.npz"])


def test_rotate_protects_explicit_input_frame(bundle):
    root, run, _ = bundle
    output = root / "rotated"
    frame = output / "canonical/default/canonical_frame.json"
    frame.parent.mkdir(parents=True)
    frame.write_bytes((run / "canonical/default/canonical_frame.json").read_bytes())
    _refuses_preserving(["--input", run / "data/raw_landscape.npz", "--canonical-frame", frame,
                        "--rotation-euler", "10", "5", "-7", "--output-dir", output, "--overwrite"],
                       [frame], entry=rotate_main)


def test_zero_opacity_2d_is_a_valid_public_render(bundle):
    _, run, _ = bundle
    assert main(["visualize", "--run-dir", str(run), "--opacity", "0", "--all"]) == 0
    output = run / "visualizations/raw/default"
    report = json.loads((output / "visualization_report.json").read_text())
    assert report["point_alpha"] == 0
    assert (output / "euler_3view_projection.png").is_file()


@pytest.mark.parametrize("command", ["canonicalize", "rotate"])
@pytest.mark.parametrize("link_kind", ["symlink", "hardlink"])
def test_existing_output_file_alias_cannot_modify_external_source(bundle, command, link_kind):
    root, run, _ = bundle
    output = root / "aliased-output"
    if command == "canonicalize":
        source = run / "data/raw_landscape.npz"
        alias = output / "canonical_landscape.npz"
        args = ["canonicalize", "--run-dir", run, "--no-visualize"]
        entry = main
    else:
        source = run / "canonical/default/canonical_frame.json"
        alias = output / "canonical/default/canonical_frame.json"
        args = ["--input", run / "data/raw_landscape.npz", "--canonical-frame", source,
                "--rotation-euler", "10", "5", "-7"]
        entry = rotate_main
    alias.parent.mkdir(parents=True)
    if link_kind == "symlink":
        alias.symlink_to(source)
    else:
        os.link(source, alias)
    _refuses_preserving([*args, "--output-dir", output, "--overwrite"], [source], entry=entry)


@pytest.mark.parametrize("command", ["animate", "canonical-views", "canonicalize"])
def test_separate_derived_output_can_be_overwritten(bundle, command):
    root, run, _ = bundle
    before = {path: path.read_bytes() for path in run.rglob("*") if path.is_file()}
    output = run / "derived" / command
    args = (_animate(root, run) if command == "animate" else _views(root, run) if command == "canonical-views"
            else ["canonicalize", "--run-dir", run, "--no-visualize"])
    args = [str(value) for value in [*args, "--output-dir", output]]
    assert main(args) == 0
    assert main([*args, "--overwrite"]) == 0
    assert before == {path: path.read_bytes() for path in before}
