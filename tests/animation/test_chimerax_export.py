from __future__ import annotations

import ast
import csv
import json
import subprocess
import sys
import types
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from cryorole.animation.chimerax import (
    execute_chimerax,
    generate_chimerax_scripts,
    validate_matching_session_baselines,
    validate_chimerax_completion,
    validate_structure_frames,
    write_scene_transforms_csv,
)


def test_generated_script_uses_exact_ids_and_absolute_left_composition(tmp_path):
    session = tmp_path / "scene.cxs"
    session.write_text("placeholder", encoding="utf-8")
    transforms = np.repeat(np.eye(4)[None, :, :], 2, axis=0)
    transform_csv = write_scene_transforms_csv(
        transforms,
        tmp_path / "chimerax" / "frame_transforms.csv",
    )
    paths = generate_chimerax_scripts(
        output_dir=tmp_path,
        session_path=session,
        transform_csv=transform_csv,
        reference_model_id="#1",
        moving_model_id="#2.1",
        structure_width=640,
        structure_height=480,
    )
    script = paths.python_script.read_text(encoding="utf-8")
    assert "model.id_string == requested" in script
    assert "delta_place * moving_baseline" in script
    assert "moving.scene_position = delta_place * moving_baseline" in script
    assert "reference scene transform changed" in script
    assert "moving.scene_position = moving_baseline" in script
    assert "chimerax_render_status.json" in script
    assert '"status": "failed"' in script
    assert '"status": "success"' in script
    assert "scene.cxs" in script
    compile(script, str(paths.python_script), "exec")
    assert paths.cxc_script.is_file()


def test_dual_view_scripts_use_layered_paths_and_record_initial_transforms(tmp_path):
    session = tmp_path / "scene.cxs"
    session.write_text("placeholder", encoding="utf-8")
    transform_csv = write_scene_transforms_csv(
        np.eye(4)[None, :, :],
        tmp_path / "frame_transforms.csv",
    )
    paths = generate_chimerax_scripts(
        output_dir=tmp_path,
        session_path=session,
        transform_csv=transform_csv,
        reference_model_id="#1",
        moving_model_id="#2",
        structure_width=640,
        structure_height=480,
        view_name="secondary",
    )
    assert paths.python_script == tmp_path / "chimerax" / "secondary" / "render_structure.py"
    assert paths.cxc_script == tmp_path / "chimerax" / "secondary" / "run_render.cxc"
    assert paths.status_path == tmp_path / "logs" / "chimerax_render_secondary_status.json"
    script = paths.python_script.read_text(encoding="utf-8")
    config_line = next(line for line in script.splitlines() if line.startswith("CONFIG ="))
    config_text = config_line.removeprefix("CONFIG = json.loads(").removesuffix(")")
    config = json.loads(ast.literal_eval(config_text))
    assert Path(config["structure_dir"]) == (
        tmp_path / "frames" / "structure" / "secondary"
    ).resolve()
    assert 'STATE["reference_baseline_matrix"]' in script
    assert 'STATE["moving_baseline_matrix"]' in script


def test_session_baseline_parity_requires_both_model_transforms_to_match():
    matrix = np.eye(4)[:3].tolist()
    primary = {
        "reference_baseline_matrix": matrix,
        "moving_baseline_matrix": matrix,
    }
    result = validate_matching_session_baselines(primary, dict(primary))
    assert result["matched"] is True
    changed = json.loads(json.dumps(primary))
    changed["moving_baseline_matrix"][0][3] = 1.0
    with pytest.raises(RuntimeError, match="initial scene transforms differ"):
        validate_matching_session_baselines(primary, changed)


def test_execute_uses_argument_list_and_handles_nonzero(tmp_path, monkeypatch):
    cxc = tmp_path / "run.cxc"
    cxc.write_text("", encoding="utf-8")
    executable = tmp_path / "ChimeraX.exe"
    executable.write_bytes(b"placeholder")
    log = tmp_path / "logs" / "chimerax_render.log"
    seen = {}

    def fake_run(arguments, **kwargs):
        seen["arguments"] = arguments
        seen["kwargs"] = kwargs
        return subprocess.CompletedProcess(arguments, 9, "stdout", "stderr")

    monkeypatch.setattr(subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="code 9"):
        execute_chimerax(
            chimerax_bin=executable,
            cxc_script=cxc,
            log_path=log,
        )
    assert isinstance(seen["arguments"], list)
    assert seen["kwargs"]["shell"] is False
    assert "stdout" in log.read_text(encoding="utf-8")


def test_linux_execute_uses_offscreen_nocolor_and_script(tmp_path, monkeypatch):
    cxc = tmp_path / "run.cxc"
    cxc.write_text("", encoding="utf-8")
    executable = tmp_path / "chimerax"
    executable.write_bytes(b"placeholder")
    seen = {}

    def fake_run(arguments, **kwargs):
        seen["arguments"] = arguments
        return subprocess.CompletedProcess(arguments, 0, "", "")

    monkeypatch.setattr(subprocess, "run", fake_run)
    execute_chimerax(
        chimerax_bin=executable,
        cxc_script=cxc,
        log_path=tmp_path / "render.log",
        platform_name="linux",
    )
    assert seen["arguments"] == [
        str(executable),
        "--offscreen",
        "--exit",
        "--nocolor",
        "--script",
        str(cxc.resolve()),
    ]


def test_completion_artifact_is_required_even_after_zero_return_code(tmp_path):
    status_path = tmp_path / "logs" / "chimerax_render_status.json"
    with pytest.raises(RuntimeError, match="missing"):
        validate_chimerax_completion(status_path, expected_count=2)
    status_path.parent.mkdir(parents=True)
    status_path.write_text(
        '{"status":"failed","error":"OpenGL unavailable"}',
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="OpenGL unavailable"):
        validate_chimerax_completion(status_path, expected_count=2)
    status_path.write_text(
        '{"status":"success","frame_count":2,"reference_unchanged":true,'
        '"moving_restored":true}',
        encoding="utf-8",
    )
    payload = validate_chimerax_completion(status_path, expected_count=2)
    assert payload["status"] == "success"


def test_renderer_exception_writes_failed_completion_status(tmp_path, monkeypatch):
    session_path = tmp_path / "scene.cxs"
    session_path.write_text("placeholder", encoding="utf-8")
    transform_csv = write_scene_transforms_csv(
        np.eye(4)[None, :, :],
        tmp_path / "transforms.csv",
    )
    scripts = generate_chimerax_scripts(
        output_dir=tmp_path,
        session_path=session_path,
        transform_csv=transform_csv,
        reference_model_id="#1",
        moving_model_id="#2",
        structure_width=32,
        structure_height=32,
    )

    class DummyPlace:
        def __init__(self, matrix=None):
            self.matrix = np.asarray(matrix if matrix is not None else np.eye(4)[:3], dtype=float)

        def __mul__(self, other):
            return DummyPlace(self.matrix)

    class DummyModel:
        def __init__(self, model_id):
            self.id_string = model_id
            self.name = model_id
            self.scene_position = DummyPlace()

    models = [DummyModel("1"), DummyModel("2")]
    fake_session = types.SimpleNamespace(
        models=types.SimpleNamespace(list=lambda: models),
    )

    def fake_run(session, command):
        if command.startswith("save "):
            raise RuntimeError("synthetic renderer failure")

    chimerax = types.ModuleType("chimerax")
    chimerax_core = types.ModuleType("chimerax.core")
    chimerax_commands = types.ModuleType("chimerax.core.commands")
    chimerax_geometry = types.ModuleType("chimerax.geometry")
    chimerax_commands.run = fake_run
    chimerax_geometry.Place = DummyPlace
    monkeypatch.setitem(sys.modules, "chimerax", chimerax)
    monkeypatch.setitem(sys.modules, "chimerax.core", chimerax_core)
    monkeypatch.setitem(sys.modules, "chimerax.core.commands", chimerax_commands)
    monkeypatch.setitem(sys.modules, "chimerax.geometry", chimerax_geometry)

    namespace = {"session": fake_session}
    with pytest.raises(RuntimeError, match="synthetic renderer failure"):
        exec(
            compile(
                scripts.python_script.read_text(encoding="utf-8"),
                str(scripts.python_script),
                "exec",
            ),
            namespace,
        )
    status = json.loads(scripts.status_path.read_text(encoding="utf-8"))
    assert status["status"] == "failed"
    assert "synthetic renderer failure" in status["error"]
    assert status["moving_restored"] is True


def test_structure_frame_validation(tmp_path):
    directory = tmp_path / "frames"
    directory.mkdir()
    for index in range(2):
        plt.imsave(directory / f"frame_{index:06d}.png", np.zeros((8, 12, 3)))
    result = validate_structure_frames(directory, expected_count=2)
    assert result.frame_count == 2
    assert result.dimensions == (12, 8)
    (directory / "frame_000001.png").unlink()
    (directory / "frame_000002.png").write_bytes(b"bad")
    with pytest.raises(RuntimeError, match="indices"):
        validate_structure_frames(directory, expected_count=2)


def test_transform_csv_has_contiguous_indices(tmp_path):
    output = write_scene_transforms_csv(
        np.repeat(np.eye(4)[None, :, :], 3, axis=0),
        tmp_path / "transforms.csv",
    )
    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [int(row["frame_index"]) for row in rows] == [0, 1, 2]
    assert "m03" in rows[0]
