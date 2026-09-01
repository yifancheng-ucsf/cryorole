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
    normalize_model_id_groups,
    validate_matching_session_baselines,
    validate_chimerax_completion,
    validate_structure_frames,
    write_scene_transforms_csv,
)


class _HierarchyPlace:
    def __init__(self, matrix=None):
        self.matrix = np.asarray(
            matrix if matrix is not None else np.eye(4)[:3],
            dtype=float,
        )

    def __mul__(self, other):
        left = np.eye(4)
        right = np.eye(4)
        left[:3] = self.matrix
        right[:3] = other.matrix
        return _HierarchyPlace((left @ right)[:3])


class _HierarchyModel:
    def __init__(self, model_id, *, parent=None):
        self.id_string = model_id
        self.name = model_id
        self.parent = parent
        self._local_position = _HierarchyPlace()

    @property
    def scene_position(self):
        if self.parent is None:
            return self._local_position
        return self.parent.scene_position * self._local_position

    @scene_position.setter
    def scene_position(self, value):
        if self.parent is None:
            self._local_position = value
            return
        parent = np.eye(4)
        target = np.eye(4)
        parent[:3] = self.parent.scene_position.matrix
        target[:3] = value.matrix
        self._local_position = _HierarchyPlace((np.linalg.inv(parent) @ target)[:3])


def _install_fake_chimerax(monkeypatch, run_command):
    chimerax = types.ModuleType("chimerax")
    chimerax_core = types.ModuleType("chimerax.core")
    chimerax_commands = types.ModuleType("chimerax.core.commands")
    chimerax_geometry = types.ModuleType("chimerax.geometry")
    chimerax_commands.run = run_command
    chimerax_geometry.Place = _HierarchyPlace
    monkeypatch.setitem(sys.modules, "chimerax", chimerax)
    monkeypatch.setitem(sys.modules, "chimerax.core", chimerax_core)
    monkeypatch.setitem(sys.modules, "chimerax.core.commands", chimerax_commands)
    monkeypatch.setitem(sys.modules, "chimerax.geometry", chimerax_geometry)


def test_model_id_groups_normalize_and_reject_duplicates_or_overlap():
    assert normalize_model_id_groups(["#1", "2"], ["#3"]) == (
        ("1", "2"),
        ("3",),
    )
    with pytest.raises(ValueError, match="duplicate reference"):
        normalize_model_id_groups(["#1", "1"], ["#2"])
    with pytest.raises(ValueError, match="both reference and moving"):
        normalize_model_id_groups(["#1"], ["1"])


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


def test_generated_script_supports_rigid_multi_model_groups(tmp_path):
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
        reference_model_id=["#1", "#2"],
        moving_model_id=["#3", "#4"],
        structure_width=640,
        structure_height=480,
    )
    script = paths.python_script.read_text(encoding="utf-8")
    config_line = next(line for line in script.splitlines() if line.startswith("CONFIG ="))
    config_text = config_line.removeprefix("CONFIG = json.loads(").removesuffix(")")
    config = json.loads(ast.literal_eval(config_text))
    assert config["reference_model_ids"] == ["1", "2"]
    assert config["moving_model_ids"] == ["3", "4"]
    assert "for moving, moving_baseline, _matrix in moving_baselines" in script
    assert "moving_group_relative_transforms_preserved" in script
    assert "stationary_unchanged_by_model" in script


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


def test_tertiary_view_script_uses_layered_paths(tmp_path):
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
        view_name="tertiary",
    )
    assert paths.python_script == tmp_path / "chimerax" / "tertiary" / "render_structure.py"
    assert paths.cxc_script == tmp_path / "chimerax" / "tertiary" / "run_render.cxc"
    assert paths.status_path == tmp_path / "logs" / "chimerax_render_tertiary_status.json"


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


def test_session_baseline_parity_validates_tertiary_transform():
    matrix = np.eye(4)[:3].tolist()
    primary = {
        "reference_baseline_matrix": matrix,
        "moving_baseline_matrix": matrix,
    }
    result = validate_matching_session_baselines(
        primary,
        dict(primary),
        dict(primary),
    )
    assert result["matched"] is True
    assert result["view_count"] == 3
    changed = json.loads(json.dumps(primary))
    changed["reference_baseline_matrix"][1][3] = 1.0
    with pytest.raises(RuntimeError, match="tertiary.*initial scene transforms differ"):
        validate_matching_session_baselines(primary, dict(primary), changed)


def test_session_baseline_parity_compares_every_group_model():
    identity = np.eye(4)[:3].tolist()
    shifted = np.eye(4)[:3]
    shifted[0, 3] = 3.0
    primary = {
        "reference_model_ids": ["1", "2"],
        "moving_model_ids": ["3", "4"],
        "reference_baseline_matrices": {"1": identity, "2": shifted.tolist()},
        "moving_baseline_matrices": {"3": identity, "4": shifted.tolist()},
    }
    result = validate_matching_session_baselines(primary, json.loads(json.dumps(primary)))
    assert result["reference_model_ids"] == ["1", "2"]
    assert result["moving_model_ids"] == ["3", "4"]
    changed = json.loads(json.dumps(primary))
    changed["moving_baseline_matrices"]["4"][1][3] = 2.0
    with pytest.raises(RuntimeError, match=r"moving model #4"):
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


def test_multi_model_completion_requires_per_model_checks(tmp_path):
    status_path = tmp_path / "status.json"
    identity = np.eye(4)[:3].tolist()
    payload = {
        "status": "success",
        "frame_count": 2,
        "reference_unchanged": True,
        "moving_restored": True,
        "stationary_unchanged": True,
        "moving_group_relative_transforms_preserved": True,
        "reference_model_ids": ["1", "2"],
        "moving_model_ids": ["3"],
        "reference_baseline_matrices": {"1": identity, "2": identity},
        "moving_baseline_matrices": {"3": identity},
        "reference_unchanged_by_model": {"1": True, "2": True},
        "stationary_unchanged_by_model": {"1": True, "2": True},
        "moving_restored_by_model": {"3": True},
        "moving_group_relative_transform_checks": {},
    }
    status_path.write_text(json.dumps(payload), encoding="utf-8")
    result = validate_chimerax_completion(
        status_path,
        expected_count=2,
        expected_reference_model_ids=["#1", "#2"],
        expected_moving_model_ids=["#3"],
    )
    assert result["reference_unchanged_by_model"]["2"] is True
    payload["moving_restored_by_model"]["3"] = False
    status_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(RuntimeError, match=r"moving model #3"):
        validate_chimerax_completion(
            status_path,
            expected_count=2,
            expected_reference_model_ids=["#1", "#2"],
            expected_moving_model_ids=["#3"],
        )


def test_multi_moving_completion_requires_every_relative_transform_check(tmp_path):
    status_path = tmp_path / "status.json"
    identity = np.eye(4)[:3].tolist()
    payload = {
        "status": "success",
        "frame_count": 1,
        "reference_unchanged": True,
        "stationary_unchanged": True,
        "moving_restored": True,
        "moving_group_relative_transforms_preserved": True,
        "reference_model_ids": ["1"],
        "moving_model_ids": ["2", "3"],
        "reference_baseline_matrices": {"1": identity},
        "moving_baseline_matrices": {"2": identity, "3": identity},
        "reference_unchanged_by_model": {"1": True},
        "stationary_unchanged_by_model": {"1": True},
        "moving_restored_by_model": {"2": True, "3": True},
        "moving_group_relative_transform_checks": {},
    }
    status_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(RuntimeError, match="every moving-group relative"):
        validate_chimerax_completion(
            status_path,
            expected_count=1,
            expected_reference_model_ids=["#1"],
            expected_moving_model_ids=["#2", "#3"],
        )


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


def test_renderer_moves_two_models_as_one_rigid_group_and_restores_them(
    tmp_path,
    monkeypatch,
):
    session_path = tmp_path / "scene.cxs"
    session_path.write_text("placeholder", encoding="utf-8")
    delta = np.eye(4)
    delta[0, 3] = 5.0
    transform_csv = write_scene_transforms_csv(
        delta[None, :, :],
        tmp_path / "transforms.csv",
    )
    scripts = generate_chimerax_scripts(
        output_dir=tmp_path,
        session_path=session_path,
        transform_csv=transform_csv,
        reference_model_id=["#1", "#2"],
        moving_model_id=["#3", "#4"],
        structure_width=32,
        structure_height=32,
    )

    class DummyPlace:
        def __init__(self, matrix=None):
            self.matrix = np.asarray(
                matrix if matrix is not None else np.eye(4)[:3],
                dtype=float,
            )

        def __mul__(self, other):
            left = np.eye(4)
            right = np.eye(4)
            left[:3] = self.matrix
            right[:3] = other.matrix
            return DummyPlace((left @ right)[:3])

    class DummyModel:
        def __init__(self, model_id, translation):
            self.id_string = model_id
            self.name = model_id
            matrix = np.eye(4)[:3]
            matrix[:, 3] = translation
            self.scene_position = DummyPlace(matrix)

    models = [
        DummyModel("1", [0, 0, 0]),
        DummyModel("2", [0, 1, 0]),
        DummyModel("3", [1, 0, 0]),
        DummyModel("4", [2, 3, 0]),
        DummyModel("5", [0, 0, 4]),
    ]
    original = {model.id_string: model.scene_position.matrix.copy() for model in models}
    rendered = []
    fake_session = types.SimpleNamespace(
        models=types.SimpleNamespace(list=lambda: models),
    )

    def fake_run(session, command):
        if command.startswith("save "):
            rendered.append(
                {model.id_string: model.scene_position.matrix.copy() for model in models}
            )

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

    exec(
        compile(
            scripts.python_script.read_text(encoding="utf-8"),
            str(scripts.python_script),
            "exec",
        ),
        {"session": fake_session},
    )

    assert len(rendered) == 1
    assert np.allclose(rendered[0]["3"][:, 3], [6, 0, 0])
    assert np.allclose(rendered[0]["4"][:, 3], [7, 3, 0])
    assert np.allclose(rendered[0]["1"], original["1"])
    assert np.allclose(rendered[0]["2"], original["2"])
    assert np.allclose(rendered[0]["5"], original["5"])
    assert all(np.allclose(model.scene_position.matrix, original[model.id_string]) for model in models)
    status = json.loads(scripts.status_path.read_text(encoding="utf-8"))
    assert status["moving_group_relative_transforms_preserved"] is True
    assert status["moving_restored_by_model"] == {"3": True, "4": True}


def test_renderer_excludes_moving_descendants_from_stationary_checks(
    tmp_path,
    monkeypatch,
):
    session_path = tmp_path / "scene.cxs"
    session_path.write_text("placeholder", encoding="utf-8")
    delta = np.eye(4)
    delta[0, 3] = 5.0
    transform_csv = write_scene_transforms_csv(
        delta[None, :, :],
        tmp_path / "transforms.csv",
    )
    scripts = generate_chimerax_scripts(
        output_dir=tmp_path,
        session_path=session_path,
        transform_csv=transform_csv,
        reference_model_id=["#1", "#2"],
        moving_model_id="#3",
        structure_width=32,
        structure_height=32,
    )
    reference_one = _HierarchyModel("1")
    reference_two = _HierarchyModel("2")
    moving = _HierarchyModel("3")
    models = [
        reference_one,
        _HierarchyModel("1.1", parent=reference_one),
        reference_two,
        _HierarchyModel("2.1", parent=reference_two),
        moving,
        _HierarchyModel("3.1", parent=moving),
    ]
    original = {model.id_string: model.scene_position.matrix.copy() for model in models}
    rendered = []
    fake_session = types.SimpleNamespace(
        models=types.SimpleNamespace(list=lambda: models),
    )

    def fake_run(session, command):
        if command.startswith("save "):
            rendered.append(
                {model.id_string: model.scene_position.matrix.copy() for model in models}
            )

    _install_fake_chimerax(monkeypatch, fake_run)
    exec(
        compile(
            scripts.python_script.read_text(encoding="utf-8"),
            str(scripts.python_script),
            "exec",
        ),
        {"session": fake_session},
    )

    assert np.allclose(rendered[0]["3"][:, 3], [5, 0, 0])
    assert np.allclose(rendered[0]["3.1"][:, 3], [5, 0, 0])
    assert all(
        np.allclose(model.scene_position.matrix, original[model.id_string])
        for model in models
    )
    status = json.loads(scripts.status_path.read_text(encoding="utf-8"))
    assert status["stationary_unchanged_by_model"] == {
        "1": True,
        "1.1": True,
        "2": True,
        "2.1": True,
    }


def test_renderer_still_rejects_true_stationary_descendant_drift(
    tmp_path,
    monkeypatch,
):
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
        moving_model_id="#3",
        structure_width=32,
        structure_height=32,
    )
    reference = _HierarchyModel("1")
    stationary = _HierarchyModel("2")
    moving = _HierarchyModel("3")
    stationary_child = _HierarchyModel("2.1", parent=stationary)
    models = [
        reference,
        _HierarchyModel("1.1", parent=reference),
        stationary,
        stationary_child,
        moving,
        _HierarchyModel("3.1", parent=moving),
    ]
    fake_session = types.SimpleNamespace(
        models=types.SimpleNamespace(list=lambda: models),
    )

    def fake_run(session, command):
        if command.startswith("save "):
            changed = stationary_child.scene_position.matrix.copy()
            changed[0, 3] = 1.0
            stationary_child.scene_position = _HierarchyPlace(changed)

    _install_fake_chimerax(monkeypatch, fake_run)
    with pytest.raises(RuntimeError, match="undeclared stationary"):
        exec(
            compile(
                scripts.python_script.read_text(encoding="utf-8"),
                str(scripts.python_script),
                "exec",
            ),
            {"session": fake_session},
        )
    status = json.loads(scripts.status_path.read_text(encoding="utf-8"))
    assert status["reference_unchanged"] is True
    assert status["stationary_unchanged_by_model"]["2.1"] is False


@pytest.mark.parametrize(
    ("reference_ids", "moving_ids", "message"),
    [
        (["#1", "#3.1"], ["#3"], "reference model #3.1"),
        (["#1"], ["#3", "#3.1"], "both parent #3 and descendant #3.1"),
    ],
)
def test_renderer_rejects_conflicting_parent_child_declarations(
    tmp_path,
    monkeypatch,
    reference_ids,
    moving_ids,
    message,
):
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
        reference_model_id=reference_ids,
        moving_model_id=moving_ids,
        structure_width=32,
        structure_height=32,
    )
    reference = _HierarchyModel("1")
    moving = _HierarchyModel("3")
    models = [reference, moving, _HierarchyModel("3.1", parent=moving)]
    fake_session = types.SimpleNamespace(
        models=types.SimpleNamespace(list=lambda: models),
    )
    _install_fake_chimerax(monkeypatch, lambda session, command: None)

    with pytest.raises(RuntimeError, match=message):
        exec(
            compile(
                scripts.python_script.read_text(encoding="utf-8"),
                str(scripts.python_script),
                "exec",
            ),
            {"session": fake_session},
        )


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
