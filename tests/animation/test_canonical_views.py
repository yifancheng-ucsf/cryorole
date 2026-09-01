from __future__ import annotations

import json
import subprocess

import numpy as np
import pytest
from PIL import Image
from scipy.spatial.transform import Rotation

from cryorole.animation.canonical_views import (
    generate_canonical_view_scripts,
    resolve_canonical_view_bases,
)
from cryorole.cli.main import build_parser, main


def _fixture_files(tmp_path):
    run_dir = tmp_path / "run"
    frame_dir = run_dir / "canonical" / "default"
    frame_dir.mkdir(parents=True)
    transform = Rotation.from_euler("z", 30, degrees=True).as_matrix()
    (frame_dir / "canonical_frame.json").write_text(
        json.dumps(
            {
                "canonical_transform": transform.tolist(),
                "transform_direction": "canonical_rv = raw_rv @ canonical_transform",
                "coordinate_space": "rotvec_ro_radians",
                "source_canonical_id": "default",
            }
        ),
        encoding="utf-8",
    )
    session = tmp_path / "scene.cxs"
    session.write_bytes(b"source-session")
    return run_dir, session, transform


def _arguments(
    run_dir,
    session,
    output_dir,
    *,
    render_mode="script-only",
    save_sessions=False,
):
    arguments = [
        "canonical-views",
        "--run-dir",
        str(run_dir),
        "--canonical-id",
        "default",
        "--chimerax-session",
        str(session),
        "--map-frame",
        "raw",
        "--render-mode",
        render_mode,
        "--width",
        "320",
        "--height",
        "240",
        "--output-dir",
        str(output_dir),
    ]
    if save_sessions:
        arguments.append("--save-sessions")
    return arguments


def test_canonical_view_bases_follow_recorded_row_vector_contract():
    canonical = Rotation.from_euler("z", 30, degrees=True).as_matrix()
    scene = Rotation.from_euler("x", 20, degrees=True).as_matrix()
    views = resolve_canonical_view_bases(canonical, scene)
    axes = scene @ canonical

    assert [view.file_name for view in views] == [
        "canonical_x_plus.png",
        "canonical_y_plus.png",
        "canonical_z_plus.png",
    ]
    assert np.allclose(views[0].outward, axes[:, 0])
    assert np.allclose(views[0].right, axes[:, 1])
    assert np.allclose(views[0].up, axes[:, 2])
    assert np.allclose(views[1].outward, axes[:, 1])
    assert np.allclose(views[1].right, axes[:, 0])
    assert np.allclose(views[1].up, -axes[:, 2])
    assert np.allclose(views[2].outward, axes[:, 2])
    assert np.allclose(views[2].right, axes[:, 0])
    assert np.allclose(views[2].up, axes[:, 1])
    assert all(np.linalg.det(view.camera_rotation) == pytest.approx(1.0) for view in views)


def test_generated_script_changes_only_camera_and_restores_it(tmp_path):
    session = tmp_path / "scene.cxs"
    session.write_bytes(b"session")
    views = resolve_canonical_view_bases(np.eye(3), np.eye(3))
    scripts = generate_canonical_view_scripts(
        output_dir=tmp_path,
        session_path=session,
        views=views,
        width=320,
        height=240,
        save_sessions=True,
    )

    script = scripts.python_script.read_text(encoding="utf-8")
    assert "camera.position = camera_place" in script
    assert "camera.position = camera_baseline" in script
    assert "model.scene_position =" not in script
    assert "models_unchanged" in script
    assert "camera_restored" in script
    assert "canonical_y_plus.png" in script
    assert "canonical_y_plus.cxs" in script
    assert 'CONFIG["save_sessions"]' in script
    compile(script, str(scripts.python_script), "exec")
    assert scripts.cxc_script.is_file()


def test_canonical_views_cli_script_only_is_non_destructive(tmp_path):
    parser = build_parser()
    parsed = parser.parse_args(
        [
            "canonical-views",
            "--run-dir",
            "run",
            "--chimerax-session",
            "scene.cxs",
            "--map-frame",
            "raw",
            "--output-dir",
            "out",
        ]
    )
    assert parsed.render_mode == "script-only"
    assert parsed.width == parsed.height == 900
    assert parsed.save_sessions is False

    run_dir, session, transform = _fixture_files(tmp_path)
    before_frame = (run_dir / "canonical" / "default" / "canonical_frame.json").read_bytes()
    output = tmp_path / "views"
    assert main(_arguments(run_dir, session, output)) == 0

    manifest = json.loads((output / "canonical_views.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "scripts_ready"
    assert manifest["map_frame"] == "raw"
    assert np.allclose(manifest["canonical_transform"], transform)
    assert len(manifest["views"]) == 3
    assert (output / "chimerax" / "render_canonical_views.py").is_file()
    assert not (output / "views").exists()
    assert session.read_bytes() == b"source-session"
    assert (
        run_dir / "canonical" / "default" / "canonical_frame.json"
    ).read_bytes() == before_frame


def test_canonical_views_save_sessions_flag_is_script_only_safe(tmp_path):
    run_dir, session, _ = _fixture_files(tmp_path)
    output = tmp_path / "views"

    assert main(_arguments(run_dir, session, output, save_sessions=True)) == 0

    manifest = json.loads((output / "canonical_views.json").read_text(encoding="utf-8"))
    script = (output / "chimerax" / "render_canonical_views.py").read_text(
        encoding="utf-8"
    )
    assert manifest["status"] == "scripts_ready"
    assert manifest["save_sessions"] is True
    assert "canonical_x_plus.cxs" in script
    assert not (output / "sessions").exists()
    assert session.read_bytes() == b"source-session"


def test_canonical_views_execute_validates_completion_and_pngs(tmp_path, monkeypatch):
    run_dir, session, _ = _fixture_files(tmp_path)
    output = tmp_path / "views"
    executable = tmp_path / "chimerax"
    executable.write_bytes(b"placeholder")

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        view_dir = output / "views"
        view_dir.mkdir(parents=True)
        for name in (
            "canonical_x_plus.png",
            "canonical_y_plus.png",
            "canonical_z_plus.png",
        ):
            Image.new("RGB", (320, 240), "white").save(view_dir / name)
        status = output / "logs" / "chimerax_render_status.json"
        status.parent.mkdir(parents=True, exist_ok=True)
        status.write_text(
            json.dumps(
                {
                    "status": "success",
                    "frame_count": 3,
                    "models_unchanged": True,
                    "camera_restored": True,
                }
            ),
            encoding="utf-8",
        )
        log_path.write_text("{}", encoding="utf-8")
        return subprocess.CompletedProcess([], 0, "", "")

    monkeypatch.setattr(
        "cryorole.animation.canonical_views.execute_chimerax",
        fake_execute,
    )
    arguments = _arguments(run_dir, session, output, render_mode="execute")
    arguments.extend(["--chimerax-bin", str(executable)])
    assert main(arguments) == 0

    manifest = json.loads((output / "canonical_views.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "rendered"
    assert manifest["completion"]["models_unchanged"] is True
    assert manifest["completion"]["camera_restored"] is True
    assert manifest["image_validation"]["dimensions"] == [320, 240]
    assert len(list((output / "views").glob("*.png"))) == 3


def test_canonical_views_execute_saves_and_validates_camera_sessions(tmp_path, monkeypatch):
    run_dir, session, _ = _fixture_files(tmp_path)
    output = tmp_path / "views"
    executable = tmp_path / "chimerax"
    executable.write_bytes(b"placeholder")

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        view_dir = output / "views"
        session_dir = output / "sessions"
        view_dir.mkdir(parents=True)
        session_dir.mkdir(parents=True)
        names = ("canonical_x_plus", "canonical_y_plus", "canonical_z_plus")
        for name in names:
            Image.new("RGB", (320, 240), "white").save(view_dir / f"{name}.png")
            (session_dir / f"{name}.cxs").write_bytes(b"derived-session")
        status = output / "logs" / "chimerax_render_status.json"
        status.parent.mkdir(parents=True, exist_ok=True)
        status.write_text(
            json.dumps(
                {
                    "status": "success",
                    "frame_count": 3,
                    "session_count": 3,
                    "sessions_saved": True,
                    "session_files": [str(session_dir / f"{name}.cxs") for name in names],
                    "models_unchanged": True,
                    "camera_restored": True,
                }
            ),
            encoding="utf-8",
        )
        log_path.write_text("{}", encoding="utf-8")
        return subprocess.CompletedProcess([], 0, "", "")

    monkeypatch.setattr(
        "cryorole.animation.canonical_views.execute_chimerax",
        fake_execute,
    )
    arguments = _arguments(
        run_dir,
        session,
        output,
        render_mode="execute",
        save_sessions=True,
    )
    arguments.extend(["--chimerax-bin", str(executable)])
    assert main(arguments) == 0

    manifest = json.loads((output / "canonical_views.json").read_text(encoding="utf-8"))
    assert manifest["completion"]["sessions_saved"] is True
    assert manifest["session_validation"]["count"] == 3
    assert manifest["output_paths"]["sessions"] == str(output / "sessions")
    assert session.read_bytes() == b"source-session"


def test_canonical_views_zero_return_without_completion_marks_failed(tmp_path, monkeypatch):
    run_dir, session, _ = _fixture_files(tmp_path)
    output = tmp_path / "views"
    executable = tmp_path / "chimerax"
    executable.write_bytes(b"placeholder")

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("{}", encoding="utf-8")
        return subprocess.CompletedProcess([], 0, "", "")

    monkeypatch.setattr(
        "cryorole.animation.canonical_views.execute_chimerax",
        fake_execute,
    )
    arguments = _arguments(run_dir, session, output, render_mode="execute")
    arguments.extend(["--chimerax-bin", str(executable)])
    with pytest.raises(SystemExit) as error:
        main(arguments)
    assert error.value.code == 2
    manifest = json.loads((output / "canonical_views.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert manifest["stage_statuses"]["render"] == "failed"
    assert session.read_bytes() == b"source-session"
