from __future__ import annotations

import csv
import json
import subprocess
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pytest

from cryorole.cli.main import build_parser, main


def _fixture_files(tmp_path):
    run_dir = tmp_path / "run"
    data_dir = run_dir / "data"
    data_dir.mkdir(parents=True)
    coordinates = np.asarray([[0, 0, 0], [0.1, 0.2, 0.3], [-0.2, 0.1, 0.4]])
    np.savez(
        data_dir / "raw_landscape.npz",
        artifact_type=np.asarray("raw_landscape"),
        schema_version=np.asarray("1"),
        particle_key=np.asarray(["a", "b", "c"]),
        coordinates_analysis=coordinates,
        sld_display=np.asarray([1.0, 2.0, 3.0]),
        sld_display_is_outlier=np.asarray([False, False, False]),
    )
    (run_dir / "run_summary.json").write_text(
        json.dumps({"euler_convention": "extrinsic_zyx"}),
        encoding="utf-8",
    )
    path_csv = tmp_path / "path.csv"
    with path_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["label", "rv_x_rad", "rv_y_rad", "rv_z_rad"],
        )
        writer.writeheader()
        writer.writerows(
            [
                {"label": "start", "rv_x_rad": 0, "rv_y_rad": 0, "rv_z_rad": 0},
                {"label": "end", "rv_x_rad": 0, "rv_y_rad": 0, "rv_z_rad": 0.5},
            ]
        )
    session = tmp_path / "scene.cxs"
    session.write_bytes(b"session")
    return run_dir, path_csv, session


def _arguments(run_dir, path_csv, session, output_dir, render_mode="script-only"):
    values = [
        "animate",
        "--run-dir",
        str(run_dir),
        "--coordinate-set",
        "raw",
        "--path-csv",
        str(path_csv),
        "--path-space",
        "rv",
        "--chimerax-session",
        str(session),
        "--reference-model-id",
        "#1",
        "--moving-model-id",
        "#2",
        "--pivot",
        "0",
        "0",
        "0",
        "--baseline-ro",
        "identity",
        "--map-frame",
        "raw",
        "--frames-per-segment",
        "3",
        "--canvas-width",
        "600",
        "--canvas-height",
        "300",
        "--structure-width",
        "120",
        "--structure-height",
        "80",
        "--render-mode",
        render_mode,
        "--output-dir",
        str(output_dir),
    ]
    return values


def test_animate_cli_registration_and_script_only_bundle(tmp_path):
    parser = build_parser()
    assert "animate" in parser.format_help()
    run_dir, path_csv, session = _fixture_files(tmp_path)
    before = {
        path.relative_to(run_dir): path.read_bytes()
        for path in run_dir.rglob("*")
        if path.is_file()
    }
    output = tmp_path / "animation"
    assert main(_arguments(run_dir, path_csv, session, output)) == 0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "scripts_ready"
    assert manifest["frame_count"] == 3
    assert manifest["quaternion_storage_order"] == "wxyz"
    assert manifest["stage_statuses"]["compositor"] == "not_requested"
    assert manifest["stage_statuses"]["encoder"] == "not_requested"
    assert manifest["output_paths"]["chimerax_render_status"] is None
    assert len(list((output / "frames" / "landscape").glob("frame_*.png"))) == 3
    assert not (output / "frames" / "structure").exists()
    assert not (output / "frames" / "composite").exists()
    assert not list(output.glob("*.mp4"))
    assert not (run_dir / "selections").exists()
    after = {
        path.relative_to(run_dir): path.read_bytes()
        for path in run_dir.rglob("*")
        if path.is_file()
    }
    assert after == before
    assert session.read_bytes() == b"session"


def test_animate_display_cli_aliases_and_style_controls():
    parser = build_parser()
    common = [
        "animate", "--run-dir", "run", "--path-csv", "path.csv", "--path-space", "ea",
        "--chimerax-session", "scene.cxs", "--reference-model-id", "#1",
        "--moving-model-id", "#2", "--pivot", "0", "0", "0", "--baseline-ro", "identity",
        "--map-frame", "raw", "--output-dir", "out",
    ]
    threshold = parser.parse_args(
        common + ["--threshold", "2", "--colormap", "viridis", "--vmin", "1",
                  "--vmax", "4", "--range", "alpha:-30:30"]
    )
    alias = parser.parse_args(common + ["--sld-threshold", "2"])
    defaults = parser.parse_args(common)
    assert threshold.threshold == alias.threshold == 2.0
    assert threshold.colormap == "viridis"
    assert threshold.range_bound == [("alpha", (-30.0, 30.0))]
    assert defaults.layout == "stacked"
    assert defaults.dual_structure_horizontal_crop == 0.0
    with pytest.raises(SystemExit):
        parser.parse_args(common + ["--threshold", "2", "--top-fraction", "0.5"])
    phase4 = parser.parse_args(
        common
        + [
            "--no-encode",
            "--composite-width",
            "1920",
            "--composite-height",
            "1080",
            "--layout",
            "side_by_side",
            "--background-color",
            "#ffffff",
        ]
    )
    assert phase4.no_encode is True
    assert phase4.composite_width == 1920
    assert phase4.layout == "side_by_side"
    dual = parser.parse_args(
        common
        + [
            "--secondary-chimerax-session",
            "scene-2.cxs",
            "--dual-structure-horizontal-crop",
            "0.10",
        ]
    )
    assert dual.secondary_chimerax_session == "scene-2.cxs"
    assert dual.dual_structure_horizontal_crop == pytest.approx(0.1)


def test_dual_view_script_only_uses_layered_paths_and_keeps_sessions_unchanged(tmp_path):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    secondary = tmp_path / "scene-secondary.cxs"
    secondary.write_bytes(b"secondary session")
    output = tmp_path / "animation"
    arguments = _arguments(run_dir, path_csv, session, output)
    arguments.extend(["--secondary-chimerax-session", str(secondary)])
    assert main(arguments) == 0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "scripts_ready"
    assert manifest["schema_version"] == "5"
    assert manifest["structure_view_count"] == 2
    assert manifest["session_transform_parity"]["status"] == "not_checked_script_only"
    assert (output / "chimerax" / "primary" / "render_structure.py").is_file()
    assert (output / "chimerax" / "secondary" / "render_structure.py").is_file()
    assert not (output / "frames" / "structure").exists()
    assert session.read_bytes() == b"session"
    assert secondary.read_bytes() == b"secondary session"


def test_dual_view_rejects_side_by_side_before_rendering(tmp_path):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    secondary = tmp_path / "scene-secondary.cxs"
    secondary.write_bytes(b"secondary session")
    arguments = _arguments(run_dir, path_csv, session, tmp_path / "animation")
    arguments.extend(
        [
            "--secondary-chimerax-session",
            str(secondary),
            "--layout",
            "side_by_side",
        ]
    )
    with pytest.raises(SystemExit):
        main(arguments)


def test_nonzero_dual_crop_rejects_single_view_before_rendering(tmp_path):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    arguments = _arguments(run_dir, path_csv, session, tmp_path / "animation")
    arguments.extend(["--dual-structure-horizontal-crop", "0.1"])
    with pytest.raises(SystemExit):
        main(arguments)


def test_dual_execute_requires_both_views_and_matching_session_transforms(
    tmp_path,
    monkeypatch,
):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    secondary = tmp_path / "scene-secondary.cxs"
    secondary.write_bytes(b"secondary session")
    output = tmp_path / "animation"
    executable = tmp_path / "ChimeraX.exe"
    executable.write_bytes(b"placeholder")

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        view = cxc_script.parent.name
        structure = output / "frames" / "structure" / view
        structure.mkdir(parents=True)
        for index in range(3):
            color = np.zeros((80, 120, 3))
            color[..., 0 if view == "primary" else 1] = 1
            plt.imsave(structure / f"frame_{index:06d}.png", color)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("{}", encoding="utf-8")
        matrix = np.eye(4)[:3].tolist()
        (output / "logs" / f"chimerax_render_{view}_status.json").write_text(
            json.dumps(
                {
                    "status": "success",
                    "frame_count": 3,
                    "reference_unchanged": True,
                    "moving_restored": True,
                    "reference_baseline_matrix": matrix,
                    "moving_baseline_matrix": matrix,
                }
            ),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess([], 0, "", "")

    monkeypatch.setattr("cryorole.animation.workflow.execute_chimerax", fake_execute)
    arguments = _arguments(run_dir, path_csv, session, output, render_mode="execute")
    arguments.extend(
        [
            "--secondary-chimerax-session",
            str(secondary),
            "--dual-structure-horizontal-crop",
            "0.10",
            "--chimerax-bin",
            str(executable),
            "--no-encode",
        ]
    )
    assert main(arguments) == 0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "composite_rendered"
    assert manifest["structure_view_count"] == 2
    assert manifest["session_transform_parity"]["matched"] is True
    assert manifest["stage_statuses"]["structure_render_primary"] == "complete"
    assert manifest["stage_statuses"]["structure_render_secondary"] == "complete"
    assert manifest["compositor"]["structure_view_count"] == 2
    assert len(manifest["compositor"]["structure_destinations"]) == 2
    assert manifest["dual_structure_horizontal_crop_fraction"] == pytest.approx(0.1)
    assert manifest["compositor"]["dual_structure_horizontal_crop_fraction"] == pytest.approx(0.1)
    assert manifest["compositor"]["structure_source_crop_rectangles"] == [
        [12, 0, 108, 80],
        [12, 0, 108, 80],
    ]
    assert manifest["compositor"]["structure_cropped_dimensions"] == [
        [96, 80],
        [96, 80],
    ]
    assert manifest["compositor"]["structure_alignments"] == ["right", "left"]
    assert len(list((output / "frames" / "composite").glob("frame_*.png"))) == 3


def test_dual_execute_transform_mismatch_prevents_composition(tmp_path, monkeypatch):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    secondary = tmp_path / "scene-secondary.cxs"
    secondary.write_bytes(b"secondary session")
    output = tmp_path / "animation"
    executable = tmp_path / "ChimeraX.exe"
    executable.write_bytes(b"placeholder")
    compositor_called = False

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        view = cxc_script.parent.name
        structure = output / "frames" / "structure" / view
        structure.mkdir(parents=True)
        for index in range(3):
            plt.imsave(
                structure / f"frame_{index:06d}.png",
                np.zeros((80, 120, 3)),
            )
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("{}", encoding="utf-8")
        reference = np.eye(4)[:3]
        moving = np.eye(4)[:3]
        if view == "secondary":
            moving = moving.copy()
            moving[0, 3] = 1
        (output / "logs" / f"chimerax_render_{view}_status.json").write_text(
            json.dumps(
                {
                    "status": "success",
                    "frame_count": 3,
                    "reference_unchanged": True,
                    "moving_restored": True,
                    "reference_baseline_matrix": reference.tolist(),
                    "moving_baseline_matrix": moving.tolist(),
                }
            ),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess([], 0, "", "")

    def fake_compose(**kwargs):
        nonlocal compositor_called
        compositor_called = True

    monkeypatch.setattr("cryorole.animation.workflow.execute_chimerax", fake_execute)
    monkeypatch.setattr(
        "cryorole.animation.workflow.compose_frame_sequences",
        fake_compose,
    )
    arguments = _arguments(run_dir, path_csv, session, output, render_mode="execute")
    arguments.extend(
        [
            "--secondary-chimerax-session",
            str(secondary),
            "--chimerax-bin",
            str(executable),
            "--no-encode",
        ]
    )
    with pytest.raises(SystemExit):
        main(arguments)
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert manifest["session_transform_parity"]["matched"] is False
    assert manifest["stage_statuses"]["compositor"] == "not_started"
    assert compositor_called is False


def test_dual_execute_secondary_completion_failure_prevents_composition(
    tmp_path,
    monkeypatch,
):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    secondary = tmp_path / "scene-secondary.cxs"
    secondary.write_bytes(b"secondary session")
    output = tmp_path / "animation"
    executable = tmp_path / "ChimeraX.exe"
    executable.write_bytes(b"placeholder")
    compositor_called = False

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        view = cxc_script.parent.name
        structure = output / "frames" / "structure" / view
        structure.mkdir(parents=True)
        for index in range(3):
            plt.imsave(
                structure / f"frame_{index:06d}.png",
                np.zeros((80, 120, 3)),
            )
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("{}", encoding="utf-8")
        if view == "primary":
            matrix = np.eye(4)[:3].tolist()
            (output / "logs" / "chimerax_render_primary_status.json").write_text(
                json.dumps(
                    {
                        "status": "success",
                        "frame_count": 3,
                        "reference_unchanged": True,
                        "moving_restored": True,
                        "reference_baseline_matrix": matrix,
                        "moving_baseline_matrix": matrix,
                    }
                ),
                encoding="utf-8",
            )
        return subprocess.CompletedProcess([], 0, "", "")

    def fake_compose(**kwargs):
        nonlocal compositor_called
        compositor_called = True

    monkeypatch.setattr("cryorole.animation.workflow.execute_chimerax", fake_execute)
    monkeypatch.setattr(
        "cryorole.animation.workflow.compose_frame_sequences",
        fake_compose,
    )
    arguments = _arguments(run_dir, path_csv, session, output, render_mode="execute")
    arguments.extend(
        [
            "--secondary-chimerax-session",
            str(secondary),
            "--chimerax-bin",
            str(executable),
            "--no-encode",
        ]
    )
    with pytest.raises(SystemExit):
        main(arguments)
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage_statuses"]["structure_render_primary"] == "complete"
    assert manifest["stage_statuses"]["structure_render_secondary"] == "failed"
    assert manifest["stage_statuses"]["compositor"] == "not_started"
    assert compositor_called is False


def test_execute_mode_updates_manifest_and_validates_frames(tmp_path, monkeypatch):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    output = tmp_path / "animation"
    executable = tmp_path / "ChimeraX.exe"
    executable.write_bytes(b"placeholder")

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        structure = output / "frames" / "structure"
        structure.mkdir(parents=True)
        for index in range(3):
            plt.imsave(structure / f"frame_{index:06d}.png", np.zeros((80, 120, 3)))
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("{}", encoding="utf-8")
        status = output / "logs" / "chimerax_render_status.json"
        status.write_text(
            json.dumps(
                {
                    "status": "success",
                    "frame_count": 3,
                    "reference_unchanged": True,
                    "moving_restored": True,
                }
            ),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess([], 0, "", "")

    monkeypatch.setattr("cryorole.animation.workflow.execute_chimerax", fake_execute)
    arguments = _arguments(run_dir, path_csv, session, output, render_mode="execute")
    arguments.extend(["--chimerax-bin", str(executable), "--no-encode"])
    assert main(arguments) == 0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "composite_rendered"
    assert manifest["ChimeraX_return_code"] == 0
    assert manifest["structure_validation"]["dimensions"] == [120, 80]
    assert manifest["stage_statuses"]["compositor"] == "complete"
    assert manifest["stage_statuses"]["encoder"] == "not_requested"
    assert manifest["compositor"]["layout"] == "stacked"
    assert manifest["compositor"]["layout_version"] == "4"
    assert manifest["compositor"]["structure_fraction"] == pytest.approx(0.64)
    assert manifest["compositor"]["landscape_fraction"] == pytest.approx(0.36)
    display_policy = manifest["resolved_visualize_display_policy"]
    assert len(display_policy["projection_panel_rectangles"]) == 3
    assert display_policy["annotation_policy"] == "coordinate_only_once_when_enabled"
    assert display_policy["shared_colorbar_policy"] == "single_horizontal_centered"
    assert display_policy["visible_colorbar_label"] == "SLD"
    assert display_policy["underlying_color_field"] == "sld_display"
    assert len(list((output / "frames" / "composite").glob("frame_*.png"))) == 3
    assert not list(output.glob("*.mp4"))


def test_execute_zero_return_without_completion_artifact_marks_failed(tmp_path, monkeypatch):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    output = tmp_path / "animation"
    executable = tmp_path / "ChimeraX.exe"
    executable.write_bytes(b"placeholder")

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("{}", encoding="utf-8")
        return subprocess.CompletedProcess([], 0, "", "")

    monkeypatch.setattr("cryorole.animation.workflow.execute_chimerax", fake_execute)
    arguments = _arguments(run_dir, path_csv, session, output, render_mode="execute")
    arguments.extend(["--chimerax-bin", str(executable), "--no-encode"])
    with pytest.raises(SystemExit) as error:
        main(arguments)
    assert error.value.code == 2
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert manifest["ChimeraX_return_code"] == 0
    assert manifest["stage_statuses"]["compositor"] == "not_started"
    assert (output / "logs" / "chimerax_render.log").is_file()
    assert session.read_bytes() == b"session"
    assert not (output / "frames" / "composite").exists()
    assert not list(output.glob("*.mp4"))


def test_execute_encodes_and_validates_movie_after_compositor(tmp_path, monkeypatch):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    output = tmp_path / "animation"
    chimerax = tmp_path / "ChimeraX.exe"
    ffmpeg = tmp_path / "ffmpeg.exe"
    ffprobe = tmp_path / "ffprobe.exe"
    for executable in (chimerax, ffmpeg, ffprobe):
        executable.write_bytes(b"placeholder")
    events = []

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        events.append("chimerax")
        structure = output / "frames" / "structure"
        structure.mkdir(parents=True)
        for index in range(3):
            plt.imsave(structure / f"frame_{index:06d}.png", np.zeros((80, 120, 3)))
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("{}", encoding="utf-8")
        (output / "logs" / "chimerax_render_status.json").write_text(
            json.dumps(
                {
                    "status": "success",
                    "frame_count": 3,
                    "reference_unchanged": True,
                    "moving_restored": True,
                }
            ),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess([], 0, "", "")

    def fake_encode(**kwargs):
        events.append("ffmpeg")
        movie = kwargs["output_path"]
        movie.write_bytes(b"mp4")
        return SimpleNamespace(
            output_path=movie,
            arguments=("ffmpeg",),
            return_code=0,
            log_path=kwargs["log_path"],
        )

    def fake_validate(**kwargs):
        events.append("ffprobe")
        return SimpleNamespace(
            movie_path=kwargs["movie_path"],
            codec="h264",
            pixel_format="yuv420p",
            dimensions=kwargs["expected_dimensions"],
            fps=kwargs["expected_fps"],
            frame_count=kwargs["expected_frame_count"],
            duration_sec=0.3,
            validation_method="reported_frame_count",
            frame_count_tolerance=0.0,
            probe_payload={"streams": []},
            arguments=("ffprobe",),
            return_code=0,
            log_path=kwargs["log_path"],
        )

    monkeypatch.setattr("cryorole.animation.workflow.execute_chimerax", fake_execute)
    monkeypatch.setattr("cryorole.animation.workflow.encode_mp4", fake_encode)
    monkeypatch.setattr("cryorole.animation.workflow.validate_mp4", fake_validate)
    arguments = _arguments(run_dir, path_csv, session, output, render_mode="execute")
    arguments.extend(
        [
            "--chimerax-bin",
            str(chimerax),
            "--ffmpeg-bin",
            str(ffmpeg),
            "--ffprobe-bin",
            str(ffprobe),
        ]
    )
    assert main(arguments) == 0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert events == ["chimerax", "ffmpeg", "ffprobe"]
    assert manifest["status"] == "movie_encoded"
    assert manifest["stage_statuses"]["compositor"] == "complete"
    assert manifest["stage_statuses"]["encoder"] == "complete"
    assert manifest["stage_statuses"]["movie_validation"] == "complete"
    assert manifest["movie_validation"]["validation_method"] == "reported_frame_count"
    assert manifest["output_paths"]["movie"].endswith("animation.mp4")


def test_compositor_failure_prevents_encoder(tmp_path, monkeypatch):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    output = tmp_path / "animation"
    chimerax = tmp_path / "ChimeraX.exe"
    ffmpeg = tmp_path / "ffmpeg.exe"
    ffprobe = tmp_path / "ffprobe.exe"
    for executable in (chimerax, ffmpeg, ffprobe):
        executable.write_bytes(b"placeholder")
    encoder_called = False

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        structure = output / "frames" / "structure"
        structure.mkdir(parents=True)
        for index in range(3):
            plt.imsave(structure / f"frame_{index:06d}.png", np.zeros((80, 120, 3)))
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("{}", encoding="utf-8")
        (output / "logs" / "chimerax_render_status.json").write_text(
            json.dumps(
                {
                    "status": "success",
                    "frame_count": 3,
                    "reference_unchanged": True,
                    "moving_restored": True,
                }
            ),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess([], 0, "", "")

    def fake_compose(**kwargs):
        raise RuntimeError("synthetic compositor failure")

    def fake_encode(**kwargs):
        nonlocal encoder_called
        encoder_called = True

    monkeypatch.setattr("cryorole.animation.workflow.execute_chimerax", fake_execute)
    monkeypatch.setattr("cryorole.animation.workflow.compose_frame_sequences", fake_compose)
    monkeypatch.setattr("cryorole.animation.workflow.encode_mp4", fake_encode)
    arguments = _arguments(run_dir, path_csv, session, output, render_mode="execute")
    arguments.extend(
        [
            "--chimerax-bin",
            str(chimerax),
            "--ffmpeg-bin",
            str(ffmpeg),
            "--ffprobe-bin",
            str(ffprobe),
        ]
    )
    with pytest.raises(SystemExit):
        main(arguments)
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert manifest["stage_statuses"]["compositor"] == "failed"
    assert manifest["stage_statuses"]["encoder"] == "not_started"
    assert encoder_called is False


def test_structure_validation_failure_prevents_compositor(tmp_path, monkeypatch):
    run_dir, path_csv, session = _fixture_files(tmp_path)
    output = tmp_path / "animation"
    chimerax = tmp_path / "ChimeraX.exe"
    chimerax.write_bytes(b"placeholder")
    compositor_called = False

    def fake_execute(*, chimerax_bin, cxc_script, log_path):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("{}", encoding="utf-8")
        (output / "logs" / "chimerax_render_status.json").write_text(
            json.dumps(
                {
                    "status": "success",
                    "frame_count": 3,
                    "reference_unchanged": True,
                    "moving_restored": True,
                }
            ),
            encoding="utf-8",
        )
        return subprocess.CompletedProcess([], 0, "", "")

    def fake_compose(**kwargs):
        nonlocal compositor_called
        compositor_called = True

    monkeypatch.setattr("cryorole.animation.workflow.execute_chimerax", fake_execute)
    monkeypatch.setattr("cryorole.animation.workflow.compose_frame_sequences", fake_compose)
    arguments = _arguments(run_dir, path_csv, session, output, render_mode="execute")
    arguments.extend(["--chimerax-bin", str(chimerax), "--no-encode"])
    with pytest.raises(SystemExit):
        main(arguments)
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage_statuses"]["structure_render"] == "failed"
    assert manifest["stage_statuses"]["compositor"] == "not_started"
    assert compositor_called is False
