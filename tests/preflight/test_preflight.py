from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from cryorole.cli.main import build_parser, preflight_command, run_command
from cryorole.preflight import PreflightRequest, run_preflight


def _write_cs(path: Path, uids: list[int], poses: np.ndarray) -> None:
    dtype = np.dtype([("uid", "<u8"), ("alignments3D/pose", "<f8", (3,))])
    values = np.zeros(len(uids), dtype=dtype)
    values["uid"] = uids
    values["alignments3D/pose"] = poses
    with path.open("wb") as handle:
        np.save(handle, values)


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    ref = tmp_path / "ref.cs"
    mov = tmp_path / "mov.cs"
    poses = np.arange(18, dtype=float).reshape(6, 3) / 100.0
    _write_cs(ref, [1, 2, 3, 4, 5, 6], poses)
    _write_cs(mov, [3, 1, 6, 2, 5, 4], poses + 0.02)
    return ref, mov


def _write_star(path: Path, names: list[str], *, include_psi: bool = True) -> None:
    columns = ["_rlnImageName #1", "_rlnAngleRot #2", "_rlnAngleTilt #3"]
    if include_psi:
        columns.append("_rlnAnglePsi #4")
    rows = [f"{name} {index} {index + 1}" + (f" {index + 2}" if include_psi else "") for index, name in enumerate(names)]
    path.write_text("\n".join(["data_particles", "", "loop_", *columns, *rows, ""]), encoding="utf-8")


def test_preflight_ready_is_structured_deterministic_and_side_effect_free(tmp_path) -> None:
    ref, mov = _inputs(tmp_path)
    output = tmp_path / "never-created"
    request = PreflightRequest(ref=ref, mov=mov, output_dir=output)

    first = run_preflight(request)
    second = run_preflight(request)

    assert first.report["readiness"] == "READY_WITH_WARNINGS"
    assert first.report["matching"]["matched_count"] == 6
    assert first.report["matching"]["reordered"] is True
    assert first.report["source_identities"]["ref"]["sha256"]
    stable_resource_fields = {
        key: value
        for key, value in first.report["resource_estimate"].items()
        if key != "available_disk_bytes"
    }
    assert stable_resource_fields == {
        key: value
        for key, value in second.report["resource_estimate"].items()
        if key != "available_disk_bytes"
    }
    assert first.report["resolved_run_command"].startswith("cryorole run ")
    assert not output.exists()


def test_preflight_blocks_low_overlap_and_invalid_pose(tmp_path) -> None:
    ref, _ = _inputs(tmp_path)
    mov = tmp_path / "mov.cs"
    poses = np.zeros((6, 3), dtype=float)
    poses[0, 0] = np.nan
    _write_cs(mov, [1, 20, 21, 22, 23, 24], poses)

    result = run_preflight(PreflightRequest(ref=ref, mov=mov))

    assert result.report["readiness"] == "BLOCKED"
    assert result.exit_code == 2
    assert result.report["errors"]
    assert result.report["recommended_next_command"]


def test_preflight_cli_json_and_run_dry_run_have_parity_without_bundle(tmp_path) -> None:
    ref, mov = _inputs(tmp_path)
    report_path = tmp_path / "preflight.json"
    output = tmp_path / "run"
    parser = build_parser()
    preflight_args = parser.parse_args(
        ["preflight", "--ref", str(ref), "--mov", str(mov), "--json", str(report_path)]
    )
    dry_args = parser.parse_args(
        [
            "run", "--ref", str(ref), "--mov", str(mov), "--output-dir", str(output),
            "--dry-run", "--json", "-",
        ]
    )

    assert preflight_command(preflight_args) == 1
    assert run_command(dry_args) == 1
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["matching"]["matched_count"] == 6
    assert payload["readiness"] == "READY_WITH_WARNINGS"
    assert not output.exists()


def test_preflight_identity_detects_input_changed_after_report(tmp_path) -> None:
    ref, mov = _inputs(tmp_path)
    result = run_preflight(PreflightRequest(ref=ref, mov=mov))
    with ref.open("ab") as handle:
        handle.write(b"changed")

    assert result.inputs_unchanged() is False


def test_preflight_identity_detects_same_size_same_mtime_content_change(tmp_path) -> None:
    ref, mov = _inputs(tmp_path)
    result = run_preflight(PreflightRequest(ref=ref, mov=mov))
    original_stat = ref.stat()
    payload = bytearray(ref.read_bytes())
    payload[-1] ^= 1
    ref.write_bytes(payload)
    os.utime(ref, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))

    assert ref.stat().st_size == original_stat.st_size
    assert ref.stat().st_mtime_ns == original_stat.st_mtime_ns
    assert result.inputs_unchanged() is False


def test_explicit_mapping_file_has_preflight_and_run_dry_run_policy_parity(
    tmp_path,
    capsys,
) -> None:
    ref, mov = tmp_path / "ref.star", tmp_path / "mov.star"
    _write_star(ref, ["ref-a", "ref-b", "ref-c"])
    _write_star(mov, ["mov-a", "mov-b", "mov-c"])
    mapping = tmp_path / "mapping.csv"
    mapping.write_text(
        "source_row_id,particle_key\n0,p3\n1,p1\n2,p2\n",
        encoding="utf-8",
    )
    parser = build_parser()
    common = [
        "--ref", str(ref), "--mov", str(mov),
        "--identity-mode", "explicit_mapping_file",
        "--mapping-file", str(mapping),
    ]

    preflight_args = parser.parse_args(["preflight", *common, "--json", "-"])
    dry_run_args = parser.parse_args(["run", *common, "--dry-run", "--json", "-"])
    assert preflight_command(preflight_args) == 0
    preflight_payload = json.loads(capsys.readouterr().out)
    assert run_command(dry_run_args) == 0
    run_payload = json.loads(capsys.readouterr().out)

    assert preflight_payload["identity_policy"] == run_payload["identity_policy"]
    assert preflight_payload["identity_policy"]["ref"]["identity_mode"] == "explicit_mapping_file"
    assert preflight_payload["identity_policy"]["ref"]["mapping_file"] == str(mapping)
    assert preflight_payload["matching"]["matched_count"] == 3


def test_preflight_star_schema_and_convention_resolution(tmp_path) -> None:
    ref, mov = tmp_path / "ref.star", tmp_path / "mov.star"
    _write_star(ref, ["1@a.mrcs", "2@a.mrcs"])
    _write_star(mov, ["2@a.mrcs", "1@a.mrcs"])

    result = run_preflight(PreflightRequest(ref=ref, mov=mov))

    assert result.report["readiness"] == "READY_WITH_WARNINGS"
    assert result.report["convention_resolution"]["ref"]["source_euler_sequence"] == "ZYZ"
    assert result.report["matching"]["identity_key"] == "_rlnImageName"


def test_preflight_reports_missing_pose_field_and_wrong_cs_vector_shape(tmp_path) -> None:
    ref, mov = tmp_path / "ref.star", tmp_path / "mov.star"
    _write_star(ref, ["1@a.mrcs"], include_psi=False)
    _write_star(mov, ["1@a.mrcs"])
    assert run_preflight(PreflightRequest(ref=ref, mov=mov)).report["readiness"] == "BLOCKED"

    bad = tmp_path / "bad.cs"
    dtype = np.dtype([("uid", "<u8"), ("alignments3D/pose", "<f8", (2,))])
    values = np.zeros(1, dtype=dtype)
    with bad.open("wb") as handle:
        np.save(handle, values)
    good = tmp_path / "good.cs"
    _write_cs(good, [1], np.zeros((1, 3)))
    result = run_preflight(PreflightRequest(ref=good, mov=bad))
    assert result.report["readiness"] == "BLOCKED"
    assert "shape" in " ".join(result.report["errors"])


def test_preflight_duplicate_and_zero_overlap_diagnostics(tmp_path) -> None:
    ref = tmp_path / "ref.cs"
    mov = tmp_path / "mov.cs"
    _write_cs(ref, [1, 1, 2], np.zeros((3, 3)))
    _write_cs(mov, [1, 2, 3], np.zeros((3, 3)))
    duplicate = run_preflight(PreflightRequest(ref=ref, mov=mov))
    assert duplicate.report["readiness"] == "BLOCKED"
    assert duplicate.report["identity_reports"]["ref"]["duplicate_count"] == 2
    assert duplicate.report["identity_reports"]["ref"]["collision_examples"]

    _write_cs(ref, [10, 11], np.zeros((2, 3)))
    _write_cs(mov, [20, 21], np.zeros((2, 3)))
    zero = run_preflight(PreflightRequest(ref=ref, mov=mov))
    assert zero.report["matching"]["matched_count"] == 0
    assert zero.report["readiness"] == "BLOCKED"


def test_preflight_row_aligned_equal_and_unequal_counts(tmp_path) -> None:
    ref, mov = tmp_path / "ref.cs", tmp_path / "mov.cs"
    _write_cs(ref, [1, 2], np.zeros((2, 3)))
    _write_cs(mov, [9, 8], np.zeros((2, 3)))
    ready = run_preflight(PreflightRequest(ref=ref, mov=mov, row_aligned=True))
    assert ready.report["readiness"] == "READY"
    assert ready.report["matching"]["row_aligned"] is True
    _write_cs(mov, [9], np.zeros((1, 3)))
    blocked = run_preflight(PreflightRequest(ref=ref, mov=mov, row_aligned=True))
    assert blocked.report["readiness"] == "BLOCKED"


def test_preflight_relative_paths_record_original_and_resolved_paths(tmp_path, monkeypatch) -> None:
    ref, mov = _inputs(tmp_path)
    monkeypatch.chdir(tmp_path)
    result = run_preflight(PreflightRequest(ref=Path("ref.cs"), mov=Path("mov.cs")))
    assert result.report["source_identities"]["ref"]["original_path"] == "ref.cs"
    assert result.report["source_identities"]["ref"]["resolved_path"] == str(ref.resolve())


def test_preflight_low_overlap_requires_explicit_policy_and_reports_dropped_rows(tmp_path) -> None:
    ref, mov = tmp_path / "ref.cs", tmp_path / "mov.cs"
    _write_cs(ref, [1, 2, 3, 4], np.zeros((4, 3)))
    _write_cs(mov, [1, 20, 30, 40], np.zeros((4, 3)))

    blocked = run_preflight(PreflightRequest(ref=ref, mov=mov))
    allowed = run_preflight(PreflightRequest(ref=ref, mov=mov, allow_low_overlap=True))

    assert blocked.report["readiness"] == "BLOCKED"
    assert allowed.report["matching"]["matched_count"] == 1
    assert allowed.report["matching"]["ref_only_count"] == 3
    assert allowed.report["matching"]["mov_only_count"] == 3
    assert allowed.report["matching"]["low_overlap_allowed"] is True
    assert allowed.report["readiness"] == "READY_WITH_WARNINGS"
