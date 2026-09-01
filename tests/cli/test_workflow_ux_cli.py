from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from cryorole.cli.main import build_parser, guide_command, next_command, status_command


def _bundle(tmp_path: Path) -> Path:
    run = tmp_path / "run"
    (run / "data").mkdir(parents=True)
    np.savez_compressed(
        run / "data" / "raw_landscape.npz",
        artifact_type=np.asarray("raw_landscape"), schema_version=np.asarray("1"),
        particle_key=np.asarray(["p1"]), coordinates_analysis=np.zeros((1, 3)),
        sld_unfloored=np.ones(1), sld_raw=np.ones(1), sld_display=np.ones(1),
        sld_display_is_outlier=np.zeros(1, bool), sld_was_floored=np.zeros(1, bool),
        sld_local_k_mean=np.ones(1), sld_effective_local_k_mean=np.ones(1),
        sld_distance_floor=np.ones(1), ref_source_row_id=np.zeros(1, int),
        mov_source_row_id=np.zeros(1, int),
    )
    (run / "run_manifest.json").write_text(
        json.dumps({"run_id": "ux", "bundle_transaction": {"state": "completed"}}), encoding="utf-8"
    )
    (run / "run_summary.json").write_text(json.dumps({"run_id": "ux"}), encoding="utf-8")
    (run / ".cryorole_bundle_complete").write_text("{}", encoding="utf-8")
    return run


def test_public_parser_exposes_workflow_ux_commands() -> None:
    parser = build_parser()
    help_text = parser.format_help()
    for command in ("preflight", "status", "next", "guide", "explore"):
        assert command in help_text
    dry = parser.parse_args(["run", "--ref", "a.cs", "--mov", "b.cs", "--dry-run"])
    assert dry.dry_run is True


def test_workflow_ux_subcommand_help_renders_without_format_errors(capsys) -> None:
    parser = build_parser()
    for command in ("preflight", "run", "status", "next", "guide", "explore"):
        try:
            parser.parse_args([command, "--help"])
        except SystemExit as exc:
            assert exc.code == 0
        else:
            raise AssertionError(f"{command} --help did not exit normally")
    output = capsys.readouterr().out
    assert "50%" in output
    assert "--max-display-points" in output


def test_status_and_next_support_json_files(tmp_path) -> None:
    run = _bundle(tmp_path)
    parser = build_parser()
    status_path = tmp_path / "status.json"
    next_path = tmp_path / "next.json"

    assert status_command(parser.parse_args(["status", "--run-dir", str(run), "--json", str(status_path)])) == 0
    assert next_command(parser.parse_args(["next", "--run-dir", str(run), "--json", str(next_path)])) == 0

    assert json.loads(status_path.read_text(encoding="utf-8"))["bundle_status"] == "completed"
    actions = json.loads(next_path.read_text(encoding="utf-8"))["actions"]
    assert any(action["command"].startswith("cryorole explore") for action in actions)


def test_guide_noninteractive_does_not_prompt_or_mutate_bundle(tmp_path, monkeypatch) -> None:
    run = _bundle(tmp_path)
    before = sorted(str(path.relative_to(run)) for path in run.rglob("*") if path.is_file())
    monkeypatch.setattr("builtins.input", lambda _prompt: (_ for _ in ()).throw(AssertionError("prompted")))
    args = build_parser().parse_args(["guide", "--run-dir", str(run), "--non-interactive"])

    assert guide_command(args) == 0
    after = sorted(str(path.relative_to(run)) for path in run.rglob("*") if path.is_file())
    assert after == before
