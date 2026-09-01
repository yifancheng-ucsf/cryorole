from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from cryorole.workflow import build_guide_plan, derive_next_actions, inspect_run_status


def _bundle(tmp_path: Path) -> Path:
    run = tmp_path / "run"
    (run / "data").mkdir(parents=True)
    (run / "selections" / "state_1").mkdir(parents=True)
    arrays = {
        "artifact_type": np.asarray("raw_landscape"),
        "schema_version": np.asarray("1"),
        "particle_key": np.asarray(["p1", "p2"]),
        "coordinates_analysis": np.zeros((2, 3)),
        "sld_unfloored": np.ones(2), "sld_raw": np.ones(2),
        "sld_display": np.ones(2), "sld_display_is_outlier": np.zeros(2, bool),
        "sld_was_floored": np.zeros(2, bool), "sld_local_k_mean": np.ones(2),
        "sld_effective_local_k_mean": np.ones(2), "sld_distance_floor": np.ones(2),
        "ref_source_row_id": np.arange(2), "mov_source_row_id": np.arange(2),
    }
    np.savez_compressed(run / "data" / "raw_landscape.npz", **arrays)
    manifest = {
        "run_id": "run-1",
        "bundle_transaction": {"state": "completed"},
        "source_identities": {},
    }
    (run / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    (run / "run_summary.json").write_text(json.dumps({"run_id": "run-1"}), encoding="utf-8")
    (run / ".cryorole_bundle_complete").write_text("{}", encoding="utf-8")
    (run / "selections" / "state_1" / "selection.json").write_text(
        json.dumps({"selection_id": "state_1", "parent_run_id": "run-1"}), encoding="utf-8"
    )
    return run


def test_status_is_artifact_derived_and_next_never_confuses_selection_with_visualization(tmp_path) -> None:
    run = _bundle(tmp_path)
    status = inspect_run_status(run)
    actions = derive_next_actions(status)

    assert status["bundle_status"] == "completed"
    assert status["run_id"] == "run-1"
    assert status["raw_landscape"]["row_count"] == 2
    assert status["selections"][0]["selection_id"] == "state_1"
    assert any("cryorole export" in action["command"] for action in actions)
    assert any(action["category"] == "optional" and "canonicalize" in action["command"] for action in actions)


def test_status_reports_incomplete_transaction_without_crashing(tmp_path) -> None:
    run = tmp_path / "partial"
    run.mkdir()
    (run / "bundle_state.json").write_text(json.dumps({"state": "computing_sld"}), encoding="utf-8")

    status = inspect_run_status(run)

    assert status["bundle_status"] == "incomplete"
    assert status["warnings"]


def test_noninteractive_guide_returns_plan_and_never_waits(tmp_path) -> None:
    run = _bundle(tmp_path)
    plan = build_guide_plan(run_dir=run, non_interactive=True)

    assert plan["mode"] == "resume"
    assert plan["non_interactive"] is True
    assert plan["requires_user_input"] is False
    assert plan["actions"]


def test_new_guide_reports_blocked_preflight_without_execution(tmp_path) -> None:
    plan = build_guide_plan(
        ref=tmp_path / "missing.star", mov=tmp_path / "missing.star",
        output_dir=tmp_path / "never", non_interactive=True,
    )
    assert plan["mode"] == "new"
    assert plan["preflight"]["readiness"] == "BLOCKED"
    assert plan["requires_user_input"] is False
    assert not (tmp_path / "never").exists()


def test_status_reports_selected_derived_landscape_and_completed_export(tmp_path) -> None:
    run = _bundle(tmp_path)
    selected = run / "selections" / "state_1" / "selected_landscape"
    selected.mkdir()
    (selected / "landscape.npz").write_bytes(b"present")
    exported = run / "exports" / "state_1"
    exported.mkdir(parents=True)
    (exported / "export_report.json").write_text("{}", encoding="utf-8")
    (exported / "selected.cs").write_bytes(b"metadata")

    status = inspect_run_status(run)
    actions = derive_next_actions(status)

    assert status["selections"][0]["selected_landscape_present"] is True
    assert status["exports"][0]["report_present"] is True
    assert status["exports"][0]["file_count"] == 2
    assert any(action["category"] == "optional" and "status" in action["command"] for action in actions)


def test_next_does_not_choose_between_multiple_canonical_frames(tmp_path) -> None:
    run = _bundle(tmp_path)
    for canonical_id in ("a", "b"):
        target = run / "canonical" / canonical_id
        target.mkdir(parents=True)
        (target / "canonical_frame.json").write_text("{}", encoding="utf-8")
        (target / "canonical_landscape.npz").write_bytes(b"present")
    status = inspect_run_status(run)
    actions = derive_next_actions(status)
    assert len(status["canonical_frames"]) == 2
    assert not any("--canonical-id" in action["command"] for action in actions)
    assert any("--space raw" in action["command"] for action in actions)
