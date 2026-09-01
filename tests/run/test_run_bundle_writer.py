from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from cryorole.run_bundle import RunBundleWriter, validate_completed_run_bundle


def _write_minimal_required_bundle(path: Path) -> None:
    (path / "data").mkdir(parents=True, exist_ok=True)
    (path / "reports").mkdir(parents=True, exist_ok=True)
    np.savez(
        path / "data/raw_landscape.npz",
        particle_key=np.asarray(["p1"]),
    )
    (path / "data/match_table.csv").write_text("particle_key\np1\n", encoding="utf-8")
    for name in ("match_report.json", "density_report.json"):
        (path / "reports" / name).write_text("{}\n", encoding="utf-8")
    (path / "run_summary.json").write_text("{}\n", encoding="utf-8")
    (path / "run_report.md").write_text("# Run\n", encoding="utf-8")
    (path / "run_manifest.json").write_text("{}\n", encoding="utf-8")


def test_transaction_commit_replaces_whole_old_bundle_without_stale_children(tmp_path) -> None:
    target = tmp_path / "run"
    (target / "selections/old").mkdir(parents=True)
    (target / "selections/old/selection.json").write_text("{}", encoding="utf-8")
    with RunBundleWriter(target, overwrite=True, run_id="run-new") as writer:
        _write_minimal_required_bundle(writer.path)
        writer.commit()

    assert not (target / "selections/old").exists()
    status = validate_completed_run_bundle(target, allow_legacy=False)
    assert status == {"status": "completed", "run_id": "run-new"}


def test_failed_validation_preserves_existing_bundle_and_keeps_failed_debug(tmp_path) -> None:
    target = tmp_path / "run"
    target.mkdir()
    sentinel = target / "old_bundle.txt"
    sentinel.write_text("unchanged", encoding="utf-8")

    with pytest.raises(ValueError, match="missing required artifacts"):
        with RunBundleWriter(target, overwrite=True, run_id="bad-run") as writer:
            (writer.path / "run_summary.json").write_text("{}", encoding="utf-8")
            writer.commit()

    assert sentinel.read_text(encoding="utf-8") == "unchanged"
    failed = tmp_path / ".run.failed-bad-run"
    assert failed.is_dir()
    report = json.loads((failed / "failure_report.json").read_text(encoding="utf-8"))
    assert report["error_type"] == "ValueError"


def test_incomplete_transactional_bundle_is_rejected(tmp_path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "bundle_state.json").write_text(
        json.dumps({"state": "computing_sld", "run_id": "partial"}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="not completed"):
        validate_completed_run_bundle(run_dir)


def test_final_publish_failure_rolls_back_old_bundle(tmp_path, monkeypatch) -> None:
    target = tmp_path / "run"
    target.mkdir()
    sentinel = target / "old_bundle.txt"
    sentinel.write_text("unchanged", encoding="utf-8")
    real_replace = os.replace

    with pytest.raises(OSError, match="injected final publish failure"):
        with RunBundleWriter(target, overwrite=True, run_id="commit-fail") as writer:
            _write_minimal_required_bundle(writer.path)
            staging = writer.path.resolve()
            resolved_target = target.resolve()

            def failing_replace(source, destination):
                if Path(source).resolve() == staging and Path(destination).resolve() == resolved_target:
                    raise OSError("injected final publish failure")
                return real_replace(source, destination)

            monkeypatch.setattr("cryorole.run_bundle.writer.os.replace", failing_replace)
            writer.commit()

    assert sentinel.read_text(encoding="utf-8") == "unchanged"
