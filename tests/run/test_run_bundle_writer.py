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


@pytest.mark.parametrize("error_type", [OSError, KeyboardInterrupt])
def test_final_publish_failure_rolls_back_old_bundle(tmp_path, monkeypatch, error_type) -> None:
    target = tmp_path / "run"
    target.mkdir()
    sentinel = target / "old_bundle.txt"
    sentinel.write_text("unchanged", encoding="utf-8")
    real_replace = os.replace

    with pytest.raises(error_type, match="injected final publish failure"):
        with RunBundleWriter(target, overwrite=True, run_id="commit-fail") as writer:
            _write_minimal_required_bundle(writer.path)
            staging = writer.path.resolve()
            resolved_target = target.resolve()

            def failing_replace(source, destination):
                if Path(source).resolve() == staging and Path(destination).resolve() == resolved_target:
                    raise error_type("injected final publish failure")
                return real_replace(source, destination)

            monkeypatch.setattr("cryorole.run_bundle.writer.os.replace", failing_replace)
            writer.commit()

    assert sentinel.read_text(encoding="utf-8") == "unchanged"
    assert not (tmp_path / ".run.publish-lock").exists()


def test_two_prepared_runs_cannot_overwrite_without_permission(tmp_path) -> None:
    target = tmp_path / "run"
    with RunBundleWriter(target, run_id="winner") as first:
        _write_minimal_required_bundle(first.path)
        with pytest.raises(FileExistsError, match="already exists"):
            with RunBundleWriter(target, run_id="loser") as second:
                _write_minimal_required_bundle(second.path)
                first.commit()
                second.commit()

    assert validate_completed_run_bundle(target, allow_legacy=False)["run_id"] == "winner"
    assert (tmp_path / ".run.failed-loser").is_dir()
    assert not (tmp_path / ".run.publish-lock").exists()


def test_commit_is_idempotent(tmp_path) -> None:
    with RunBundleWriter(tmp_path / "run", run_id="one-run") as writer:
        _write_minimal_required_bundle(writer.path)
        target = writer.commit()
        state = (target / "bundle_state.json").read_bytes()
        manifest = (target / "run_manifest.json").read_bytes()
        assert writer.commit() == target
        assert (target / "bundle_state.json").read_bytes() == state
        assert (target / "run_manifest.json").read_bytes() == manifest


def test_publish_lock_protects_gap_during_overwrite(tmp_path, monkeypatch) -> None:
    target = tmp_path / "run"
    target.mkdir()
    (target / "old.txt").write_text("old", encoding="utf-8")
    real_replace = os.replace
    checked = []
    with RunBundleWriter(target, overwrite=True, run_id="winner") as first:
        with RunBundleWriter(target, overwrite=True, run_id="contender") as second:
            _write_minimal_required_bundle(first.path)
            _write_minimal_required_bundle(second.path)

            def concurrent_replace(source, destination):
                result = real_replace(source, destination)
                if Path(source) == target and ".previous-" in Path(destination).name:
                    assert not target.exists()
                    with pytest.raises(FileExistsError, match="Another run is publishing"):
                        second.commit()
                    checked.append(True)
                return result

            monkeypatch.setattr("cryorole.run_bundle.writer.os.replace", concurrent_replace)
            first.commit()
    assert checked == [True]
    assert validate_completed_run_bundle(target, allow_legacy=False)["run_id"] == "winner"
    assert not (tmp_path / ".run.publish-lock").exists()


def test_backup_cleanup_failure_preserves_published_bundle(tmp_path, monkeypatch) -> None:
    target = tmp_path / "run"
    target.mkdir()
    (target / "old.txt").write_text("old", encoding="utf-8")

    def failed_cleanup(*args, **kwargs):
        raise OSError("injected backup cleanup failure")

    monkeypatch.setattr("cryorole.run_bundle.writer._safe_remove_tree", failed_cleanup)
    with pytest.raises(OSError, match="backup cleanup failure"):
        with RunBundleWriter(target, overwrite=True, run_id="new-run") as writer:
            _write_minimal_required_bundle(writer.path)
            writer.commit()

    assert validate_completed_run_bundle(target, allow_legacy=False)["run_id"] == "new-run"
    assert (tmp_path / ".run.previous-new-run/old.txt").read_text(encoding="utf-8") == "old"
    assert not (target / "failure_report.json").exists()
    assert writer.commit() == target
    assert not (tmp_path / ".run.publish-lock").exists()
