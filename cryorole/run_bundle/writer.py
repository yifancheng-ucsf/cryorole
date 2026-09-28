"""Transactional creation and validation of run bundles."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import stat
import tempfile
from typing import Any
from uuid import uuid4


BUNDLE_STATE_FILENAME = "bundle_state.json"
COMPLETION_MARKER_FILENAME = ".cryorole_bundle_complete"
REQUIRED_RUN_ARTIFACTS = (
    "data/raw_landscape.npz",
    "data/match_table.csv",
    "reports/match_report.json",
    "reports/density_report.json",
    "run_summary.json",
    "run_report.md",
    "run_manifest.json",
)
VALID_STATES = (
    "preparing",
    "matching",
    "computing_ro",
    "computing_sld",
    "writing",
    "validating",
    "completed",
    "failed",
)


class RunBundleWriter:
    """Write a complete bundle beside its destination and atomically publish it."""

    def __init__(self, target_dir: str | Path, *, overwrite: bool = False, run_id: str | None = None):
        self.target_dir = Path(target_dir).expanduser().resolve()
        self.overwrite = bool(overwrite)
        self.run_id = run_id or str(uuid4())
        self.work_dir: Path | None = None
        self._history: list[dict[str, str]] = []
        self._committed = False

    def __enter__(self) -> "RunBundleWriter":
        parent = self.target_dir.parent
        parent.mkdir(parents=True, exist_ok=True)
        if self.target_dir.exists() and not self.target_dir.is_dir():
            raise ValueError(f"Run output path exists and is not a directory: {self.target_dir}")
        if self.target_dir.exists() and not self.overwrite:
            raise FileExistsError(f"Run output directory already exists: {self.target_dir}")
        self.work_dir = Path(
            tempfile.mkdtemp(prefix=f".{self.target_dir.name}.preparing-{self.run_id[:8]}-", dir=parent)
        )
        self.set_state("preparing")
        return self

    def __exit__(self, exc_type, exc, traceback) -> bool:
        if exc is not None and not self._committed:
            self.fail(exc)
        return False

    @property
    def path(self) -> Path:
        if self.work_dir is None:
            raise RuntimeError("RunBundleWriter has not been entered")
        return self.work_dir

    def set_state(self, state: str, *, detail: str | None = None) -> None:
        if state not in VALID_STATES:
            raise ValueError(f"Unsupported run bundle state: {state}")
        timestamp = datetime.now(timezone.utc).isoformat()
        entry = {"state": state, "timestamp": timestamp}
        if detail:
            entry["detail"] = detail
        self._history.append(entry)
        payload = {
            "artifact_type": "cryorole_run_bundle_state",
            "schema_version": "1",
            "run_id": self.run_id,
            "state": state,
            "updated_at": timestamp,
            "history": self._history,
        }
        self._write_json(self.path / BUNDLE_STATE_FILENAME, payload)

    def public_path(self, path: str | Path) -> Path:
        candidate = Path(path).resolve()
        relative = candidate.relative_to(self.path.resolve())
        return self.target_dir / relative

    def publicize(self, value: Any) -> Any:
        """Replace staging paths in nested report data with final bundle paths."""

        source = str(self.path)
        target = str(self.target_dir)
        if isinstance(value, str):
            return value.replace(source, target)
        if isinstance(value, dict):
            return {key: self.publicize(item) for key, item in value.items()}
        if isinstance(value, list):
            return [self.publicize(item) for item in value]
        if isinstance(value, tuple):
            return tuple(self.publicize(item) for item in value)
        return value

    def commit(self) -> Path:
        """Validate, finalize metadata, and publish the complete bundle."""

        self.set_state("validating")
        self._validate_required_artifacts()
        self._rewrite_json_staging_paths()
        self.set_state("completed")
        self._finalize_manifest()
        (self.path / COMPLETION_MARKER_FILENAME).write_text(
            json.dumps({"run_id": self.run_id, "state": "completed"}) + "\n",
            encoding="utf-8",
        )
        self._publish()
        self._committed = True
        return self.target_dir

    def fail(self, error: BaseException) -> Path | None:
        """Preserve an explicitly failed debug bundle when possible."""

        if self.work_dir is None or not self.work_dir.exists():
            return None
        from cryorole.errors import CancelledError

        if isinstance(error, CancelledError):
            # A caller-requested cancel is not a failure to debug: discard staging.
            shutil.rmtree(self.work_dir, ignore_errors=True)
            self.work_dir = None
            return None
        try:
            self.set_state("failed", detail=f"{type(error).__name__}: {error}")
            self._write_json(
                self.path / "failure_report.json",
                {
                    "artifact_type": "cryorole_run_failure",
                    "schema_version": "1",
                    "run_id": self.run_id,
                    "failed_stage": self._history[-2]["state"] if len(self._history) > 1 else "preparing",
                    "error_type": type(error).__name__,
                    "reason": str(error),
                },
            )
            failed = self.target_dir.parent / f".{self.target_dir.name}.failed-{self.run_id[:8]}"
            if failed.exists():
                failed = self.target_dir.parent / f".{self.target_dir.name}.failed-{self.run_id}"
            os.replace(self.path, failed)
            self.work_dir = failed
            return failed
        except Exception:
            return self.work_dir

    def _validate_required_artifacts(self) -> None:
        missing = [relative for relative in REQUIRED_RUN_ARTIFACTS if not (self.path / relative).is_file()]
        if missing:
            raise ValueError(f"Run bundle is missing required artifacts: {missing}")
        npz_path = self.path / "data/raw_landscape.npz"
        if npz_path.stat().st_size <= 0:
            raise ValueError("raw_landscape.npz is empty")
        import numpy as np

        with np.load(npz_path, allow_pickle=False) as payload:
            if "particle_key" not in payload or len(payload["particle_key"]) == 0:
                raise ValueError("raw_landscape.npz has no particle rows")

    def _rewrite_json_staging_paths(self) -> None:
        for path in self.path.rglob("*.json"):
            if path.name in {BUNDLE_STATE_FILENAME, "failure_report.json"}:
                continue
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, UnicodeDecodeError):
                continue
            self._write_json(path, self.publicize(payload))

    def _finalize_manifest(self) -> None:
        summary_path = self.path / "run_summary.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        summary["run_id"] = self.run_id
        summary["bundle_state"] = "completed"
        summary["completion_marker"] = COMPLETION_MARKER_FILENAME
        self._write_json(summary_path, self.publicize(summary))
        manifest_path = self.path / "run_manifest.json"
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        payload["run_id"] = self.run_id
        payload["bundle_transaction"] = {
            "schema_version": "1",
            "state": "completed",
            "completion_marker": COMPLETION_MARKER_FILENAME,
            "state_history": self._history,
        }
        self._write_json(manifest_path, self.publicize(payload))

    def _publish(self) -> None:
        target = self.target_dir
        staging = self.path
        if not target.exists():
            os.replace(staging, target)
            self.work_dir = target
            return
        backup = target.parent / f".{target.name}.previous-{self.run_id[:8]}"
        if backup.exists():
            raise FileExistsError(f"Transactional backup path already exists: {backup}")
        os.replace(target, backup)
        try:
            os.replace(staging, target)
        except Exception:
            os.replace(backup, target)
            raise
        self.work_dir = target
        _safe_remove_tree(backup, expected_parent=target.parent)

    @staticmethod
    def _write_json(path: Path, payload: Any) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(f".{path.name}.tmp-{uuid4().hex[:8]}")
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, path)


def validate_completed_run_bundle(run_dir: str | Path, *, allow_legacy: bool = True) -> dict[str, Any]:
    """Reject transactional bundles that are not demonstrably complete."""

    path = Path(run_dir)
    manifest_path = path / "run_manifest.json"
    state_path = path / BUNDLE_STATE_FILENAME
    marker_path = path / COMPLETION_MARKER_FILENAME
    manifest: dict[str, Any] = {}
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    transaction = manifest.get("bundle_transaction") if isinstance(manifest, dict) else None
    if transaction is None and not state_path.exists():
        if allow_legacy:
            return {"status": "legacy_compat", "run_id": manifest.get("run_id")}
        raise ValueError("Run bundle has no completion state and cannot be consumed safely")
    state = transaction.get("state") if isinstance(transaction, dict) else None
    if state is None and state_path.is_file():
        state = json.loads(state_path.read_text(encoding="utf-8")).get("state")
    if state != "completed" or not marker_path.is_file():
        raise ValueError(f"Run bundle is not completed (state={state!r}); downstream command refused")
    return {"status": "completed", "run_id": manifest.get("run_id")}


def _safe_remove_tree(path: Path, *, expected_parent: Path) -> None:
    resolved = path.resolve()
    parent = expected_parent.resolve()
    if resolved.parent != parent or resolved == parent:
        raise ValueError(f"Refusing to remove transactional backup outside target parent: {resolved}")
    def make_writable_and_retry(function, target, _error_info):
        os.chmod(target, stat.S_IWRITE)
        function(target)

    shutil.rmtree(resolved, onerror=make_writable_and_retry)
