"""Read-only run status derived from immutable and downstream artifacts."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np

from cryorole.provenance import verify_source_identity
from cryorole.run_bundle.writer import COMPLETION_MARKER_FILENAME, REQUIRED_RUN_ARTIFACTS


STATUS_SCHEMA_VERSION = "1.0"


def inspect_run_status(run_dir: str | Path) -> dict[str, Any]:
    """Inspect actual files; no mutable workflow-state file is trusted alone."""

    root = Path(run_dir).expanduser().resolve()
    warnings: list[str] = []
    if not root.is_dir():
        return {
            "artifact_type": "cryorole_workflow_status",
            "schema_version": STATUS_SCHEMA_VERSION,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "run_dir": str(root),
            "bundle_status": "missing",
            "run_id": None,
            "warnings": ["Run directory does not exist"],
            "missing_required_artifacts": list(REQUIRED_RUN_ARTIFACTS),
            "raw_landscape": None,
            "visualizations": [], "canonical_frames": [], "selections": [], "exports": [],
            "source_metadata": {},
        }

    manifest = _read_json(root / "run_manifest.json", warnings)
    summary = _read_json(root / "run_summary.json", warnings)
    state_payload = _read_json(root / "bundle_state.json", warnings)
    run_id = manifest.get("run_id") or summary.get("run_id") or state_payload.get("run_id")
    transaction = manifest.get("bundle_transaction") if isinstance(manifest, dict) else None
    transaction_state = transaction.get("state") if isinstance(transaction, dict) else None
    state = transaction_state or state_payload.get("state")
    marker = (root / COMPLETION_MARKER_FILENAME).is_file()
    if state == "failed":
        bundle_status = "failed"
    elif state == "completed" and marker:
        bundle_status = "completed"
    elif transaction is None and not state_payload:
        bundle_status = "legacy"
        warnings.append("Legacy bundle has no transactional completion metadata")
    else:
        bundle_status = "incomplete"
        warnings.append(f"Transactional bundle is not complete (state={state!r}, marker={marker})")

    missing = [relative for relative in REQUIRED_RUN_ARTIFACTS if not (root / relative).is_file()]
    if missing:
        warnings.append(f"Missing required artifacts: {missing}")
    raw = _inspect_raw_npz(root / "data" / "raw_landscape.npz", warnings)
    source_identities = manifest.get("source_identities") or summary.get("source_identities") or {}
    source_metadata = {
        domain: _source_status(record)
        for domain, record in source_identities.items()
        if isinstance(record, dict)
    }
    return {
        "artifact_type": "cryorole_workflow_status",
        "schema_version": STATUS_SCHEMA_VERSION,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "run_dir": str(root),
        "bundle_status": bundle_status,
        "transaction_state": state,
        "completion_marker": marker,
        "run_id": run_id,
        "provenance_present": bool(source_identities),
        "source_metadata": source_metadata,
        "raw_landscape": raw,
        "visualizations": _scan_dirs(root / "visualizations"),
        "canonical_frames": _canonical_status(root),
        "selections": _selection_status(root, run_id=run_id, warnings=warnings),
        "exports": _export_status(root),
        "missing_required_artifacts": missing,
        "artifact_integrity": "ok" if not missing and raw and raw.get("valid") else "warnings",
        "warnings": warnings,
    }


def _read_json(path: Path, warnings: list[str]) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        warnings.append(f"Cannot read {path.name}: {exc}")
        return {}
    return payload if isinstance(payload, dict) else {}


def _inspect_raw_npz(path: Path, warnings: list[str]) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        with np.load(path, allow_pickle=False) as payload:
            count = len(payload["particle_key"])
            schema = str(payload["schema_version"].item())
            coordinates_shape = list(np.asarray(payload["coordinates_analysis"]).shape)
        return {
            "path": str(path), "valid": count > 0 and coordinates_shape == [count, 3],
            "row_count": count, "schema_version": schema,
            "coordinates_shape": coordinates_shape, "size_bytes": path.stat().st_size,
        }
    except Exception as exc:
        warnings.append(f"Raw landscape validation failed: {exc}")
        return {"path": str(path), "valid": False, "error": str(exc)}


def _source_status(record: dict[str, object]) -> dict[str, object]:
    try:
        verification = verify_source_identity(record)
        return {
            "path": verification.resolved_path,
            "available": True,
            "hash_status": "verified",
        }
    except FileNotFoundError as exc:
        return {"path": record.get("resolved_path"), "available": False, "hash_status": "missing", "detail": str(exc)}
    except ValueError as exc:
        return {"path": record.get("resolved_path"), "available": True, "hash_status": "mismatch_or_unverified", "detail": str(exc)}


def _scan_dirs(path: Path) -> list[dict[str, object]]:
    if not path.is_dir():
        return []
    return [
        {"id": str(child.relative_to(path)).replace("\\", "/"), "path": str(child)}
        for child in sorted(path.rglob("*")) if child.is_dir() and any(item.is_file() for item in child.iterdir())
    ]


def _canonical_status(root: Path) -> list[dict[str, object]]:
    canonical_root = root / "canonical"
    if not canonical_root.is_dir():
        return []
    results = []
    for child in sorted(canonical_root.iterdir()):
        if child.is_dir():
            results.append({
                "canonical_id": child.name,
                "path": str(child),
                "frame_present": (child / "canonical_frame.json").is_file(),
                "landscape_present": (child / "canonical_landscape.npz").is_file(),
            })
    return results

def _selection_status(root: Path, *, run_id: object, warnings: list[str]) -> list[dict[str, object]]:
    selection_root = root / "selections"
    if not selection_root.is_dir():
        return []
    results = []
    for child in sorted(selection_root.iterdir()):
        payload = _read_json(child / "selection.json", warnings)
        parent_run_id = payload.get("parent_run_id")
        if run_id and parent_run_id and parent_run_id != run_id:
            warnings.append(f"Selection {child.name} belongs to different run_id {parent_run_id}")
        results.append({
            "selection_id": payload.get("selection_id", child.name),
            "path": str(child), "parent_run_id": parent_run_id,
            "selected_count": payload.get("selected_count"),
            "selected_landscape_present": (child / "selected_landscape" / "landscape.npz").is_file(),
        })
    return results


def _export_status(root: Path) -> list[dict[str, object]]:
    export_root = root / "exports"
    if not export_root.is_dir():
        return []
    results = []
    for child in sorted(export_root.iterdir()):
        if child.is_dir():
            report = child / "export_report.json"
            results.append({
                "selection_id": child.name, "path": str(child),
                "report_present": report.is_file(),
                "file_count": sum(1 for item in child.rglob("*") if item.is_file()),
            })
    return results
