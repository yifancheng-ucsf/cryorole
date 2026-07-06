"""Run manifest writing for cryoROLE provenance."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from cryorole import __version__
from cryorole.export.serialization import to_json_safe
from cryorole.models.manifest_report import ManifestReport
from cryorole.models.policies import RunManifestPolicy


def write_run_manifest(
    *,
    policy: RunManifestPolicy,
    input_paths: Sequence[str | Path] | None = None,
    source_types: Mapping[str, str] | None = None,
    row_counts: Mapping[str, int] | None = None,
    import_reports: Any = None,
    active_policies: Mapping[str, Any] | None = None,
    match_report: Any = None,
    density_report: Any = None,
    canonicalization_report: Any = None,
    selection: Any = None,
    export_report: Any = None,
    additional_reports: Mapping[str, Any] | None = None,
    additional_results: Mapping[str, Any] | None = None,
    output_artifacts: Mapping[str, Any] | None = None,
) -> ManifestReport:
    """Write a non-destructive JSON manifest for a full or partial workflow."""

    output_path = Path(policy.output_path)
    if output_path.is_dir():
        raise IsADirectoryError(f"Manifest output path is a directory: {output_path}")
    _validate_manifest_output_path(output_path)
    if policy.compute_file_hashes:
        _new_hash(policy.hash_algorithm)
    if output_path.exists() and not policy.overwrite:
        raise FileExistsError(f"Manifest output path already exists: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now(timezone.utc).isoformat()
    manifest = build_run_manifest_payload(
        policy=policy,
        timestamp=timestamp,
        input_paths=input_paths,
        source_types=source_types,
        row_counts=row_counts,
        import_reports=import_reports,
        active_policies=active_policies,
        match_report=match_report,
        density_report=density_report,
        canonicalization_report=canonicalization_report,
        selection=selection,
        export_report=export_report,
        additional_reports=additional_reports,
        additional_results=additional_results,
        output_artifacts=output_artifacts,
    )
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(to_json_safe(manifest), handle, indent=2, sort_keys=True)
        handle.write("\n")

    return ManifestReport(
        output_path=str(output_path),
        size_bytes=output_path.stat().st_size,
        manifest_type=policy.manifest_type,
        schema_version=policy.schema_version,
        timestamp=timestamp,
        manifest_policy=to_json_safe(policy),
    )


def build_run_manifest_payload(
    *,
    policy: RunManifestPolicy,
    timestamp: str,
    input_paths: Sequence[str | Path] | None = None,
    source_types: Mapping[str, str] | None = None,
    row_counts: Mapping[str, int] | None = None,
    import_reports: Any = None,
    active_policies: Mapping[str, Any] | None = None,
    match_report: Any = None,
    density_report: Any = None,
    canonicalization_report: Any = None,
    selection: Any = None,
    export_report: Any = None,
    additional_reports: Mapping[str, Any] | None = None,
    additional_results: Mapping[str, Any] | None = None,
    output_artifacts: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a JSON-safe manifest payload without writing it."""

    input_provenance = {
        "inputs": _input_records(input_paths or (), policy=policy),
        "source_types": to_json_safe(source_types or {}),
        "row_counts": to_json_safe(row_counts or {}),
        "import_reports": to_json_safe(import_reports),
    }
    reports = {
        "match_report": to_json_safe(match_report),
        "density_report": to_json_safe(density_report),
        "canonicalization_report": to_json_safe(canonicalization_report),
    }
    if additional_reports:
        reports.update(to_json_safe(additional_reports))

    results = {
        "selection_summary": _selection_summary(selection, policy=policy),
    }
    if additional_results:
        results.update(to_json_safe(additional_results))

    artifacts = dict(to_json_safe(output_artifacts or {}))
    if export_report is not None:
        safe_export_report = to_json_safe(export_report)
        reports["export_report"] = safe_export_report
        artifacts["export_report_output_paths"] = safe_export_report.get("output_paths", {})
        artifacts["export_report_row_counts"] = safe_export_report.get("row_counts", {})
    else:
        reports["export_report"] = None

    return {
        "manifest_type": policy.manifest_type,
        "schema_version": policy.schema_version,
        "timestamp": timestamp,
        "cryorole_version": __version__,
        "workflow_name": policy.workflow_name,
        "command": policy.command,
        "input_provenance": input_provenance,
        "active_policies": to_json_safe(active_policies or {}),
        "reports": reports,
        "results": results,
        "output_artifacts": artifacts,
        "manifest_policy": to_json_safe(policy),
    }


def _input_records(
    input_paths: Sequence[str | Path],
    *,
    policy: RunManifestPolicy,
) -> list[dict[str, Any]]:
    records = []
    for input_path in input_paths:
        path = Path(input_path)
        record: dict[str, Any] = {"path": str(path)}
        if path.exists():
            record["exists"] = True
            record["size_bytes"] = path.stat().st_size
            if policy.compute_file_hashes:
                record["hash_algorithm"] = policy.hash_algorithm
                record["hash"] = _file_hash(path, policy.hash_algorithm)
        else:
            record["exists"] = False
            if policy.compute_file_hashes:
                record["hash_algorithm"] = policy.hash_algorithm
                record["hash"] = None
        records.append(record)
    return records


def _validate_manifest_output_path(output_path: Path) -> None:
    if output_path.suffix.lower() != ".json":
        raise ValueError(f"Manifest output path must have a .json suffix: {output_path}")


def _new_hash(algorithm: str):
    try:
        return hashlib.new(algorithm)
    except ValueError as exc:
        raise ValueError(f"Unsupported manifest hash_algorithm: {algorithm}") from exc


def _file_hash(path: Path, algorithm: str) -> str:
    hasher = _new_hash(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def _selection_summary(selection: Any, *, policy: RunManifestPolicy) -> Any:
    if selection is None:
        return None
    safe_selection = to_json_safe(selection)
    summary = {
        "selection_id": safe_selection.get("selection_id"),
        "parent_landscape_id": safe_selection.get("parent_landscape_id"),
        "selection_mode": safe_selection.get("selection_mode"),
        "selection_basis": safe_selection.get("selection_basis"),
        "metric": safe_selection.get("metric"),
        "selected_count": safe_selection.get("selected_count"),
        "total_count": safe_selection.get("total_count"),
        "selected_particle_keys_included": policy.include_selected_particle_keys,
        "active_policy": safe_selection.get("active_policy"),
    }
    if policy.include_selected_particle_keys:
        summary["selected_particle_keys"] = safe_selection.get("selected_particle_keys")
    return summary
