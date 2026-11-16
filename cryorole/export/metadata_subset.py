"""Direct source metadata subset export for selected cryoROLE particles."""

from __future__ import annotations

from dataclasses import asdict, replace
import json
import os
import shutil
import tempfile
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from cryorole.export.cs_group import check_cs_group, resolve_group, write_cs_group_subset
from cryorole.export.cs_subset import write_cryosparc_cs_subset
from cryorole.export.serialization import to_json_safe
from cryorole.export.star_subset import write_relion_star_subset
from cryorole.io.writers.landscape_store import (
    LANDSCAPE_NPZ_SCHEMA_VERSION,
    resolve_landscape_path,
)
from cryorole.models.policies import SelectionMetadataExportPolicy
from cryorole.models.selection import Selection
from cryorole.provenance import SourceIdentityGuard, stream_sha256, verify_source_identity
from cryorole.provenance.output_paths import protect_input_directory
from cryorole.run_bundle import validate_completed_run_bundle

RELION_POSE_HEADERS = ("_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi")

DOMAIN_ROW_ID_COLUMNS = {
    "ref": "ref_source_row_id",
    "mov": "mov_source_row_id",
}


def export_selection_metadata_subset(
    selection: Selection,
    *,
    policy: SelectionMetadataExportPolicy,
    selection_path: str | Path | None = None,
) -> dict[str, Any]:
    """Export source STAR/CS metadata subsets for a Selection.

    This is a bridge-layer operation. It never recomputes selections, never
    reinterprets source orientations, and uses source row IDs as the export
    authority.
    """

    _validate_policy(policy)
    run_dir = Path(policy.run_dir).resolve() if policy.run_dir is not None else None
    if run_dir is not None:
        validate_completed_run_bundle(run_dir)
    selection_path_obj = Path(selection_path).resolve() if selection_path else None
    selection_dir = selection_path_obj.parent if selection_path_obj else None
    output_dir = _resolve_output_dir(
        policy=policy,
        run_dir=run_dir,
        selection=selection,
    )

    domains = _domains(policy.domain)
    source_info, source_verification, run_id = (
        _resolve_and_verify_sources(run_dir, policy=policy, domains=domains)
        if run_dir is not None else ({}, {}, None)
    )
    source_guard = SourceIdentityGuard({domain: {'resolved_path': source_info[domain]['path'],
                                               'sha256': source_info[domain]['sha256'],
                                               'input_path_absolute': source_info[domain]['export_input_path_absolute']}
                                        for domain in domains if domain in source_verification})
    if run_id is not None and selection.parent_run_id not in (None, run_id):
        raise ValueError(
            f"Selection parent_run_id {selection.parent_run_id!r} does not match run_id {run_id!r}"
        )
    selected_rows = _load_selected_rows(
        selection=selection, selection_dir=selection_dir, run_dir=run_dir,
        domains=domains, source_info=source_info,
    )
    # Resolve and check every CryoSPARC result group before anything is written, so a group
    # problem stops the export with nothing staged or published.
    formats = {
        domain: _resolve_format_for_domain(requested=policy.format, domain=domain, source_info=source_info)
        for domain in domains
    }
    group_plans = {
        domain: _plan_cryosparc_group(policy, domain, source_info)
        for domain in domains
        if formats[domain] == "cryosparc_cs"
    }
    inputs = [run_dir, selection_dir] + [Path(v["path"]) for v in source_info.values() if v.get("path")]
    inputs += [path for info in source_info.values() for path in info.get('protected_input_paths', ())]
    protect_input_directory(output_dir, inputs)
    _validate_output_dir(output_dir, overwrite=policy.overwrite)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.export-", dir=output_dir.parent))
    staging = scratch / "export"
    backup = scratch / "previous"
    publish_lock = output_dir.parent / f".{output_dir.name}.export-publish"
    staging.mkdir()
    published = recovery_required = False
    try:
        report = _write_metadata_export(
            staging, output_dir, selection=selection, selected_rows=selected_rows, policy=policy, domains=domains,
            run_dir=run_dir, run_id=run_id, selection_path_obj=selection_path_obj, source_info=source_info,
            source_verification=source_verification, formats=formats, group_plans=group_plans,
        )
        try:
            publish_lock.mkdir()
        except FileExistsError as exc:
            raise FileExistsError(f"Another export is publishing to {output_dir}") from exc
        try:
            _validate_output_dir(output_dir, overwrite=policy.overwrite)
            if policy.overwrite and output_dir.exists():
                _keep_existing_files(output_dir, staging)
            protect_input_directory(output_dir, inputs)
            try:
                source_guard.assert_unchanged()
            except (OSError, ValueError) as error:
                raise ValueError('Source metadata changed or became unavailable during export; no metadata was published.') from error
            try:
                if policy.overwrite and output_dir.exists():
                    os.replace(output_dir, backup)
                os.replace(staging, output_dir)
            except BaseException:
                try:
                    if backup.exists():
                        if output_dir.exists():
                            shutil.rmtree(output_dir)
                        os.replace(backup, output_dir)
                    elif not staging.exists() and output_dir.exists():
                        shutil.rmtree(output_dir)
                except BaseException as error:
                    recovery_required = True
                    raise RuntimeError(f'Export rollback failed; recover preserved files in {scratch}') from error
                raise
            published = True
        finally:
            if not recovery_required:
                try:
                    publish_lock.rmdir()
                except OSError as error:
                    if not published:
                        raise
                    warnings.warn(f'Export published; leftover lock needs cleanup: {error}', RuntimeWarning)
        return report
    finally:
        if not recovery_required and (published or not backup.exists()):
            shutil.rmtree(scratch, ignore_errors=True)


def _plan_cryosparc_group(
    policy: SelectionMetadataExportPolicy,
    domain: str,
    source_info: Mapping[str, Mapping[str, Any]],
) -> tuple[Any, str, list[str]]:
    """(group or None, resolved_by, candidates) for one CryoSPARC domain, checked; raises before writing."""

    explicit = policy.cryosparc_group_ref if domain == "ref" else policy.cryosparc_group_mov
    if not policy.mirror_cryosparc_group:
        if explicit is not None:
            raise ValueError("--no-cryosparc-group cannot be combined with --cryosparc-group-ref/-mov.")
        return None, "disabled", []
    source_path = source_info.get(domain, {}).get("path")
    if not source_path:
        raise ValueError(f"Cannot locate original {domain} source metadata path")
    group, resolved_by, candidates = resolve_group(source_path, explicit)
    if group is not None:
        expected = source_info.get(domain, {}).get("row_count")
        try:
            check_cs_group(group, source_path, int(expected) if expected is not None else None)
        except ValueError as exc:
            hint = (
                "" if resolved_by == "explicit"
                else " Choose the right group with --cryosparc-group-ref/-mov, or export the .cs subset "
                "alone with --no-cryosparc-group."
            )
            raise ValueError(f"{exc}{hint}") from exc
    return group, resolved_by, candidates


def _write_metadata_export(
    staging: Path,
    output_dir: Path,
    *,
    selection: Selection,
    selected_rows: pd.DataFrame,
    policy: SelectionMetadataExportPolicy,
    domains: tuple[str, ...],
    run_dir: Path | None,
    run_id: str | None,
    selection_path_obj: Path | None,
    source_info: Mapping,
    source_verification: Mapping,
    formats: Mapping[str, str],
    group_plans: Mapping[str, tuple[Any, str, list[str]]],
) -> dict[str, Any]:
    source_row_id_fields = {domain: DOMAIN_ROW_ID_COLUMNS[domain] for domain in domains}
    source_row_id_stats = {
        domain: _source_row_id_stats(selected_rows, row_id_field)
        for domain, row_id_field in source_row_id_fields.items()
    }
    timestamp = datetime.now(timezone.utc).isoformat()

    selected_keys_path = output_dir / "selected_particle_keys.txt"
    selected_rows_path = output_dir / "selected_landscape_rows.csv"
    _write_selected_keys_txt(staging / selected_keys_path.name, selection.selected_particle_keys)
    selected_rows.to_csv(staging / selected_rows_path.name, index=False)

    format_resolved: dict[str, str | None] = {"ref": None, "mov": None}
    outputs: dict[str, str | None] = {
        "selected_particle_keys_txt": str(selected_keys_path),
        "selected_landscape_rows_csv": str(selected_rows_path),
        "ref_metadata": None,
        "mov_metadata": None,
    }
    domain_reports: dict[str, str | None] = {"ref": None, "mov": None}
    next_steps: dict[str, list[str]] = {}
    warnings: list[str] = []

    for domain in domains:
        resolved_format = formats[domain]
        format_resolved[domain] = resolved_format
        if resolved_format == "keys":
            continue
        domain_dir = staging / domain
        domain_dir.mkdir(parents=True, exist_ok=True)
        row_id_column = DOMAIN_ROW_ID_COLUMNS[domain]
        row_ids = _validated_domain_row_ids(selected_rows, row_id_column)
        domain_report = _export_domain_metadata(
            domain=domain,
            resolved_format=resolved_format,
            source_info=source_info,
            row_ids=row_ids,
            row_id_column=row_id_column,
            domain_dir=domain_dir,
            overwrite=policy.overwrite,
            group_plan=group_plans.get(domain),
            selection_id=selection.selection_id,
            run_id=run_id,
        )
        domain_report = _public_paths(domain_report, staging, output_dir)
        next_steps[domain] = list(domain_report.get("next_steps", ()))
        report_path = domain_dir / f"export_{domain}_report.json"
        _write_json(report_path, domain_report)
        domain_reports[domain] = str(output_dir / domain / report_path.name)
        outputs[f"{domain}_metadata"] = domain_report.get("output_path")
        warnings.extend(str(warning) for warning in domain_report.get("warnings", ()))

    export_report = {
        "artifact_type": "cryorole_metadata_subset_export",
        "schema_version": "2.0",
        "timestamp": timestamp,
        "run_dir": str(run_dir) if run_dir is not None else None,
        "selection_id": selection.selection_id,
        "parent_run_id": selection.parent_run_id,
        "run_id": run_id,
        "selection_path": str(selection_path_obj) if selection_path_obj is not None else None,
        "resolved_by": dict(policy.resolved_by or {}),
        "domain": policy.domain,
        "format_requested": policy.format,
        "format_resolved": format_resolved,
        "output_dir": str(output_dir),
        "overwrite": bool(policy.overwrite),
        "selected_count": selection.selected_count,
        "outputs": outputs,
        "domain_reports": domain_reports,
        "source_row_id_fields": source_row_id_fields,
        "source_row_id_stats": source_row_id_stats,
        "source_files": {
            domain: source_info.get(domain, {}).get("path")
            for domain in ("ref", "mov")
        },
        "source_verification": source_verification,
        "allow_unverified_source": bool(policy.allow_unverified_source),
        "export_policy": to_json_safe(policy),
        "next_steps": next_steps,
        "warnings": warnings,
    }
    _write_json(staging / "export_report.json", export_report)
    _write_next_steps(staging / "NEXT_STEPS.txt", next_steps)
    return export_report


def preview_metadata_export(run_dir: str | Path, *, domain: str = "both") -> dict[str, Any]:
    """What an export of this run would use, per domain, without writing or hashing anything.

    Reports the recorded source file and format and, for CryoSPARC sources, the
    result group (.csg) that would be mirrored (only candidate next to the source).
    """

    run = Path(run_dir)
    validate_completed_run_bundle(run)
    info = _resolve_source_info(run)
    preview: dict[str, Any] = {}
    for name in _domains(domain):
        entry = dict(info.get(name, {}))
        entry["cryosparc_group"] = None
        entry["cryosparc_group_note"] = None
        if entry.get("source_type") == "cryosparc_cs" and entry.get("path"):
            try:
                group, resolved_by, candidates = resolve_group(entry["path"], None)
            except ValueError as exc:
                entry["cryosparc_group_note"] = str(exc)
            else:
                entry["cryosparc_group"] = str(group.path) if group else None
                entry["cryosparc_group_resolved_by"] = resolved_by
                if group is not None:
                    try:
                        count = entry.get("row_count")
                        check_cs_group(group, entry["path"], int(count) if count is not None else None)
                    except (OSError, ValueError) as exc:
                        entry["cryosparc_group_usable"] = False
                        entry["cryosparc_group_note"] = f"This result group would stop before writing: {exc}"
                    else:
                        entry["cryosparc_group_usable"] = True
                if group is None:
                    entry["cryosparc_group_note"] = _NO_GROUP_NOTE.format(name=Path(entry["path"]).name)
                else:
                    try:  # the checks an export makes before writing, so a problem shows here first
                        check_cs_group(group, entry["path"], entry.get("row_count"))
                    except ValueError as exc:
                        entry["cryosparc_group_note"] = f"{exc} The export would stop before writing."
                        entry["cryosparc_group_usable"] = False
        preview[name] = entry
    return preview


_NO_GROUP_NOTE = (
    "No CryoSPARC result group (.csg) names {name} next to it, so the .cs subset cannot go back through "
    "Import Result Group. Point to the source job's particle group, e.g. <project>/J75/J75_particles.csg "
    "(--cryosparc-group-ref / --cryosparc-group-mov)."
)


def _write_next_steps(path: Path, next_steps: Mapping[str, list[str]]) -> None:
    lines = []
    for domain, steps in next_steps.items():
        if steps:
            lines.append(f"[{domain}]")
            lines.extend(steps)
            lines.append("")
    if lines:
        path.write_text("\n".join(lines), encoding="utf-8")


def _validate_policy(policy: SelectionMetadataExportPolicy) -> None:
    if policy.domain not in {"ref", "mov", "both"}:
        raise ValueError("--domain must be ref, mov, or both")
    if policy.format not in {"auto", "relion_star", "cryosparc_cs", "keys"}:
        raise ValueError("--format must be auto, relion_star, cryosparc_cs, or keys")


def _domains(domain: str) -> tuple[str, ...]:
    if domain == "both":
        return ("ref", "mov")
    return (domain,)


def _resolve_output_dir(
    *,
    policy: SelectionMetadataExportPolicy,
    run_dir: Path | None,
    selection: Selection,
) -> Path:
    if policy.output_dir is not None:
        return Path(policy.output_dir)
    if run_dir is None:
        raise ValueError("Direct metadata export requires --output-dir when --run-dir is absent")
    return run_dir / "exports" / selection.selection_id


def _validate_output_dir(path: Path, *, overwrite: bool) -> None:
    if path.is_symlink() or (path.exists() and not path.is_dir()):
        raise ValueError(f"Export output path must be a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(f"Export output directory already exists: {path}")


def _keep_existing_files(source: Path, staging: Path) -> None:
    """Preserve existing entries that this export does not replace."""
    for entry in source.iterdir():
        target = staging / entry.name
        if target.exists() or target.is_symlink():
            if entry.is_dir() and not entry.is_symlink() and target.is_dir() and not target.is_symlink():
                _keep_existing_files(entry, target)
        elif entry.is_dir() and not entry.is_symlink():
            shutil.copytree(entry, target, symlinks=True)
        else:
            shutil.copy2(entry, target, follow_symlinks=False)


def _public_paths(value: Any, staging: Path, output_dir: Path) -> Any:
    if isinstance(value, str):
        return value.replace(str(staging), str(output_dir))
    if isinstance(value, dict):
        return {key: _public_paths(item, staging, output_dir) for key, item in value.items()}
    if isinstance(value, list):
        return [_public_paths(item, staging, output_dir) for item in value]
    return value


def _load_selected_rows(
    *,
    selection: Selection,
    selection_dir: Path | None,
    run_dir: Path | None,
    domains: tuple[str, ...],
    source_info: Mapping,
) -> pd.DataFrame:
    selected_rows_path = selection_dir / "selected_landscape_rows.csv" if selection_dir else None
    companion = None
    if selected_rows_path is not None and selected_rows_path.exists():
        rows = pd.read_csv(selected_rows_path, dtype={"particle_key": str}, keep_default_na=False)
        companion = _reorder_selected_rows(rows, selection)
    if run_dir is None:
        if companion is None:
            raise ValueError(
                "Selected landscape rows are required for metadata export; provide "
                "--run-dir or a selection directory containing selected_landscape_rows.csv"
            )
        return companion

    trusted = _selected_rows_from_run_landscape(run_dir, selection)
    if companion is None:
        return trusted
    for domain, field in DOMAIN_ROW_ID_COLUMNS.items():
        if domain not in domains and field not in companion:
            continue
        observed = _validated_domain_row_ids(companion, field)
        expected = _validated_domain_row_ids(trusted, field)
        count = source_info.get(domain, {}).get("row_count")
        if count is not None and any(i >= count for i in observed):
            raise ValueError(f"{selected_rows_path}: {field} has out-of-bounds source row IDs")
        if observed != expected:
            raise ValueError(
                f"{selected_rows_path}: {field} differs from the parent run's particle-to-source-row mapping"
            )
    return companion


def _selected_rows_from_run_landscape(run_dir: Path, selection: Selection) -> pd.DataFrame:
    landscape_path = resolve_landscape_path(run_dir, space="raw")
    if landscape_path.suffix.lower() == ".npz":
        # Export needs only the key and source-row arrays, even for a multi-million-particle run.
        with np.load(landscape_path, allow_pickle=False) as data:
            if str(data["schema_version"].item()) != LANDSCAPE_NPZ_SCHEMA_VERSION:
                raise ValueError(f"Unsupported parent landscape schema: {landscape_path}")
            keys = data["particle_key"]
            columns = {"particle_key": keys.astype(str)}
            for field in DOMAIN_ROW_ID_COLUMNS.values():
                if field in data:
                    ids = data[field]
                    if ids.shape != keys.shape or ids.ndim != 1 or ids.dtype.kind not in "iu":
                        raise ValueError(f"Invalid {field} in parent landscape: {landscape_path}")
                    columns[field] = ids
        table = pd.DataFrame(columns)
    elif landscape_path.suffix.lower() == ".csv":
        fields = {"particle_key", *DOMAIN_ROW_ID_COLUMNS.values()}
        table = pd.read_csv(landscape_path, usecols=lambda name: name in fields,
                            dtype={"particle_key": str}, keep_default_na=False)
    else:
        raise ValueError(
            "Metadata subset export can reconstruct selected rows from raw "
            "NPZ/CSV landscapes only"
        )
    return _reorder_selected_rows(table, selection)


def _reorder_selected_rows(rows: pd.DataFrame, selection: Selection) -> pd.DataFrame:
    if "particle_key" not in rows.columns:
        raise ValueError("selected_landscape_rows.csv must include particle_key")
    duplicated = rows["particle_key"].astype(str).duplicated()
    if duplicated.any():
        duplicates = rows.loc[duplicated, "particle_key"].astype(str).head().tolist()
        raise ValueError(f"Selected landscape rows contain duplicate particle_key values: {duplicates}")
    indexed = rows.assign(particle_key=rows["particle_key"].astype(str)).set_index("particle_key")
    missing = [str(key) for key in selection.selected_particle_keys if str(key) not in indexed.index]
    if missing:
        raise ValueError(
            "Selected landscape rows are missing selected particle keys: "
            f"{missing[:5]}"
        )
    ordered = indexed.loc[[str(key) for key in selection.selected_particle_keys]].reset_index()
    return ordered


def _validated_domain_row_ids(rows: pd.DataFrame, row_id_column: str) -> list[int]:
    if row_id_column not in rows.columns:
        raise ValueError(f"Selected rows lack required {row_id_column} column")
    values = rows[row_id_column]
    if values.isna().any():
        raise ValueError(f"{row_id_column} contains missing source row IDs")
    row_ids: list[int] = []
    for value in values:
        if isinstance(value, bool):
            raise ValueError(f"{row_id_column} must contain integer source row IDs")
        try:
            numeric = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{row_id_column} must contain integer source row IDs") from exc
        if not np.isfinite(numeric) or not numeric.is_integer():
            raise ValueError(f"{row_id_column} must contain integer source row IDs")
        row_ids.append(int(numeric))
    if any(row_id < 0 for row_id in row_ids):
        raise ValueError(f"{row_id_column} contains negative source row IDs")
    if len(set(row_ids)) != len(row_ids):
        raise ValueError(f"{row_id_column} contains duplicate source row IDs")
    return row_ids


def _source_row_id_stats(rows: pd.DataFrame, row_id_column: str) -> dict[str, int | None]:
    row_ids = _validated_domain_row_ids(rows, row_id_column)
    if not row_ids:
        return {
            "count": 0,
            "min": None,
            "max": None,
            "unique_count": 0,
        }
    return {
        "count": len(row_ids),
        "min": min(row_ids),
        "max": max(row_ids),
        "unique_count": len(set(row_ids)),
    }


def _resolve_source_info(run_dir: Path) -> dict[str, dict[str, str | None]]:
    summary_path = run_dir / "run_summary.json"
    manifest_path = run_dir / "run_manifest.json"
    summary = _read_json_if_exists(summary_path)
    manifest = _read_json_if_exists(manifest_path)
    input_paths = _domain_mapping(summary.get("input_paths")) or _manifest_input_paths(manifest)
    source_types = _domain_mapping(summary.get("source_types")) or _manifest_source_types(
        manifest,
        input_paths,
    )
    info: dict[str, dict[str, str | None]] = {}
    for domain in ("ref", "mov"):
        path = input_paths.get(domain)
        source_type = source_types.get(domain) or _source_type_from_path(path)
        if path is None:
            info[domain] = {"path": None, "source_type": None}
            continue
        resolved_path = _resolve_source_path(path, run_dir)
        info[domain] = {
            "path": str(resolved_path),
            "source_type": _metadata_format_from_source_type(source_type, resolved_path),
        }
    return info


def _resolve_and_verify_sources(
    run_dir: Path,
    *,
    policy: SelectionMetadataExportPolicy,
    domains: tuple[str, ...],
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]], str | None]:
    summary = _read_json_if_exists(run_dir / "run_summary.json")
    manifest = _read_json_if_exists(run_dir / "run_manifest.json")
    run_id_value = summary.get("run_id") or manifest.get("run_id")
    run_id = str(run_id_value) if run_id_value is not None else None
    identities = summary.get("source_identities") or manifest.get("source_identities") or {}
    if not isinstance(identities, Mapping):
        identities = {}
    legacy_info = _resolve_source_info(run_dir)
    source_info: dict[str, dict[str, Any]] = dict(legacy_info)
    verification: dict[str, dict[str, Any]] = {}
    relocated = {"ref": policy.relocated_ref, "mov": policy.relocated_mov}
    for domain in ("ref", "mov"):
        record = identities.get(domain)
        paths = []
        if isinstance(record, Mapping):
            for field in ('original_path', 'input_path_absolute', 'resolved_path', 'path'):
                value = record.get(field)
                # Old relative spellings have no recorded cwd. Do not invent a
                # historical location from the current export working directory.
                if value and Path(str(value)).expanduser().is_absolute():
                    paths.append(str(Path(str(value)).expanduser()))
        current = relocated[domain]
        if current is not None:
            paths.append(str(Path(current).expanduser().absolute()))
        source_info.setdefault(domain, {})['protected_input_paths'] = paths
    for domain in domains:
        record = identities.get(domain) if isinstance(identities, Mapping) else None
        if record is None:
            legacy = legacy_info.get(domain, {})
            legacy_path = legacy.get("path")
            record = {
                "path": legacy_path,
                "resolved_path": legacy_path,
                "source_type": legacy.get("source_type"),
            }
        result = verify_source_identity(
            record,
            relocated_path=relocated[domain],
            allow_unverified_source=policy.allow_unverified_source,
        )
        source_type = (
            str(record.get("source_type"))
            if isinstance(record, Mapping) and record.get("source_type") is not None
            else legacy_info.get(domain, {}).get("source_type")
        )
        resolved_path = Path(result.resolved_path)
        recorded_row_count = record.get("row_count") if isinstance(record, Mapping) else None
        expected_hash = record.get('sha256') or record.get('hash')
        # Explicit legacy trust cannot prove old bytes. Bind the current read
        # without relabeling it as historically verified.
        export_hash = str(expected_hash) if expected_hash else stream_sha256(resolved_path)
        source_info[domain] = {
            "path": result.resolved_path,
            "source_type": _metadata_format_from_source_type(source_type, resolved_path),
            "row_count": int(recorded_row_count) if recorded_row_count is not None else None,
            "sha256": export_hash,
            "protected_input_paths": source_info[domain]['protected_input_paths'],
            "export_input_path_absolute": (str(Path(relocated[domain]).expanduser().absolute())
                                           if relocated[domain] is not None else None),
        }
        verification[domain] = {
            "resolved_path": result.resolved_path,
            "verification_status": result.verification_status,
            "relocated": result.relocated,
            "sha256_verified": result.sha256_verified,
            "legacy_unverified_allowed": result.legacy_unverified_allowed,
            "export_source_sha256": export_hash,
        }
    return source_info, verification, run_id


def _read_json_if_exists(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Run provenance JSON must be an object: {path}")
    return payload


def _domain_mapping(value: Any) -> dict[str, str]:
    if not isinstance(value, Mapping):
        return {}
    result = {}
    for domain in ("ref", "mov"):
        if domain in value and value[domain] is not None:
            result[domain] = str(value[domain])
    return result


def _manifest_input_paths(manifest: dict[str, Any]) -> dict[str, str]:
    input_files = _domain_mapping(manifest.get("input_files"))
    if input_files:
        return input_files
    records = manifest.get("input_provenance", {}).get("inputs", [])
    if isinstance(records, list) and len(records) >= 2:
        paths = []
        for record in records[:2]:
            if isinstance(record, Mapping) and record.get("path") is not None:
                paths.append(str(record["path"]))
        if len(paths) == 2:
            return {"ref": paths[0], "mov": paths[1]}
    return {}


def _manifest_source_types(
    manifest: dict[str, Any],
    input_paths: Mapping[str, str],
) -> dict[str, str]:
    domain_source_types = _domain_mapping(manifest.get("source_types"))
    if domain_source_types:
        return domain_source_types
    source_types = manifest.get("input_provenance", {}).get("source_types", {})
    if not isinstance(source_types, Mapping):
        return {}
    by_domain = {}
    for domain, path in input_paths.items():
        if path in source_types:
            by_domain[domain] = str(source_types[path])
    return by_domain


def _resolve_source_path(path: str, run_dir: Path) -> Path:
    candidate = Path(path)
    if candidate.exists():
        return candidate
    relative = run_dir / candidate
    if relative.exists():
        return relative
    return candidate


def _source_type_from_path(path: str | None) -> str | None:
    if path is None:
        return None
    suffix = Path(path).suffix.lower()
    if suffix == ".star":
        return "relion_star"
    if suffix == ".cs":
        return "cryosparc_cs"
    return None


def _metadata_format_from_source_type(
    source_type: str | None,
    path: Path,
) -> str | None:
    normalized = (source_type or "").lower()
    if normalized in {"relion", "relion_star", "star"}:
        return "relion_star"
    if normalized in {"cryosparc", "cryosparc_cs", "cs"}:
        return "cryosparc_cs"
    return _source_type_from_path(str(path))


def _resolve_format_for_domain(
    *,
    requested: str,
    domain: str,
    source_info: Mapping[str, Mapping[str, str | None]],
) -> str:
    if requested == "keys":
        return "keys"
    source_format = source_info.get(domain, {}).get("source_type")
    if requested == "auto":
        if source_format in {"relion_star", "cryosparc_cs"}:
            return source_format
        raise ValueError(
            f"Could not resolve source format for {domain}; provide run provenance "
            "or use --format keys"
        )
    if requested != source_format:
        raise ValueError(
            f"Forced --format {requested} is incompatible with {domain} source "
            f"format {source_format!r}"
        )
    return requested


def _export_domain_metadata(
    *,
    domain: str,
    resolved_format: str,
    source_info: Mapping[str, Mapping[str, str | None]],
    row_ids: list[int],
    row_id_column: str,
    domain_dir: Path,
    overwrite: bool,
    group_plan: tuple[Any, str, list[str]] | None = None,
    selection_id: str = "",
    run_id: str | None = None,
) -> dict[str, Any]:
    source_path = source_info.get(domain, {}).get("path")
    if not source_path:
        raise ValueError(f"Cannot locate original {domain} source metadata path")
    expected_row_count = source_info.get(domain, {}).get("row_count")
    if resolved_format == "relion_star":
        output_path = domain_dir / f"selected_{domain}.star"
        result = write_relion_star_subset(
            source_path,
            output_path,
            row_ids,
            row_id_field=row_id_column,
            overwrite=overwrite,
            expected_source_row_count=expected_row_count,
            required_headers=RELION_POSE_HEADERS,
        )
        uid_output_path = None
        group_report = None
        next_steps = [
            f"{output_path} is a verbatim subset of {Path(source_path).name} (optics and CTF unchanged); "
            "use it directly in RELION, for example:",
            f"  relion_reconstruct --i {output_path} --o {domain}_subset.mrc --ctf",
        ]
    elif resolved_format == "cryosparc_cs":
        output_path = domain_dir / f"selected_{domain}.cs"
        uid_path = domain_dir / f"selected_{domain}_uids.txt"
        result = write_cryosparc_cs_subset(
            source_path,
            output_path,
            row_ids,
            row_id_field=row_id_column,
            overwrite=overwrite,
            uid_output_path=uid_path,
            expected_source_row_count=expected_row_count,
        )
        uid_output_path = result.uid_output_path
        group, resolved_by, candidates = group_plan if group_plan is not None else resolve_group(source_path, None)
        if group is None and resolved_by == "disabled":
            group_report = {"resolved_by": "disabled", "candidates": []}
            next_steps = [f"{output_path} was written without a CryoSPARC result group (--no-cryosparc-group)."]
        elif group is None:
            group_report = {"resolved_by": "none", "candidates": candidates}
            note = _NO_GROUP_NOTE.format(name=Path(source_path).name)
            next_steps = [note]
            result = replace(result, warnings=(*result.warnings, note))
        else:
            written = write_cs_group_subset(
                group,
                resolved_by=resolved_by,
                source_path=source_path,
                source_row_count=result.source_row_count,
                main_subset_path=output_path,
                domain=domain,
                domain_dir=domain_dir,
                selection_id=selection_id,
                run_id=run_id,
                overwrite=overwrite,
            )
            group_report = asdict(written)
            next_steps = list(written.import_steps)
    else:
        raise ValueError(f"Unsupported resolved export format: {resolved_format}")

    return {
        "artifact_type": "cryorole_metadata_subset_domain_export",
        "schema_version": "2.0",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "domain": domain,
        "source_path": result.source_path,
        "output_path": result.output_path,
        "source_type": resolved_format,
        "row_id_column": row_id_column,
        "selected_row_count": result.selected_row_count,
        "source_row_count": result.source_row_count,
        "duplicate_row_id_count": 0,
        "missing_row_ids": [],
        "out_of_bounds_row_ids": [],
        "uid_output_path": uid_output_path,
        "source_row_count_verified": expected_row_count is not None,
        "source_file_unchanged": True,
        "cryosparc_group": group_report,
        "next_steps": next_steps,
        "warnings": list(result.warnings),
    }


def _write_selected_keys_txt(path: Path, keys: tuple[object, ...]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for key in keys:
            handle.write(f"{key}\n")


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(to_json_safe(payload), handle, indent=2, sort_keys=True)
        handle.write("\n")
