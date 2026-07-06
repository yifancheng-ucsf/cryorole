"""Direct source metadata subset export for selected cryoROLE particles."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from cryorole.export.cs_subset import write_cryosparc_cs_subset
from cryorole.export.serialization import to_json_safe
from cryorole.export.star_subset import write_relion_star_subset
from cryorole.io.writers.landscape_store import (
    read_landscape_csv,
    read_landscape_npz_arrays,
    resolve_landscape_path,
)
from cryorole.models.policies import SelectionMetadataExportPolicy
from cryorole.models.selection import Selection


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
    selection_path_obj = Path(selection_path).resolve() if selection_path else None
    selection_dir = selection_path_obj.parent if selection_path_obj else None
    output_dir = _resolve_output_dir(
        policy=policy,
        run_dir=run_dir,
        selection=selection,
    )
    _prepare_output_dir(output_dir, overwrite=policy.overwrite)

    selected_rows = _load_selected_rows(
        selection=selection,
        selection_dir=selection_dir,
        run_dir=run_dir,
    )
    domains = _domains(policy.domain)
    source_row_id_fields = {
        domain: DOMAIN_ROW_ID_COLUMNS[domain]
        for domain in domains
        if DOMAIN_ROW_ID_COLUMNS[domain] in selected_rows.columns
    }
    source_row_id_stats = {
        domain: _source_row_id_stats(selected_rows, DOMAIN_ROW_ID_COLUMNS[domain])
        for domain in domains
    }
    for domain in domains:
        _validated_domain_row_ids(selected_rows, DOMAIN_ROW_ID_COLUMNS[domain])

    source_info = _resolve_source_info(run_dir) if run_dir is not None else {}
    timestamp = datetime.now(timezone.utc).isoformat()

    selected_keys_path = output_dir / "selected_particle_keys.txt"
    selected_rows_path = output_dir / "selected_landscape_rows.csv"
    _write_selected_keys_txt(selected_keys_path, selection.selected_particle_keys)
    selected_rows.to_csv(selected_rows_path, index=False)

    format_resolved: dict[str, str | None] = {"ref": None, "mov": None}
    outputs: dict[str, str | None] = {
        "selected_particle_keys_txt": str(selected_keys_path),
        "selected_landscape_rows_csv": str(selected_rows_path),
        "ref_metadata": None,
        "mov_metadata": None,
    }
    domain_reports: dict[str, str | None] = {"ref": None, "mov": None}
    warnings: list[str] = []

    for domain in domains:
        resolved_format = _resolve_format_for_domain(
            requested=policy.format,
            domain=domain,
            source_info=source_info,
        )
        format_resolved[domain] = resolved_format
        if resolved_format == "keys":
            continue
        domain_dir = output_dir / domain
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
        )
        report_path = domain_dir / f"export_{domain}_report.json"
        _write_json(report_path, domain_report)
        domain_reports[domain] = str(report_path)
        outputs[f"{domain}_metadata"] = domain_report.get("output_path")
        warnings.extend(str(warning) for warning in domain_report.get("warnings", ()))

    export_report = {
        "artifact_type": "cryorole_metadata_subset_export",
        "schema_version": "2.0",
        "timestamp": timestamp,
        "run_dir": str(run_dir) if run_dir is not None else None,
        "selection_id": selection.selection_id,
        "selection_path": str(selection_path_obj) if selection_path_obj is not None else None,
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
        "export_policy": to_json_safe(policy),
        "warnings": warnings,
    }
    _write_json(output_dir / "export_report.json", export_report)
    return export_report


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


def _prepare_output_dir(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not path.is_dir():
        raise ValueError(f"Export output path exists and is not a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(f"Export output directory already exists: {path}")
    path.mkdir(parents=True, exist_ok=True)


def _load_selected_rows(
    *,
    selection: Selection,
    selection_dir: Path | None,
    run_dir: Path | None,
) -> pd.DataFrame:
    selected_rows_path = (
        selection_dir / "selected_landscape_rows.csv"
        if selection_dir is not None
        else None
    )
    if selected_rows_path is not None and selected_rows_path.exists():
        rows = pd.read_csv(selected_rows_path)
        return _reorder_selected_rows(rows, selection)
    if run_dir is None:
        raise ValueError(
            "Selected landscape rows are required for metadata export; provide "
            "--run-dir or a selection directory containing selected_landscape_rows.csv"
        )
    return _selected_rows_from_run_landscape(run_dir, selection)


def _selected_rows_from_run_landscape(run_dir: Path, selection: Selection) -> pd.DataFrame:
    landscape_path = resolve_landscape_path(run_dir, space="raw")
    if landscape_path.suffix.lower() == ".npz":
        arrays = read_landscape_npz_arrays(landscape_path)
        table = pd.DataFrame(
            {
                "particle_key": arrays.particle_key.astype(str),
                "ref_source_row_id": (
                    arrays.ref_source_row_id
                    if arrays.ref_source_row_id is not None
                    else np.full(arrays.n_points, -1, dtype=np.int64)
                ),
                "mov_source_row_id": (
                    arrays.mov_source_row_id
                    if arrays.mov_source_row_id is not None
                    else np.full(arrays.n_points, -1, dtype=np.int64)
                ),
            }
        )
    elif landscape_path.suffix.lower() == ".csv":
        landscape = read_landscape_csv(landscape_path)
        table = landscape.data[
            [
                column
                for column in ("particle_key", "ref_source_row_id", "mov_source_row_id")
                if column in landscape.data.columns
            ]
        ].copy()
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
) -> dict[str, Any]:
    source_path = source_info.get(domain, {}).get("path")
    if not source_path:
        raise ValueError(f"Cannot locate original {domain} source metadata path")
    if resolved_format == "relion_star":
        output_path = domain_dir / f"selected_{domain}.star"
        result = write_relion_star_subset(
            source_path,
            output_path,
            row_ids,
            row_id_field=row_id_column,
            overwrite=overwrite,
        )
        uid_output_path = None
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
        )
        uid_output_path = result.uid_output_path
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
        "source_file_unchanged": True,
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
