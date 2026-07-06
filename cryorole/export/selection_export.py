"""Non-destructive export of Selection objects."""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from json import JSONDecodeError
from pathlib import Path
from typing import Iterable

from cryorole.export.serialization import to_json_safe
from cryorole.models.export_report import ExportReport
from cryorole.models.policies import SelectionExportPolicy, SelectionPolicy
from cryorole.models.selection import Selection


REQUIRED_SELECTION_JSON_FIELDS = (
    "selection_id",
    "selection_mode",
    "selection_basis",
    "metric",
    "selected_particle_keys",
    "selected_count",
    "total_count",
)


def export_selection(
    selection: Selection,
    *,
    policy: SelectionExportPolicy,
) -> ExportReport:
    """Export selected particle keys and selection metadata without source mutation."""

    _validate_export_policy(policy)
    output_dir = Path(policy.output_dir)
    _validate_simple_filename(policy.selected_keys_filename, "selected_keys_filename")
    _validate_simple_filename(policy.selection_json_filename, "selection_json_filename")
    _validate_simple_filename(policy.export_report_filename, "export_report_filename")
    selected_keys_path = output_dir / policy.selected_keys_filename
    selection_json_path = output_dir / policy.selection_json_filename
    export_report_path = output_dir / policy.export_report_filename
    _prepare_output_paths(
        (selected_keys_path, selection_json_path, export_report_path),
        overwrite=policy.overwrite,
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    _write_selected_keys_csv(selected_keys_path, selection)
    _write_selection_json(selection_json_path, selection, policy)

    timestamp = datetime.now(timezone.utc).isoformat()
    report = ExportReport(
        output_paths={
            "selected_particle_keys_csv": str(selected_keys_path),
            "selection_json": str(selection_json_path),
            "export_report_json": str(export_report_path),
        },
        row_counts={
            "selected_particle_keys_csv": selection.selected_count,
            "selection_json": 1,
            "export_report_json": 1,
        },
        source_selection_id=selection.selection_id,
        timestamp=timestamp,
        export_policy=to_json_safe(policy),
    )
    _write_export_report_json(export_report_path, report)
    return report


def read_selection_json(path: str | Path) -> Selection:
    """Read a previously exported Selection JSON artifact.

    Re-export requires ``selected_particle_keys`` to be present in the JSON.
    Selection JSON files written with
    ``include_selected_particle_keys_in_json=False`` are not re-exportable until
    a future paired selected-particle-key CSV reader is implemented.
    """

    input_path = Path(path)
    if not input_path.exists():
        raise ValueError(f"Selection JSON file does not exist: {input_path}")
    try:
        with input_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except JSONDecodeError as exc:
        raise ValueError(f"Malformed selection JSON: {input_path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("Selection JSON must be an object")
    _validate_selection_artifact_header(payload)
    missing = [field for field in REQUIRED_SELECTION_JSON_FIELDS if field not in payload]
    if missing:
        raise ValueError(f"Selection JSON missing required fields: {missing}")

    selected_particle_keys = payload["selected_particle_keys"]
    if not isinstance(selected_particle_keys, list):
        raise ValueError("Selection JSON selected_particle_keys must be a list")
    selected_count = _non_negative_int(payload["selected_count"], "selected_count")
    total_count = _non_negative_int(payload["total_count"], "total_count")
    if len(selected_particle_keys) != selected_count:
        raise ValueError(
            "Selection JSON selected_particle_keys length must equal selected_count"
        )
    if selected_count > total_count:
        raise ValueError("Selection JSON selected_count must be <= total_count")
    parent_landscape_metadata = payload.get("parent_landscape_metadata")
    if parent_landscape_metadata is None:
        parent_landscape_metadata = {}
    elif not isinstance(parent_landscape_metadata, dict):
        raise ValueError("Selection JSON parent_landscape_metadata must be an object")
    parent_landscape_id = payload.get("parent_landscape_id")
    if parent_landscape_id is not None and not isinstance(parent_landscape_id, str):
        raise ValueError("Selection JSON parent_landscape_id must be a string or null")

    active_policy_payload = payload.get("active_policy")
    active_policy = None
    if active_policy_payload is not None:
        if not isinstance(active_policy_payload, dict):
            raise ValueError("Selection JSON active_policy must be an object when present")
        active_policy = _selection_policy_from_payload(active_policy_payload)

    return Selection(
        selection_id=payload["selection_id"],
        parent_landscape_id=parent_landscape_id,
        parent_landscape_metadata=parent_landscape_metadata,
        selection_mode=payload["selection_mode"],
        selection_basis=payload["selection_basis"],
        metric=payload["metric"],
        center_input=_optional_tuple(payload.get("center_input")),
        center_input_representation=payload.get("center_input_representation"),
        center_input_space=payload.get("center_input_space"),
        center_euler_sequence=payload.get("center_euler_sequence"),
        center_euler_convention=payload.get("center_euler_convention"),
        center_scipy_euler_sequence=payload.get("center_scipy_euler_sequence"),
        center_degrees=payload.get("center_degrees"),
        center_evaluated=_optional_tuple(payload.get("center_evaluated")),
        evaluation_space=payload.get("evaluation_space"),
        resolved_evaluation_space=payload.get("resolved_evaluation_space"),
        radius=payload.get("radius"),
        radius_unit=payload.get("radius_unit"),
        threshold=payload.get("threshold"),
        sld_min=payload.get("sld_min"),
        sld_max=payload.get("sld_max"),
        threshold_operator=payload.get("threshold_operator"),
        density_support_field=payload.get("density_support_field"),
        density_artifact_policy=_density_artifact_policy_from_payload(
            payload.get("density_artifact_policy", "include_all")
        ),
        density_artifact_flag_field=payload.get("density_artifact_flag_field"),
        density_artifact_candidate_count_before=payload.get(
            "density_artifact_candidate_count_before"
        ),
        density_artifact_candidate_count_after=payload.get(
            "density_artifact_candidate_count_after"
        ),
        density_artifact_excluded_count=payload.get("density_artifact_excluded_count", 0),
        top_fraction=payload.get("top_fraction"),
        random_fraction=payload.get("random_fraction"),
        random_seed=payload.get("random_seed"),
        random_candidate_count=payload.get("random_candidate_count"),
        metadata_domain=payload.get("metadata_domain"),
        metadata_source_file=payload.get("metadata_source_file"),
        metadata_column=payload.get("metadata_column"),
        metadata_values=tuple(payload.get("metadata_values") or ()),
        metadata_source_row_id_field=payload.get("metadata_source_row_id_field"),
        metadata_candidate_count=payload.get("metadata_candidate_count"),
        metadata_missing_count=payload.get("metadata_missing_count"),
        metadata_invalid_count=payload.get("metadata_invalid_count"),
        tie_break_rule=payload.get("tie_break_rule"),
        range_coordinate_source=payload.get("range_coordinate_source"),
        resolved_range_coordinate_source=payload.get("resolved_range_coordinate_source"),
        range_representation=payload.get("range_representation"),
        range_euler_sequence=payload.get("range_euler_sequence"),
        range_euler_convention=payload.get("range_euler_convention"),
        range_scipy_euler_sequence=payload.get("range_scipy_euler_sequence"),
        range_degrees=payload.get("range_degrees"),
        range_bounds=_range_bounds_from_payload(payload.get("range_bounds")),
        selected_particle_keys=tuple(selected_particle_keys),
        selected_count=selected_count,
        total_count=total_count,
        active_policy=active_policy,
    )


def _validate_export_policy(policy: SelectionExportPolicy) -> None:
    if policy.export_reconstructed_maps:
        raise ValueError("Selection export does not support reconstructed map export")
    if policy.invoke_external_tools:
        raise ValueError("Selection export does not invoke external tools")


def _validate_simple_filename(filename: str, field_name: str) -> None:
    path = Path(filename)
    if path.is_absolute() or path.name != filename or ".." in path.parts:
        raise ValueError(f"{field_name} must be a simple filename")
    if filename in {"", ".", ".."}:
        raise ValueError(f"{field_name} must be a simple filename")


def _prepare_output_paths(paths: Iterable[Path], *, overwrite: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(f"Output path already exists: {existing[0]}")


def _write_selected_keys_csv(path: Path, selection: Selection) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["selection_rank", "particle_key"])
        for rank, particle_key in enumerate(selection.selected_particle_keys, start=1):
            writer.writerow([rank, particle_key])


def _write_selection_json(
    path: Path,
    selection: Selection,
    policy: SelectionExportPolicy,
) -> None:
    payload = _selection_json_payload(selection)
    if not policy.include_selected_particle_keys_in_json:
        payload.pop("selected_particle_keys", None)
    payload["active_policy"] = to_json_safe(selection.active_policy)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(to_json_safe(payload), handle, indent=2, sort_keys=True)
        handle.write("\n")


def _write_export_report_json(path: Path, report: ExportReport) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(to_json_safe(report), handle, indent=2, sort_keys=True)
        handle.write("\n")


def _selection_json_payload(selection: Selection) -> dict[str, object]:
    return {
        "artifact_type": "selection",
        "schema_version": "1",
        "selection_id": selection.selection_id,
        "parent_landscape_id": selection.parent_landscape_id,
        "parent_landscape_metadata": selection.parent_landscape_metadata,
        "selection_mode": selection.selection_mode,
        "selection_basis": selection.selection_basis,
        "metric": selection.metric,
        "density_support_field": selection.density_support_field,
        "density_artifact_policy": selection.density_artifact_policy,
        "density_artifact_flag_field": selection.density_artifact_flag_field,
        "density_artifact_candidate_count_before": (
            selection.density_artifact_candidate_count_before
        ),
        "density_artifact_candidate_count_after": (
            selection.density_artifact_candidate_count_after
        ),
        "density_artifact_excluded_count": selection.density_artifact_excluded_count,
        "threshold": selection.threshold,
        "sld_min": selection.sld_min,
        "sld_max": selection.sld_max,
        "threshold_operator": selection.threshold_operator,
        "top_fraction": selection.top_fraction,
        "random_fraction": selection.random_fraction,
        "random_seed": selection.random_seed,
        "random_candidate_count": selection.random_candidate_count,
        "metadata_domain": selection.metadata_domain,
        "metadata_source_file": selection.metadata_source_file,
        "metadata_column": selection.metadata_column,
        "metadata_values": selection.metadata_values,
        "metadata_source_row_id_field": selection.metadata_source_row_id_field,
        "metadata_candidate_count": selection.metadata_candidate_count,
        "metadata_missing_count": selection.metadata_missing_count,
        "metadata_invalid_count": selection.metadata_invalid_count,
        "tie_break_rule": selection.tie_break_rule,
        "center_input": selection.center_input,
        "center_input_representation": selection.center_input_representation,
        "center_input_space": selection.center_input_space,
        "center_euler_sequence": selection.center_euler_sequence,
        "center_euler_convention": selection.center_euler_convention,
        "center_scipy_euler_sequence": selection.center_scipy_euler_sequence,
        "center_degrees": selection.center_degrees,
        "center_evaluated": selection.center_evaluated,
        "evaluation_space": selection.evaluation_space,
        "resolved_evaluation_space": selection.resolved_evaluation_space,
        "radius": selection.radius,
        "radius_unit": selection.radius_unit,
        "range_coordinate_source": selection.range_coordinate_source,
        "resolved_range_coordinate_source": selection.resolved_range_coordinate_source,
        "range_representation": selection.range_representation,
        "range_euler_sequence": selection.range_euler_sequence,
        "range_euler_convention": selection.range_euler_convention,
        "range_scipy_euler_sequence": selection.range_scipy_euler_sequence,
        "range_degrees": selection.range_degrees,
        "range_bounds": selection.range_bounds,
        "selected_count": selection.selected_count,
        "total_count": selection.total_count,
        "selected_particle_keys": selection.selected_particle_keys,
    }


def _validate_selection_artifact_header(payload: dict[str, object]) -> None:
    artifact_type = payload.get("artifact_type")
    if artifact_type is not None and artifact_type != "selection":
        raise ValueError(
            f"Selection JSON artifact_type must be 'selection', got {artifact_type!r}"
        )
    schema_version = payload.get("schema_version")
    if schema_version is not None and schema_version != "1":
        raise ValueError(
            f"Selection JSON schema_version must be '1', got {schema_version!r}"
        )


def _selection_policy_from_payload(payload: dict[str, object]) -> SelectionPolicy:
    data = dict(payload)
    data["center_input"] = _optional_tuple(data.get("center_input"))
    data["range_bounds"] = _range_bounds_from_payload(data.get("range_bounds")) or {}
    data["metadata_values"] = tuple(data.get("metadata_values") or ())
    data["density_artifact_policy"] = _density_artifact_policy_from_payload(
        data.get("density_artifact_policy", "include_all")
    )
    try:
        return SelectionPolicy(**data)
    except TypeError as exc:
        raise ValueError(
            "Selection JSON active_policy is incompatible with SelectionPolicy schema"
        ) from exc


def _non_negative_int(value, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"Selection JSON {field_name} must be a non-negative integer")
    if value < 0:
        raise ValueError(f"Selection JSON {field_name} must be a non-negative integer")
    return value


def _optional_tuple(value):
    if value is None:
        return None
    if not isinstance(value, list):
        raise ValueError("Selection JSON vector fields must be lists when present")
    return tuple(value)


def _range_bounds_from_payload(value):
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError("Selection JSON range_bounds must be an object when present")
    bounds = {}
    for axis_name, axis_bounds in value.items():
        if axis_bounds is None:
            bounds[axis_name] = None
            continue
        if not isinstance(axis_bounds, list) or len(axis_bounds) != 2:
            raise ValueError("Selection JSON range_bounds entries must be length-2 lists")
        bounds[axis_name] = tuple(axis_bounds)
    return bounds


def _density_artifact_policy_from_payload(value):
    if value == "exclude_display_clipped":
        return "exclude_display_outliers"
    return value
