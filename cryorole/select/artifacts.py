"""Standard writer for scientific Selection artifact directories."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from cryorole.export.selection_export import export_selection
from cryorole.export.serialization import to_json_safe
from cryorole.export.landscape import write_json_artifact
from cryorole.models.policies import SelectionExportPolicy
from cryorole.models.selection import Selection


@dataclass(frozen=True)
class SelectedRowProvenance:
    """Selected keys and source-row backtracking in stable selection order."""

    particle_key: Sequence[object]
    ref_source_row_id: Sequence[object] | None = None
    mov_source_row_id: Sequence[object] | None = None

    def __post_init__(self) -> None:
        count = len(self.particle_key)
        for values in (self.ref_source_row_id, self.mov_source_row_id):
            if values is not None and len(values) != count:
                raise ValueError("Selected-row provenance arrays must have equal lengths")


@dataclass(frozen=True)
class SelectionArtifactRequest:
    """Inputs required to persist one standard Selection artifact."""

    selection: Selection
    output_dir: str | Path
    selected_rows: SelectedRowProvenance
    overwrite: bool = False
    summary: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class SelectionArtifactResult:
    """Paths and reports written for one standard Selection artifact."""

    output_dir: Path
    selection_json: Path
    selected_particle_keys_csv: Path
    selection_csv: Path
    selected_landscape_rows_csv: Path
    selection_summary_json: Path
    export_report: object


def write_selection_artifact(request: SelectionArtifactRequest) -> SelectionArtifactResult:
    """Write the common CLI/interactive Selection artifact contract."""

    selection = request.selection
    rows = request.selected_rows
    if tuple(str(value) for value in rows.particle_key) != tuple(
        str(value) for value in selection.selected_particle_keys
    ):
        raise ValueError("Selected-row provenance order must match Selection particle keys")
    output_dir = Path(request.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    export_report = export_selection(
        selection,
        policy=SelectionExportPolicy(output_dir=output_dir, overwrite=request.overwrite),
    )
    selection_csv = _write_ranked_keys(
        output_dir / "selection.csv",
        selection.selected_particle_keys,
        overwrite=request.overwrite,
    )
    selected_rows_csv = _write_selected_rows(
        output_dir / "selected_landscape_rows.csv",
        rows,
        overwrite=request.overwrite,
    )
    summary_payload = {
        "artifact_type": "selection_summary",
        "schema_version": "1",
        "timestamp": selection.created_at or datetime.now(timezone.utc).isoformat(),
        "selection_id": selection.selection_id,
        "parent_run_id": selection.parent_run_id,
        "parent_landscape_metadata": to_json_safe(selection.parent_landscape_metadata),
        "selection_mode": selection.selection_mode,
        "selection_basis": selection.selection_basis,
        "metric": selection.metric,
        "selected_count": selection.selected_count,
        "total_count": selection.total_count,
        "interaction_provenance": to_json_safe(selection.interaction_provenance),
        "source_row_provenance_fields": [
            name
            for name, values in (
                ("particle_key", rows.particle_key),
                ("ref_source_row_id", rows.ref_source_row_id),
                ("mov_source_row_id", rows.mov_source_row_id),
            )
            if values is not None
        ],
        "output_paths": to_json_safe(export_report.output_paths),
        "export_report_output_paths": to_json_safe(export_report.output_paths),
        "export_report_row_counts": to_json_safe(export_report.row_counts),
        "export_report": to_json_safe(export_report),
        "selection_csv": str(selection_csv),
        "selected_particle_keys_path": to_json_safe(
            export_report.output_paths.get("selected_particle_keys_csv")
        ),
        "selected_landscape_rows_csv": str(selected_rows_csv),
        **dict(request.summary),
    }
    summary_path = write_json_artifact(
        summary_payload,
        output_dir / "selection_summary.json",
        overwrite=request.overwrite,
    )
    return SelectionArtifactResult(
        output_dir=output_dir,
        selection_json=output_dir / "selection.json",
        selected_particle_keys_csv=output_dir / "selected_particle_keys.csv",
        selection_csv=selection_csv,
        selected_landscape_rows_csv=selected_rows_csv,
        selection_summary_json=Path(summary_path),
        export_report=export_report,
    )


def _write_ranked_keys(path: Path, keys: Sequence[object], *, overwrite: bool) -> Path:
    _prepare_path(path, overwrite=overwrite)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["selection_rank", "particle_key"])
        for rank, particle_key in enumerate(keys, start=1):
            writer.writerow([rank, particle_key])
    return path


def _write_selected_rows(
    path: Path,
    rows: SelectedRowProvenance,
    *,
    overwrite: bool,
) -> Path:
    _prepare_path(path, overwrite=overwrite)
    columns = ["particle_key"]
    if rows.ref_source_row_id is not None:
        columns.append("ref_source_row_id")
    if rows.mov_source_row_id is not None:
        columns.append("mov_source_row_id")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for index, particle_key in enumerate(rows.particle_key):
            row = [particle_key]
            if rows.ref_source_row_id is not None:
                row.append(rows.ref_source_row_id[index])
            if rows.mov_source_row_id is not None:
                row.append(rows.mov_source_row_id[index])
            writer.writerow(row)
    return path


def _prepare_path(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
