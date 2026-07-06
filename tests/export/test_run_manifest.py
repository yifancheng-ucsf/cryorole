from __future__ import annotations

import json
from dataclasses import replace

import pytest

from cryorole.export import write_run_manifest
from cryorole.models.export_report import ExportReport
from cryorole.models.policies import (
    DensityPolicy,
    RunManifestPolicy,
    SelectionExportPolicy,
    SelectionPolicy,
)
from cryorole.models.selection import Selection


def _make_selection() -> Selection:
    policy = SelectionPolicy(
        selection_mode="threshold_by_density",
        threshold=2.0,
        selection_id="sel-1",
    )
    return Selection(
        selection_id="sel-1",
        parent_landscape_id="landscape-1",
        selection_mode="threshold_by_density",
        selection_basis="density:sld_raw",
        metric="euclidean",
        threshold=2.0,
        threshold_operator=">=",
        top_fraction=None,
        selected_particle_keys=("p1", "p2"),
        selected_count=2,
        total_count=3,
        active_policy=policy,
    )


def _make_export_report(tmp_path) -> ExportReport:
    export_policy = SelectionExportPolicy(output_dir=tmp_path)
    return ExportReport(
        output_paths={
            "selected_particle_keys_csv": str(tmp_path / "selected_particle_keys.csv"),
            "selection_json": str(tmp_path / "selection.json"),
            "export_report_json": str(tmp_path / "export_report.json"),
        },
        row_counts={
            "selected_particle_keys_csv": 2,
            "selection_json": 1,
            "export_report_json": 1,
        },
        source_selection_id="sel-1",
        timestamp="2026-04-26T00:00:00+00:00",
        export_policy={
            "output_dir": str(tmp_path),
            "overwrite": False,
        },
    )


def test_manifest_json_is_written(tmp_path):
    output_path = tmp_path / "manifest.json"

    report = write_run_manifest(policy=RunManifestPolicy(output_path=output_path))

    assert output_path.exists()
    assert report.output_path == str(output_path)
    assert report.size_bytes == output_path.stat().st_size
    assert report.size_bytes > 0


def test_manifest_contains_schema_version_and_timestamp(tmp_path):
    output_path = tmp_path / "manifest.json"

    write_run_manifest(policy=RunManifestPolicy(output_path=output_path))

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["manifest_type"] == "cryorole_run_manifest"
    assert payload["schema_version"] == "1"
    assert payload["timestamp"]
    assert payload["cryorole_version"]


def test_manifest_serializes_policies_reports_and_results_json_safe(tmp_path):
    input_file = tmp_path / "input.star"
    input_file.write_text("data_particles\n", encoding="utf-8")
    selection = _make_selection()
    density_policy = DensityPolicy(k_neighbors=3)

    write_run_manifest(
        policy=RunManifestPolicy(
            output_path=tmp_path / "manifest.json",
            compute_file_hashes=True,
        ),
        input_paths=(input_file,),
        source_types={str(input_file): "relion"},
        row_counts={str(input_file): 2},
        active_policies={"density_policy": density_policy},
        density_report={"n_points": 2},
        selection=selection,
    )

    payload = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert payload["input_provenance"]["inputs"][0]["path"] == str(input_file)
    assert payload["input_provenance"]["inputs"][0]["hash_algorithm"] == "sha256"
    assert payload["input_provenance"]["inputs"][0]["hash"]
    assert payload["active_policies"]["density_policy"]["k_neighbors"] == 3
    assert payload["reports"]["density_report"]["n_points"] == 2
    assert payload["results"]["selection_summary"]["selection_id"] == "sel-1"
    assert payload["results"]["selection_summary"]["selected_count"] == 2
    assert payload["results"]["selection_summary"]["total_count"] == 3
    assert payload["results"]["selection_summary"]["selected_particle_keys_included"] is False
    assert "selected_particle_keys" not in payload["results"]["selection_summary"]
    assert payload["results"]["selection_summary"]["active_policy"]["threshold"] == 2.0


def test_selected_particle_keys_are_omitted_from_manifest_by_default(tmp_path):
    write_run_manifest(
        policy=RunManifestPolicy(output_path=tmp_path / "manifest.json"),
        selection=_make_selection(),
    )

    payload = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    summary = payload["results"]["selection_summary"]
    assert summary["selected_particle_keys_included"] is False
    assert summary["selected_count"] == 2
    assert summary["total_count"] == 3
    assert "selected_particle_keys" not in summary


def test_include_selected_particle_keys_true_includes_keys(tmp_path):
    write_run_manifest(
        policy=RunManifestPolicy(
            output_path=tmp_path / "manifest.json",
            include_selected_particle_keys=True,
        ),
        selection=_make_selection(),
    )

    payload = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    summary = payload["results"]["selection_summary"]
    assert summary["selected_particle_keys_included"] is True
    assert summary["selected_particle_keys"] == ["p1", "p2"]


def test_manifest_includes_export_report_output_paths(tmp_path):
    export_report = _make_export_report(tmp_path)

    write_run_manifest(
        policy=RunManifestPolicy(output_path=tmp_path / "manifest.json"),
        export_report=export_report,
    )

    payload = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert payload["reports"]["export_report"]["source_selection_id"] == "sel-1"
    assert payload["output_artifacts"]["export_report_output_paths"] == export_report.output_paths
    assert payload["output_artifacts"]["export_report_row_counts"] == export_report.row_counts


def test_manifest_supports_partial_runs_with_missing_optional_sections(tmp_path):
    write_run_manifest(policy=RunManifestPolicy(output_path=tmp_path / "manifest.json"))

    payload = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert payload["input_provenance"]["inputs"] == []
    assert payload["active_policies"] == {}
    assert payload["reports"]["match_report"] is None
    assert payload["reports"]["export_report"] is None
    assert payload["results"]["selection_summary"] is None
    assert payload["output_artifacts"] == {}


def test_existing_manifest_path_fails_by_default(tmp_path):
    output_path = tmp_path / "manifest.json"
    write_run_manifest(policy=RunManifestPolicy(output_path=output_path))

    with pytest.raises(FileExistsError, match="Manifest output path already exists"):
        write_run_manifest(policy=RunManifestPolicy(output_path=output_path))


def test_invalid_hash_algorithm_fails_even_without_input_paths(tmp_path):
    with pytest.raises(ValueError, match="Unsupported manifest hash_algorithm"):
        write_run_manifest(
            policy=RunManifestPolicy(
                output_path=tmp_path / "manifest.json",
                compute_file_hashes=True,
                hash_algorithm="not-a-real-hash",
            )
        )


def test_manifest_output_path_as_directory_fails_clearly(tmp_path):
    with pytest.raises(IsADirectoryError, match="Manifest output path is a directory"):
        write_run_manifest(policy=RunManifestPolicy(output_path=tmp_path))


def test_non_json_manifest_output_path_fails_clearly(tmp_path):
    with pytest.raises(ValueError, match="must have a .json suffix"):
        write_run_manifest(policy=RunManifestPolicy(output_path=tmp_path / "manifest.txt"))


def test_overwrite_true_permits_overwrite(tmp_path):
    output_path = tmp_path / "manifest.json"
    write_run_manifest(
        policy=RunManifestPolicy(output_path=output_path, workflow_name="first")
    )

    write_run_manifest(
        policy=RunManifestPolicy(
            output_path=output_path,
            workflow_name="second",
            overwrite=True,
        )
    )

    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["workflow_name"] == "second"


def test_manifest_writing_does_not_mutate_input_objects(tmp_path):
    selection = _make_selection()
    export_report = _make_export_report(tmp_path)
    manifest_policy = RunManifestPolicy(output_path=tmp_path / "manifest.json")
    before_selection = replace(selection)
    before_export_report = replace(export_report)

    write_run_manifest(
        policy=manifest_policy,
        selection=selection,
        export_report=export_report,
    )

    assert selection == before_selection
    assert export_report == before_export_report


def test_manifest_does_not_write_source_metadata_map_or_reconstruction_files(tmp_path):
    write_run_manifest(policy=RunManifestPolicy(output_path=tmp_path / "manifest.json"))

    output_names = {path.name for path in tmp_path.iterdir()}
    assert output_names == {"manifest.json"}
