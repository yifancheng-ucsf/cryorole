from __future__ import annotations

import csv
import json
from dataclasses import replace

import pytest

from cryorole.export import export_selection, read_selection_json
from cryorole.models.policies import SelectionExportPolicy, SelectionPolicy
from cryorole.models.selection import Selection


def _make_selection() -> Selection:
    policy = SelectionPolicy(
        selection_mode="threshold_by_density",
        threshold=2.0,
        parent_landscape_id="landscape-1",
        parent_landscape_metadata={"dataset": "synthetic"},
        selection_id="sel-1",
    )
    return Selection(
        selection_id="sel-1",
        parent_landscape_id="landscape-1",
        parent_landscape_metadata={"dataset": "synthetic"},
        selection_mode="threshold_by_density",
        selection_basis="density:sld_raw",
        metric="euclidean",
        threshold=2.0,
        threshold_operator=">=",
        density_support_field="sld_raw",
        top_fraction=None,
        selected_particle_keys=("p3", "p1", "p2"),
        selected_count=3,
        total_count=5,
        active_policy=policy,
    )


def test_selected_particle_keys_csv_preserves_selected_order(tmp_path):
    selection = _make_selection()

    report = export_selection(
        selection,
        policy=SelectionExportPolicy(output_dir=tmp_path),
    )

    with open(report.output_paths["selected_particle_keys_csv"], newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert set(rows[0].keys()) == {"selection_rank", "particle_key"}
    assert [row["selection_rank"] for row in rows] == ["1", "2", "3"]
    assert [row["particle_key"] for row in rows] == ["p3", "p1", "p2"]


def test_selection_json_contains_core_selection_metadata(tmp_path):
    selection = _make_selection()

    report = export_selection(
        selection,
        policy=SelectionExportPolicy(output_dir=tmp_path),
    )

    with open(report.output_paths["selection_json"], encoding="utf-8") as handle:
        payload = json.load(handle)
    assert payload["artifact_type"] == "selection"
    assert payload["schema_version"] == "1"
    assert payload["selection_id"] == "sel-1"
    assert payload["parent_landscape_id"] == "landscape-1"
    assert payload["parent_landscape_metadata"] == {"dataset": "synthetic"}
    assert payload["selection_mode"] == "threshold_by_density"
    assert payload["selection_basis"] == "density:sld_raw"
    assert payload["metric"] == "euclidean"
    assert payload["density_support_field"] == "sld_raw"
    assert payload["threshold"] == 2.0
    assert payload["threshold_operator"] == ">="
    assert payload["selected_count"] == 3
    assert payload["total_count"] == 5
    assert payload["selected_particle_keys"] == ["p3", "p1", "p2"]


def test_active_policy_is_serialized(tmp_path):
    selection = _make_selection()

    report = export_selection(
        selection,
        policy=SelectionExportPolicy(output_dir=tmp_path),
    )

    with open(report.output_paths["selection_json"], encoding="utf-8") as handle:
        payload = json.load(handle)
    assert payload["active_policy"]["selection_mode"] == "threshold_by_density"
    assert payload["active_policy"]["threshold"] == 2.0
    assert payload["active_policy"]["parent_landscape_metadata"] == {"dataset": "synthetic"}


def test_read_selection_json_rehydrates_selection_and_active_policy(tmp_path):
    selection = _make_selection()
    report = export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))

    loaded = read_selection_json(report.output_paths["selection_json"])

    assert loaded.selection_id == "sel-1"
    assert loaded.selected_particle_keys == ("p3", "p1", "p2")
    assert loaded.selected_count == 3
    assert loaded.total_count == 5
    assert loaded.parent_landscape_metadata == {"dataset": "synthetic"}
    assert isinstance(loaded.active_policy, SelectionPolicy)
    assert loaded.active_policy.selection_mode == "threshold_by_density"
    assert loaded.active_policy.threshold == 2.0


def test_read_selection_json_rehydrates_random_policy_fields(tmp_path):
    policy = SelectionPolicy(
        selection_mode="random",
        random_fraction=0.25,
        random_seed=123,
        parent_landscape_metadata={"space": "raw"},
        selection_id="random-sel",
    )
    selection = Selection(
        selection_id="random-sel",
        parent_landscape_id="landscape-1",
        parent_landscape_metadata={"space": "raw"},
        selection_mode="random",
        selection_basis="random_fraction",
        metric="random_without_replacement",
        density_support_field="sld_raw",
        top_fraction=None,
        random_fraction=0.25,
        random_seed=123,
        random_candidate_count=8,
        selected_particle_keys=("p1", "p4"),
        selected_count=2,
        total_count=8,
        active_policy=policy,
    )
    report = export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))

    loaded = read_selection_json(report.output_paths["selection_json"])

    assert loaded.selection_mode == "random"
    assert loaded.random_fraction == 0.25
    assert loaded.random_seed == 123
    assert loaded.random_candidate_count == 8
    assert loaded.active_policy.random_fraction == 0.25
    assert loaded.active_policy.random_seed == 123


def test_selection_json_round_trips_metadata_selection_fields(tmp_path):
    policy = SelectionPolicy(
        selection_mode="metadata_value",
        metadata_domain="ref",
        metadata_column="rlnClassNumber",
        metadata_values=("2", "3"),
        metadata_source_file="ref.star",
        metadata_source_row_id_field="ref_source_row_id",
        selection_id="ref_classes_2_3",
    )
    selection = Selection(
        selection_id="ref_classes_2_3",
        parent_landscape_id="raw",
        selection_mode="metadata_value",
        selection_basis="source_metadata:ref:rlnClassNumber",
        metric="metadata_value_match",
        density_support_field="sld_raw",
        top_fraction=None,
        metadata_domain="ref",
        metadata_source_file="ref.star",
        metadata_column="rlnClassNumber",
        metadata_values=("2", "3"),
        metadata_source_row_id_field="ref_source_row_id",
        metadata_candidate_count=10,
        metadata_missing_count=1,
        metadata_invalid_count=0,
        selected_particle_keys=("p2", "p3"),
        selected_count=2,
        total_count=10,
        active_policy=policy,
    )

    report = export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    loaded = read_selection_json(report.output_paths["selection_json"])

    assert loaded.selection_mode == "metadata_value"
    assert loaded.metric == "metadata_value_match"
    assert loaded.metadata_domain == "ref"
    assert loaded.metadata_column == "rlnClassNumber"
    assert loaded.metadata_values == ("2", "3")
    assert loaded.metadata_source_row_id_field == "ref_source_row_id"
    assert loaded.metadata_candidate_count == 10
    assert loaded.metadata_missing_count == 1
    assert loaded.active_policy.metadata_values == ("2", "3")


def test_read_selection_json_rehydrates_tuple_like_fields(tmp_path):
    selection = Selection(
        selection_id="sel-radius",
        parent_landscape_id="landscape-1",
        selection_mode="radius_around_center",
        selection_basis="radius:analysis",
        metric="euclidean",
        center_input=(1.0, 2.0, 3.0),
        center_input_representation="rotvec",
        center_input_space="evaluation",
        center_evaluated=(1.0, 2.0, 3.0),
        evaluation_space="analysis",
        resolved_evaluation_space="analysis",
        radius=2.0,
        top_fraction=None,
        selected_particle_keys=("p1",),
        selected_count=1,
        total_count=2,
        active_policy=SelectionPolicy(
            selection_mode="radius_around_center",
            center_input=(1.0, 2.0, 3.0),
            evaluation_space="analysis",
            radius=2.0,
        ),
    )
    report = export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))

    loaded = read_selection_json(report.output_paths["selection_json"])

    assert loaded.center_input == (1.0, 2.0, 3.0)
    assert loaded.center_evaluated == (1.0, 2.0, 3.0)
    assert loaded.active_policy.center_input == (1.0, 2.0, 3.0)


def test_read_selection_json_rehydrates_range_bounds(tmp_path):
    selection = Selection(
        selection_id="sel-range",
        parent_landscape_id="landscape-1",
        selection_mode="range_by_coordinates",
        selection_basis="range:euler:analysis",
        metric="coordinate_box",
        range_coordinate_source="analysis",
        resolved_range_coordinate_source="analysis",
        range_representation="euler",
        range_euler_sequence="ZYX",
        range_degrees=True,
        range_bounds={"alpha": (60.0, 120.0), "beta": None},
        top_fraction=None,
        selected_particle_keys=("p2",),
        selected_count=1,
        total_count=3,
        active_policy=SelectionPolicy(
            selection_mode="range_by_coordinates",
            range_coordinate_source="analysis",
            range_bounds={"alpha": (60.0, 120.0), "beta": None},
        ),
    )
    report = export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))

    loaded = read_selection_json(report.output_paths["selection_json"])

    assert loaded.range_bounds == {"alpha": (60.0, 120.0), "beta": None}
    assert loaded.active_policy.range_bounds == {"alpha": (60.0, 120.0), "beta": None}


def test_read_selection_json_supports_legacy_payload_without_header(tmp_path):
    selection = _make_selection()
    report = export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    path = tmp_path / "selection.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload.pop("artifact_type")
    payload.pop("schema_version")
    path.write_text(json.dumps(payload), encoding="utf-8")

    loaded = read_selection_json(report.output_paths["selection_json"])

    assert loaded.selection_id == "sel-1"
    assert loaded.selected_particle_keys == ("p3", "p1", "p2")


def test_read_selection_json_missing_file_fails_clearly(tmp_path):
    with pytest.raises(ValueError, match="does not exist"):
        read_selection_json(tmp_path / "missing_selection.json")


def test_read_selection_json_malformed_json_fails_clearly(tmp_path):
    path = tmp_path / "selection.json"
    path.write_text("{not json", encoding="utf-8")

    with pytest.raises(ValueError, match="Malformed selection JSON"):
        read_selection_json(path)


def test_read_selection_json_requires_object(tmp_path):
    path = tmp_path / "selection.json"
    path.write_text("[]", encoding="utf-8")

    with pytest.raises(ValueError, match="must be an object"):
        read_selection_json(path)


def test_read_selection_json_requires_selected_particle_keys(tmp_path):
    selection = _make_selection()
    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    path = tmp_path / "selection.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload.pop("selected_particle_keys")
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="missing required fields"):
        read_selection_json(path)


def test_read_selection_json_requires_selected_particle_keys_list(tmp_path):
    selection = _make_selection()
    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    path = tmp_path / "selection.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["selected_particle_keys"] = "p1"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="selected_particle_keys must be a list"):
        read_selection_json(path)


def test_read_selection_json_selected_count_mismatch_fails_clearly(tmp_path):
    selection = _make_selection()
    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    path = tmp_path / "selection.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["selected_count"] = 2
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="length must equal selected_count"):
        read_selection_json(path)


def test_read_selection_json_selected_count_greater_than_total_count_fails(tmp_path):
    selection = _make_selection()
    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    path = tmp_path / "selection.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["total_count"] = 2
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="selected_count must be <= total_count"):
        read_selection_json(path)


@pytest.mark.parametrize("field_name", ["selected_count", "total_count"])
def test_read_selection_json_negative_counts_fail_clearly(tmp_path, field_name):
    selection = _make_selection()
    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    path = tmp_path / "selection.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload[field_name] = -1
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match=f"{field_name} must be a non-negative integer"):
        read_selection_json(path)


@pytest.mark.parametrize("field_name", ["selected_count", "total_count"])
def test_read_selection_json_non_integer_counts_fail_clearly(tmp_path, field_name):
    selection = _make_selection()
    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    path = tmp_path / "selection.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload[field_name] = 1.5
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match=f"{field_name} must be a non-negative integer"):
        read_selection_json(path)


def test_read_selection_json_parent_landscape_metadata_must_be_object(tmp_path):
    selection = _make_selection()
    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    path = tmp_path / "selection.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["parent_landscape_metadata"] = ["not", "an", "object"]
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="parent_landscape_metadata must be an object"):
        read_selection_json(path)


def test_read_selection_json_parent_landscape_id_must_be_string_or_null(tmp_path):
    selection = _make_selection()
    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    path = tmp_path / "selection.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["parent_landscape_id"] = 123
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="parent_landscape_id must be a string or null"):
        read_selection_json(path)


def test_read_selection_json_incompatible_active_policy_fails_clearly(tmp_path):
    selection = _make_selection()
    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))
    path = tmp_path / "selection.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["active_policy"]["unexpected_field"] = "surprise"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="active_policy is incompatible"):
        read_selection_json(path)


def test_export_report_records_output_paths_and_counts(tmp_path):
    selection = _make_selection()
    policy = SelectionExportPolicy(output_dir=tmp_path)

    report = export_selection(selection, policy=policy)

    assert report.export_type == "selection_export"
    assert report.schema_version == "1"
    assert set(report.output_paths) == {
        "selected_particle_keys_csv",
        "selection_json",
        "export_report_json",
    }
    assert report.row_counts == {
        "selected_particle_keys_csv": 3,
        "selection_json": 1,
        "export_report_json": 1,
    }
    assert report.source_selection_id == "sel-1"
    assert report.timestamp
    assert report.export_policy["output_dir"] == str(tmp_path)
    assert report.export_policy["overwrite"] is False


def test_export_report_json_is_written(tmp_path):
    selection = _make_selection()

    report = export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))

    with open(report.output_paths["export_report_json"], encoding="utf-8") as handle:
        payload = json.load(handle)
    assert payload["export_type"] == "selection_export"
    assert payload["schema_version"] == "1"
    assert payload["source_selection_id"] == "sel-1"
    assert payload["output_paths"]["selected_particle_keys_csv"].endswith(
        "selected_particle_keys.csv"
    )
    assert payload["row_counts"]["selected_particle_keys_csv"] == 3


def test_export_policy_is_json_safe_in_export_report_json(tmp_path):
    selection = _make_selection()

    report = export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))

    with open(report.output_paths["export_report_json"], encoding="utf-8") as handle:
        payload = json.load(handle)
    assert payload["export_policy"]["output_dir"] == str(tmp_path)
    assert payload["export_policy"]["selected_keys_filename"] == "selected_particle_keys.csv"
    assert payload["export_policy"]["export_report_filename"] == "export_report.json"


def test_existing_output_path_fails_by_default(tmp_path):
    selection = _make_selection()
    policy = SelectionExportPolicy(output_dir=tmp_path)
    export_selection(selection, policy=policy)

    with pytest.raises(FileExistsError, match="Output path already exists"):
        export_selection(selection, policy=policy)


def test_overwrite_true_permits_overwrite(tmp_path):
    selection = _make_selection()
    policy = SelectionExportPolicy(output_dir=tmp_path)
    export_selection(selection, policy=policy)
    smaller_selection = replace(
        selection,
        selected_particle_keys=("p9",),
        selected_count=1,
    )

    export_selection(
        smaller_selection,
        policy=SelectionExportPolicy(output_dir=tmp_path, overwrite=True),
    )

    with open(tmp_path / "selected_particle_keys.csv", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["particle_key"] for row in rows] == ["p9"]


def test_export_does_not_modify_selection_object(tmp_path):
    selection = _make_selection()
    before = selection

    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))

    assert selection == before


def test_export_does_not_write_source_metadata_or_map_files(tmp_path):
    selection = _make_selection()

    export_selection(selection, policy=SelectionExportPolicy(output_dir=tmp_path))

    output_names = {path.name for path in tmp_path.iterdir()}
    assert output_names == {
        "selected_particle_keys.csv",
        "selection.json",
        "export_report.json",
    }


def test_export_rejects_reconstruction_or_external_tool_flags(tmp_path):
    selection = _make_selection()

    with pytest.raises(ValueError, match="reconstructed map"):
        export_selection(
            selection,
            policy=SelectionExportPolicy(output_dir=tmp_path, export_reconstructed_maps=True),
        )
    with pytest.raises(ValueError, match="external tools"):
        export_selection(
            selection,
            policy=SelectionExportPolicy(output_dir=tmp_path, invoke_external_tools=True),
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"selected_keys_filename": "../selected.csv"},
        {"selection_json_filename": "../selection.json"},
        {"export_report_filename": "../report.json"},
    ],
)
def test_parent_traversal_filenames_are_rejected(tmp_path, kwargs):
    selection = _make_selection()

    with pytest.raises(ValueError, match="simple filename"):
        export_selection(
            selection,
            policy=SelectionExportPolicy(output_dir=tmp_path, **kwargs),
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"selected_keys_filename": "C:/tmp/selected.csv"},
        {"selection_json_filename": "C:/tmp/selection.json"},
        {"export_report_filename": "C:/tmp/report.json"},
    ],
)
def test_absolute_filenames_are_rejected(tmp_path, kwargs):
    selection = _make_selection()

    with pytest.raises(ValueError, match="simple filename"):
        export_selection(
            selection,
            policy=SelectionExportPolicy(output_dir=tmp_path, **kwargs),
        )
