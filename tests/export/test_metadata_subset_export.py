from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from cryorole.cli.main import build_parser, export_metadata_command, export_selection_command
from cryorole.export import export_selection, export_selection_metadata_subset
from cryorole.io.readers.star_reader import read_relion_star
from cryorole.io.writers.landscape_store import write_landscape_npz
from cryorole.models.landscape import Landscape
from cryorole.models.policies import SelectionExportPolicy, SelectionMetadataExportPolicy
from cryorole.models.selection import Selection


def _make_landscape() -> Landscape:
    n = 5
    coordinates = [np.array([float(index), 0.0, 0.0]) for index in range(n)]
    return Landscape(
        data=pd.DataFrame(
            {
                "particle_key": [f"p{index}" for index in range(n)],
                "coordinates_analysis": coordinates,
                "coordinates_display": coordinates,
                "sld_unfloored": np.ones(n),
                "sld_raw": np.ones(n),
                "sld_display": np.ones(n),
                "sld_was_floored": np.zeros(n, dtype=bool),
                "sld_local_k_mean": np.ones(n),
                "sld_effective_local_k_mean": np.ones(n),
                "sld_distance_floor": np.full(n, 1e-4),
                "ref_source_row_id": np.arange(n),
                "mov_source_row_id": np.arange(n - 1, -1, -1),
            }
        )
    )


def _make_selection() -> Selection:
    return Selection(
        selection_id="sel",
        parent_landscape_id="raw",
        selection_mode="threshold_by_density",
        selection_basis="density:sld_raw",
        metric="so3_geodesic",
        density_support_field="sld_raw",
        top_fraction=None,
        threshold=1.0,
        threshold_operator=">=",
        selected_particle_keys=("p3", "p1"),
        selected_count=2,
        total_count=5,
    )


def _write_star(path) -> None:
    rows = "\n".join(
        f"{index + 1}@stack.mrcs {10 + index} {20 + index} {30 + index}"
        for index in range(5)
    )
    path.write_text(
        "\n".join(
            [
                "data_optics",
                "",
                "loop_",
                "_rlnOpticsGroup #1",
                "_rlnOpticsGroupName #2",
                "1 opticsGroup1",
                "",
                "data_particles",
                "",
                "loop_",
                "_rlnImageName #1",
                "_rlnAngleRot #2",
                "_rlnAngleTilt #3",
                "_rlnAnglePsi #4",
                rows,
                "",
            ]
        ),
        encoding="utf-8",
    )


def _write_cs(path) -> np.ndarray:
    dtype = np.dtype(
        [
            ("uid", "<u8"),
            ("alignments3D/pose", "<f4", (3,)),
            ("alignments3D/shift", "<f4", (2,)),
        ]
    )
    array = np.zeros(5, dtype=dtype)
    array["uid"] = np.array([100, 101, 102, 103, 104], dtype=np.uint64)
    array["alignments3D/pose"] = np.arange(15, dtype=np.float32).reshape(5, 3)
    array["alignments3D/shift"] = np.arange(10, dtype=np.float32).reshape(5, 2)
    with path.open("wb") as handle:
        np.save(handle, array)
    return array


def _make_run_bundle(tmp_path, *, ref_suffix=".star", mov_suffix=".star"):
    run_dir = tmp_path / "run"
    data_dir = run_dir / "data"
    selection_dir = run_dir / "selections" / "sel"
    data_dir.mkdir(parents=True)
    ref_path = tmp_path / f"ref{ref_suffix}"
    mov_path = tmp_path / f"mov{mov_suffix}"
    if ref_suffix == ".star":
        _write_star(ref_path)
    else:
        _write_cs(ref_path)
    if mov_suffix == ".star":
        _write_star(mov_path)
    else:
        _write_cs(mov_path)
    write_landscape_npz(_make_landscape(), data_dir / "raw_landscape.npz")
    (run_dir / "run_summary.json").write_text(
        json.dumps(
            {
                "input_paths": {"ref": str(ref_path), "mov": str(mov_path)},
                "source_types": {
                    "ref": "relion" if ref_suffix == ".star" else "cryosparc",
                    "mov": "relion" if mov_suffix == ".star" else "cryosparc",
                },
            }
        ),
        encoding="utf-8",
    )
    report = export_selection(
        _make_selection(),
        policy=SelectionExportPolicy(output_dir=selection_dir),
    )
    return run_dir, ref_path, mov_path, report.output_paths["selection_json"]


def test_direct_cli_parser_accepts_metadata_export_options(tmp_path) -> None:
    parser = build_parser()

    args = parser.parse_args(
        [
            "export",
            "--run-dir",
            str(tmp_path / "run"),
            "--selection-id",
            "sel",
            "--domain",
            "both",
            "--format",
            "auto",
        ]
    )

    assert args.command == "export"
    assert args.handler is export_metadata_command
    assert args.domain == "both"
    assert args.format == "auto"


def test_export_help_documents_public_and_advanced_inputs(capsys) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["export", "--help"])

    help_text = capsys.readouterr().out
    for option in (
        "--run-dir",
        "--selection-id",
        "--domain",
        "--format",
        "--overwrite",
        "--selection",
        "--output-dir",
    ):
        assert option in help_text
    assert "primary export" in help_text
    assert "input with --selection-id" in help_text
    assert "Advanced: direct path to a selection.json artifact" in help_text
    assert "RUN/exports/<selection_id>/" in help_text
    assert "source, run, and selection files are unchanged" in help_text


def test_export_selection_alias_accepts_run_dir_selection_id(tmp_path) -> None:
    parser = build_parser()

    args = parser.parse_args(
        [
            "export",
            "selection",
            "--run-dir",
            str(tmp_path / "run"),
            "--selection-id",
            "sel",
        ]
    )

    assert args.handler is export_selection_command
    assert args.run_dir == str(tmp_path / "run")
    assert args.selection_id == "sel"


@pytest.mark.parametrize("domain", ["ref", "mov"])
def test_direct_cli_parser_accepts_domain_variants(tmp_path, domain) -> None:
    parser = build_parser()

    args = parser.parse_args(
        [
            "export",
            "--run-dir",
            str(tmp_path / "run"),
            "--selection-id",
            "sel",
            "--domain",
            domain,
        ]
    )

    assert args.domain == domain


def test_direct_cli_command_exports_from_run_dir_selection_id(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "sel",
            "--domain",
            "ref",
            "--format",
            "auto",
        ]
    )

    assert export_metadata_command(args) == 0

    assert (run_dir / "exports" / "sel" / "selected_particle_keys.txt").exists()
    assert (run_dir / "exports" / "sel" / "selected_landscape_rows.csv").exists()
    assert (run_dir / "exports" / "sel" / "ref" / "selected_ref.star").exists()


def test_export_output_dir_overrides_default_as_advanced_path(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)
    output_dir = tmp_path / "custom_export"
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "sel",
            "--domain",
            "ref",
            "--output-dir",
            str(output_dir),
        ]
    )

    assert export_metadata_command(args) == 0

    assert (output_dir / "selected_particle_keys.txt").exists()
    assert (output_dir / "ref" / "selected_ref.star").exists()
    assert not (run_dir / "exports" / "sel").exists()


def test_export_metadata_command_does_not_reselect(tmp_path, monkeypatch) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)
    monkeypatch.setattr(
        "cryorole.cli.main.select_particles",
        lambda *_args, **_kwargs: pytest.fail("export must not recompute selection"),
    )
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "sel",
            "--domain",
            "ref",
        ]
    )

    assert export_metadata_command(args) == 0


def test_export_missing_selection_input_errors_are_actionable(tmp_path) -> None:
    parser = build_parser()

    missing_selection_id = parser.parse_args(["export", "--run-dir", str(tmp_path / "run")])
    with pytest.raises(ValueError, match="Provide --run-dir RUN --selection-id ID"):
        export_metadata_command(missing_selection_id)

    missing_run_dir = parser.parse_args(["export", "--selection-id", "sel"])
    with pytest.raises(ValueError, match="Provide --run-dir RUN --selection-id ID"):
        export_metadata_command(missing_run_dir)


def test_export_selection_and_selection_id_conflict(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, selection_path = _make_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "sel",
            "--selection",
            selection_path,
        ]
    )

    with pytest.raises(ValueError, match="Use either --selection or --selection-id"):
        export_metadata_command(args)


def test_export_missing_selection_file_error_mentions_compact_and_advanced_inputs(tmp_path) -> None:
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "--run-dir",
            str(tmp_path / "run"),
            "--selection-id",
            "missing",
        ]
    )

    with pytest.raises(ValueError, match="--selection PATH/to/selection.json"):
        export_metadata_command(args)


def test_direct_selection_artifact_export_requires_output_dir(tmp_path) -> None:
    _run_dir, _ref_path, _mov_path, selection_path = _make_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(["export", "--selection", selection_path])

    with pytest.raises(ValueError, match="Direct --selection export requires --output-dir"):
        export_metadata_command(args)


def test_direct_selection_artifact_export_works_with_output_dir(tmp_path) -> None:
    _run_dir, _ref_path, _mov_path, selection_path = _make_run_bundle(tmp_path)
    output_dir = tmp_path / "selection_reexport"
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "--selection",
            selection_path,
            "--output-dir",
            str(output_dir),
        ]
    )

    assert export_metadata_command(args) == 0

    assert (output_dir / "selection.json").exists()
    assert (output_dir / "selected_particle_keys.csv").exists()
    assert (output_dir / "export_report.json").exists()


def test_relion_star_subset_export_preserves_columns_and_source(tmp_path) -> None:
    run_dir, ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)
    source_before = ref_path.read_text(encoding="utf-8")

    report = export_selection_metadata_subset(
        _make_selection(),
        policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="ref", format="auto"),
        selection_path=run_dir / "selections" / "sel" / "selection.json",
    )

    output_path = run_dir / "exports" / "sel" / "ref" / "selected_ref.star"
    assert output_path.exists()
    particles = read_relion_star(output_path).particles
    assert len(particles) == 2
    assert particles.columns.tolist() == [
        "_rlnImageName",
        "_rlnAngleRot",
        "_rlnAngleTilt",
        "_rlnAnglePsi",
    ]
    assert particles["_rlnImageName"].tolist() == ["4@stack.mrcs", "2@stack.mrcs"]
    assert "data_optics" in output_path.read_text(encoding="utf-8")
    assert ref_path.read_text(encoding="utf-8") == source_before
    assert report["outputs"]["ref_metadata"] == str(output_path)


def test_cryosparc_cs_subset_export_preserves_dtype_uid_and_source(tmp_path) -> None:
    run_dir, _ref_path, mov_path, _selection_path = _make_run_bundle(
        tmp_path,
        mov_suffix=".cs",
    )
    source_before = mov_path.read_bytes()
    source_array = np.load(mov_path, allow_pickle=False)

    export_selection_metadata_subset(
        _make_selection(),
        policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="mov", format="auto"),
        selection_path=run_dir / "selections" / "sel" / "selection.json",
    )

    output_path = run_dir / "exports" / "sel" / "mov" / "selected_mov.cs"
    uid_path = run_dir / "exports" / "sel" / "mov" / "selected_mov_uids.txt"
    subset = np.load(output_path, allow_pickle=False)
    assert subset.dtype == source_array.dtype
    np.testing.assert_array_equal(subset["uid"], source_array["uid"][[1, 3]])
    np.testing.assert_allclose(
        subset["alignments3D/pose"],
        source_array["alignments3D/pose"][[1, 3]],
    )
    assert uid_path.read_text(encoding="utf-8").splitlines() == ["101", "103"]
    assert mov_path.read_bytes() == source_before


def test_both_domains_export_uses_ref_and_mov_row_ids(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)

    export_selection_metadata_subset(
        _make_selection(),
        policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="both", format="auto"),
        selection_path=run_dir / "selections" / "sel" / "selection.json",
    )

    ref_particles = read_relion_star(
        run_dir / "exports" / "sel" / "ref" / "selected_ref.star"
    ).particles
    mov_particles = read_relion_star(
        run_dir / "exports" / "sel" / "mov" / "selected_mov.star"
    ).particles
    assert ref_particles["_rlnImageName"].tolist() == ["4@stack.mrcs", "2@stack.mrcs"]
    assert mov_particles["_rlnImageName"].tolist() == ["2@stack.mrcs", "4@stack.mrcs"]


def test_mixed_format_auto_exports_star_and_cs(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(
        tmp_path,
        ref_suffix=".star",
        mov_suffix=".cs",
    )

    report = export_selection_metadata_subset(
        _make_selection(),
        policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="both", format="auto"),
        selection_path=run_dir / "selections" / "sel" / "selection.json",
    )

    assert report["format_resolved"] == {"ref": "relion_star", "mov": "cryosparc_cs"}
    assert (run_dir / "exports" / "sel" / "ref" / "selected_ref.star").exists()
    assert (run_dir / "exports" / "sel" / "mov" / "selected_mov.cs").exists()


def test_export_report_records_public_audit_fields(tmp_path) -> None:
    run_dir, ref_path, mov_path, _selection_path = _make_run_bundle(tmp_path)

    report = export_selection_metadata_subset(
        _make_selection(),
        policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="both", format="auto"),
        selection_path=run_dir / "selections" / "sel" / "selection.json",
    )

    report_path = run_dir / "exports" / "sel" / "export_report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["run_dir"] == str(run_dir.resolve())
    assert payload["selection_id"] == "sel"
    assert payload["selection_path"] == str((run_dir / "selections" / "sel" / "selection.json").resolve())
    assert payload["selected_count"] == 2
    assert payload["domain"] == "both"
    assert payload["format_requested"] == "auto"
    assert payload["format_resolved"] == {"ref": "relion_star", "mov": "relion_star"}
    assert payload["output_dir"] == str(run_dir.resolve() / "exports" / "sel")
    assert payload["overwrite"] is False
    assert payload["source_row_id_fields"] == {
        "ref": "ref_source_row_id",
        "mov": "mov_source_row_id",
    }
    assert payload["source_row_id_stats"]["ref"] == {
        "count": 2,
        "min": 1,
        "max": 3,
        "unique_count": 2,
    }
    assert payload["source_row_id_stats"]["mov"] == {
        "count": 2,
        "min": 1,
        "max": 3,
        "unique_count": 2,
    }
    assert payload["source_files"] == {"ref": str(ref_path), "mov": str(mov_path)}
    assert isinstance(payload["warnings"], list)
    assert report["output_dir"] == payload["output_dir"]


def test_keys_format_writes_only_key_and_row_tables(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)

    report = export_selection_metadata_subset(
        _make_selection(),
        policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="both", format="keys"),
        selection_path=run_dir / "selections" / "sel" / "selection.json",
    )

    assert (run_dir / "exports" / "sel" / "selected_particle_keys.txt").exists()
    assert (run_dir / "exports" / "sel" / "selected_landscape_rows.csv").exists()
    assert not (run_dir / "exports" / "sel" / "ref").exists()
    assert report["format_resolved"] == {"ref": "keys", "mov": "keys"}


def test_domain_ref_exports_only_ref(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)

    export_selection_metadata_subset(
        _make_selection(),
        policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="ref", format="auto"),
        selection_path=run_dir / "selections" / "sel" / "selection.json",
    )

    assert (run_dir / "exports" / "sel" / "ref" / "selected_ref.star").exists()
    assert not (run_dir / "exports" / "sel" / "mov").exists()


def test_domain_mov_exports_only_mov(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)

    export_selection_metadata_subset(
        _make_selection(),
        policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="mov", format="auto"),
        selection_path=run_dir / "selections" / "sel" / "selection.json",
    )

    assert (run_dir / "exports" / "sel" / "mov" / "selected_mov.star").exists()
    assert not (run_dir / "exports" / "sel" / "ref").exists()


def test_missing_selected_row_id_column_fails(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)
    rows_path = run_dir / "selections" / "sel" / "selected_landscape_rows.csv"
    pd.DataFrame({"particle_key": ["p3", "p1"], "mov_source_row_id": [1, 3]}).to_csv(
        rows_path,
        index=False,
    )

    with pytest.raises(ValueError, match="ref_source_row_id"):
        export_selection_metadata_subset(
            _make_selection(),
            policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="ref", format="auto"),
            selection_path=run_dir / "selections" / "sel" / "selection.json",
        )


def test_missing_mov_row_id_column_fails(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)
    rows_path = run_dir / "selections" / "sel" / "selected_landscape_rows.csv"
    pd.DataFrame({"particle_key": ["p3", "p1"], "ref_source_row_id": [3, 1]}).to_csv(
        rows_path,
        index=False,
    )

    with pytest.raises(ValueError, match="mov_source_row_id"):
        export_selection_metadata_subset(
            _make_selection(),
            policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="mov", format="auto"),
            selection_path=run_dir / "selections" / "sel" / "selection.json",
        )


def test_out_of_bounds_row_id_fails(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)
    rows_path = run_dir / "selections" / "sel" / "selected_landscape_rows.csv"
    pd.DataFrame(
        {
            "particle_key": ["p3", "p1"],
            "ref_source_row_id": [99, 1],
            "mov_source_row_id": [1, 3],
        }
    ).to_csv(rows_path, index=False)

    with pytest.raises(ValueError, match="out-of-bounds"):
        export_selection_metadata_subset(
            _make_selection(),
            policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="ref", format="auto"),
            selection_path=run_dir / "selections" / "sel" / "selection.json",
        )


def test_forced_format_incompatible_with_source_fails(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(
        tmp_path,
        ref_suffix=".cs",
    )

    with pytest.raises(ValueError, match="incompatible"):
        export_selection_metadata_subset(
            _make_selection(),
            policy=SelectionMetadataExportPolicy(
                run_dir=run_dir,
                domain="ref",
                format="relion_star",
            ),
            selection_path=run_dir / "selections" / "sel" / "selection.json",
        )


def test_existing_output_without_overwrite_fails(tmp_path) -> None:
    run_dir, _ref_path, _mov_path, _selection_path = _make_run_bundle(tmp_path)
    (run_dir / "exports" / "sel").mkdir(parents=True)

    with pytest.raises(FileExistsError, match="Export output directory already exists"):
        export_selection_metadata_subset(
            _make_selection(),
            policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="ref", format="auto"),
            selection_path=run_dir / "selections" / "sel" / "selection.json",
        )


def test_existing_output_with_overwrite_preserves_inputs(tmp_path) -> None:
    run_dir, ref_path, mov_path, _selection_path = _make_run_bundle(tmp_path)
    raw_path = run_dir / "data" / "raw_landscape.npz"
    selection_path = run_dir / "selections" / "sel" / "selection.json"
    raw_before = raw_path.read_bytes()
    selection_before = selection_path.read_text(encoding="utf-8")
    ref_before = ref_path.read_text(encoding="utf-8")
    mov_before = mov_path.read_text(encoding="utf-8")

    export_selection_metadata_subset(
        _make_selection(),
        policy=SelectionMetadataExportPolicy(run_dir=run_dir, domain="ref", format="auto"),
        selection_path=selection_path,
    )
    report = export_selection_metadata_subset(
        _make_selection(),
        policy=SelectionMetadataExportPolicy(
            run_dir=run_dir,
            domain="ref",
            format="auto",
            overwrite=True,
        ),
        selection_path=selection_path,
    )

    assert report["overwrite"] is True
    assert raw_path.read_bytes() == raw_before
    assert selection_path.read_text(encoding="utf-8") == selection_before
    assert ref_path.read_text(encoding="utf-8") == ref_before
    assert mov_path.read_text(encoding="utf-8") == mov_before
