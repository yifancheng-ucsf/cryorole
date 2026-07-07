from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from cryorole.cli.main import (
    build_parser,
    canonicalize_command,
    export_metadata_command,
    export_selection_command,
    main,
    manifest_command,
    run_command,
    select_command,
    visualize_command,
)
from cryorole.export import export_selection, write_landscape_json, write_run_manifest
from cryorole.io.writers.landscape_store import (
    read_landscape,
    write_landscape_npz,
    write_canonical_landscape_csv_from_npz as real_write_canonical_landscape_csv_from_npz,
)
from cryorole.models.canonicalization_report import CanonicalizationReport
from cryorole.models.density_report import DensityReport
from cryorole.models.landscape import Landscape
from cryorole.models.policies import SelectionExportPolicy, SelectionPolicy
from cryorole.models.selection import Selection


class FakeRunner:
    def __init__(self) -> None:
        self.calls = []

    def run_phase1(self, *args, **kwargs):
        self.calls.append(("run_phase1", args, kwargs))
        return SimpleNamespace(
            pose_a=SimpleNamespace(data=pd.DataFrame({"particle_key": ["p1", "p2", "p3"]})),
            pose_b=SimpleNamespace(data=pd.DataFrame({"particle_key": ["p1", "p2", "p3"]})),
            match_report={
                "matched_count": 3,
                "status": "ok",
                "match_key": "uid",
                "ref_row_count": 3,
                "mov_row_count": 3,
                "matched_row_count": 3,
                "dropped_ref_only_count": 0,
                "dropped_mov_only_count": 0,
                "matched_rows_reordered": False,
                "warnings": (),
            },
        )

    def compute_ro_for_phase1_result(self, phase1, **kwargs):
        self.calls.append(("compute_ro_for_phase1_result", phase1, kwargs))
        return "phase2"

    def compute_density_for_ro_result(self, phase2, *, density_policy=None):
        self.calls.append(("compute_density_for_ro_result", phase2, density_policy))
        return SimpleNamespace(landscape=_make_landscape())

    def canonicalize_density_result(self, phase3, *, canonicalization_policy=None):
        self.calls.append(("canonicalize_density_result", phase3, canonicalization_policy))
        return SimpleNamespace(landscape=_make_canonical_landscape())

    def write_manifest(self, *, manifest_policy, **manifest_sections):
        self.calls.append(("write_manifest", manifest_policy, manifest_sections))
        report = write_run_manifest(policy=manifest_policy, **manifest_sections)
        return SimpleNamespace(manifest_report=report)


class HighSldRunner(FakeRunner):
    def compute_density_for_ro_result(self, phase2, *, density_policy=None):
        self.calls.append(("compute_density_for_ro_result", phase2, density_policy))
        return SimpleNamespace(landscape=_make_landscape(sld_values=[1.0, 20.0, 150.0]))


def test_run_help_shows_compact_public_surface(capsys) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["run", "--help"])

    help_text = capsys.readouterr().out
    for public_option in (
        "--ref",
        "--mov",
        "--ref-domain",
        "--mov-domain",
        "--output-dir",
        "--row-aligned",
        "--no-visualize",
    ):
        assert public_option in help_text
    for hidden_option in (
        "--visual-style",
        "--color-vmax",
        "--display-sld-threshold",
        "--run-backend",
        "--identity-mode",
        "--profile-time",
        "--raw-csv-chunk-size",
        "--write-debug-json",
    ):
        assert hidden_option not in help_text


def test_run_cli_rejects_removed_percentile_sld_display_arguments() -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "run",
                "--ref",
                "ref.cs",
                "--mov",
                "mov.cs",
                "--sld-display-mode",
                "robust_percentile_clip",
            ]
        )

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "run",
                "--ref",
                "ref.cs",
                "--mov",
                "mov.cs",
                "--sld-display-clip-percentile",
                "95",
            ]
        )


def _make_density_report() -> DensityReport:
    return DensityReport(
        n_points=3,
        requested_k_neighbors=2,
        effective_k_neighbors=2,
        global_local_k_mean=1.0,
        distance_floor=0.0001,
        distance_floor_fraction=1e-4,
        n_floored_points=0,
        fraction_floored_points=0.0,
        n_inf_sld_unfloored=0,
        n_inf_sld_raw=0,
        max_sld_unfloored=1.0,
        max_sld_raw=1.0,
        p99_sld_unfloored=1.0,
        p99_sld_raw=1.0,
        max_over_p99_unfloored=1.0,
        max_over_p99_raw=1.0,
        floored_particle_keys=(),
        floored_rows=(),
    )


def _make_landscape(sld_values: list[float] | None = None) -> Landscape:
    sld_values = sld_values or [1.0, 1.0, 1.0]
    data = pd.DataFrame(
        {
        "particle_key": ["p1", "p2", "p3"],
            "coordinates_analysis": [
                np.array([0.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0.0, 1.0, 0.0]),
            ],
            "coordinates_display": [
                np.array([0.0, 0.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([0.0, 1.0, 0.0]),
            ],
            "sld_unfloored": sld_values,
            "sld_raw": sld_values,
            "sld_display": sld_values,
            "sld_display_is_outlier": [False, False, False],
            "sld_was_floored": [False, False, False],
            "sld_local_k_mean": [1.0, 1.0, 1.0],
            "sld_effective_local_k_mean": [1.0, 1.0, 1.0],
            "sld_distance_floor": [0.0001, 0.0001, 0.0001],
        }
    )
    return Landscape(
        data=data,
        active_policies={},
        density_report=_make_density_report(),
    )


def _make_source_row_landscape() -> Landscape:
    landscape = _make_landscape()
    data = landscape.data.copy(deep=True)
    data["ref_source_row_id"] = [0, 1, 2]
    data["mov_source_row_id"] = [2, 1, 0]
    return Landscape(
        data=data,
        active_policies=landscape.active_policies,
        density_report=landscape.density_report,
    )


def _write_metadata_star(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "data_optics",
                "",
                "loop_",
                "_rlnOpticsGroup #1",
                "_rlnVoltage #2",
                "1 300",
                "",
                "data_particles",
                "",
                "loop_",
                "_rlnImageName #1",
                "_rlnAngleRot #2",
                "_rlnAngleTilt #3",
                "_rlnAnglePsi #4",
                "_rlnClassNumber #5",
                "_rlnOpticsGroup #6",
                "000001@stack.mrcs 0 0 0 1 10",
                "000002@stack.mrcs 0 0 0 2 20",
                "000003@stack.mrcs 0 0 0 1 10",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _make_metadata_selection_run_bundle(tmp_path):
    run_dir = tmp_path / "run"
    data_dir = run_dir / "data"
    data_dir.mkdir(parents=True)
    ref_path = tmp_path / "ref.star"
    mov_path = tmp_path / "mov.star"
    _write_metadata_star(ref_path)
    _write_metadata_star(mov_path)
    write_landscape_npz(_make_source_row_landscape(), data_dir / "raw_landscape.npz")
    (run_dir / "run_summary.json").write_text(
        json.dumps(
            {
                "input_paths": {"ref": str(ref_path), "mov": str(mov_path)},
                "source_types": {"ref": "relion", "mov": "relion"},
                "euler_convention": "extrinsic_zyx",
                "scipy_euler_sequence": "zyx",
            }
        ),
        encoding="utf-8",
    )
    return run_dir, ref_path, mov_path


def _make_select_run_bundle(tmp_path, landscape: Landscape | None = None) -> Path:
    run_dir = tmp_path / "run"
    data_dir = run_dir / "data"
    data_dir.mkdir(parents=True)
    write_landscape_npz(landscape or _make_landscape(), data_dir / "raw_landscape.npz")
    (run_dir / "run_summary.json").write_text(
        json.dumps(
            {
                "euler_convention": "extrinsic_zyx",
                "scipy_euler_sequence": "zyx",
            }
        ),
        encoding="utf-8",
    )
    return run_dir


def _make_canonical_landscape() -> Landscape:
    landscape = _make_landscape()
    data = landscape.data.copy(deep=True)
    data["coordinates_canonical"] = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
    ]
    return Landscape(
        data=data,
        canonical_transform=np.eye(3),
        active_policies={},
        density_report=landscape.density_report,
        canonicalization_report=CanonicalizationReport(
            fit_subset="all",
            n_fit_points=3,
            density_support_field="sld_raw",
            singular_values=(1.0, 0.5, 0.1),
            explained_variance_ratios=(0.7, 0.2, 0.1),
            axis_degeneracy_detected=False,
            warnings=(),
        ),
    )


EXPECTED_DEFAULT_VISUALIZATION_FILES = (
    "euler_alpha_beta.png",
    "euler_beta_gamma.png",
    "euler_alpha_gamma.png",
    "euler_3view_projection.png",
    "rotvec_xy.png",
    "rotvec_yz.png",
    "rotvec_xz.png",
    "rotvec_3view_projection.png",
    "landscape_3d_euler.png",
    "landscape_3d_rotvec.png",
    "display_table.csv",
    "visualization_report.json",
)


def _assert_default_visualization_group(group_dir: Path) -> dict:
    for filename in EXPECTED_DEFAULT_VISUALIZATION_FILES:
        assert (group_dir / filename).exists()
    return json.loads((group_dir / "visualization_report.json").read_text(encoding="utf-8"))


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


def _make_subset_selection(
    *,
    selection_id: str = "subset",
    selected_particle_keys: tuple[str, ...] = ("p3", "p1"),
    total_count: int = 3,
) -> Selection:
    policy = SelectionPolicy(
        selection_mode="threshold_by_density",
        threshold=1.0,
        parent_landscape_id="landscape-1",
        parent_landscape_metadata={"dataset": "synthetic"},
        selection_id=selection_id,
    )
    return Selection(
        selection_id=selection_id,
        parent_landscape_id="landscape-1",
        parent_landscape_metadata={"dataset": "synthetic"},
        selection_mode="threshold_by_density",
        selection_basis="density:sld_raw",
        metric="euclidean",
        threshold=1.0,
        threshold_operator=">=",
        density_support_field="sld_raw",
        top_fraction=None,
        selected_particle_keys=selected_particle_keys,
        selected_count=len(selected_particle_keys),
        total_count=total_count,
        active_policy=policy,
    )


def test_cli_command_tree_parses_expected_commands() -> None:
    parser = build_parser()

    assert parser.parse_args(["align", "--ref", "ref.star", "--mov", "mov.star"]).command == "align"
    assert parser.parse_args(["visualize", "--run-dir", "run"]).command == "visualize"
    assert (
        parser.parse_args(
            ["visualize", "--run-dir", "run", "--selection-id", "subset"]
        ).selection_id
        == "subset"
    )
    assert parser.parse_args(["canonicalize", "--landscape", "l.json", "--output", "c.json"]).command == "canonicalize"
    assert parser.parse_args(["export", "selection", "--selection", "s.json", "--output-dir", "out"]).command == "export"
    assert parser.parse_args(["manifest", "--output", "manifest.json"]).command == "manifest"


def test_run_command_writes_landscape_density_report_and_manifest(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir)])

    assert run_command(args, runner=FakeRunner()) == 0

    assert (output_dir / "data" / "raw_landscape.npz").exists()
    assert (output_dir / "data" / "raw_landscape.csv").exists()
    assert (output_dir / "data" / "match_table.csv").exists()
    assert (output_dir / "reports" / "density_report.json").exists()
    assert (output_dir / "reports" / "match_report.json").exists()
    assert (output_dir / "run_summary.json").exists()
    assert (output_dir / "visualizations" / "raw_default" / "visualization_report.json").exists()
    assert (output_dir / "run_manifest.json").exists()
    assert not (output_dir / "landscape.json").exists()
    raw = pd.read_csv(output_dir / "data" / "raw_landscape.csv")
    expected_columns = {
        "particle_key",
        "raw_rv_x_rad",
        "raw_ea_zyx_alpha_deg",
        "rotvec_x",
        "rotvec_y",
        "rotvec_z",
        "euler_alpha",
        "euler_beta",
        "euler_gamma",
        "sld_unfloored",
        "sld_raw",
        "sld_display",
        "coordinate_source",
    }
    assert expected_columns.issubset(raw.columns)
    assert len(raw) == 3
    assert not any("intrinsic" in column or "extrinsic" in column for column in raw.columns)
    visualization_report = json.loads(
        (
            output_dir
            / "visualizations"
            / "raw_default"
            / "all_particles"
            / "visualization_report.json"
        ).read_text(encoding="utf-8")
    )
    assert visualization_report["euler_convention"] == "extrinsic_zyx"
    assert visualization_report["scipy_euler_sequence"] == "zyx"


def test_run_debug_landscape_json_is_opt_in(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--write-debug-json",
        ]
    )

    run_command(args, runner=FakeRunner())

    debug_json = output_dir / "debug" / "landscape_debug.json"
    assert debug_json.exists()
    payload = json.loads(debug_json.read_text(encoding="utf-8"))
    assert payload["artifact_type"] == "landscape"


def test_run_manifest_records_resolved_cryosparc_uid_and_core_policies(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir)])

    run_command(args, runner=FakeRunner())

    manifest = json.loads((output_dir / "run_manifest.json").read_text(encoding="utf-8"))
    active_policies = manifest["active_policies"]
    assert active_policies["identity_policy_ref"]["identity_mode"] == "cryosparc_uid"
    assert active_policies["identity_policy_mov"]["identity_mode"] == "cryosparc_uid"
    assert active_policies["identity_policy_ref"]["identity_columns"] == ["uid"]
    assert active_policies["convention_policy_ref"] is None
    assert active_policies["convention_policy_mov"] is None
    assert active_policies["match_policy"]["join_type"] == "inner"
    assert active_policies["density_policy"]["k_neighbors"] == 50
    assert active_policies["representation_policy"]["euler_convention"] == "extrinsic_zyx"
    assert active_policies["representation_policy"]["scipy_euler_sequence"] == "zyx"
    assert "canonicalization_policy" not in active_policies
    assert manifest["results"]["euler_metadata"]["euler_convention"] == "extrinsic_zyx"
    assert manifest["output_artifacts"]["run_summary_json"] == str(output_dir / "run_summary.json")
    assert manifest["schema_version"] == "2.0"
    assert manifest["output_artifacts"]["raw_landscape_npz"] == str(
        output_dir / "data" / "raw_landscape.npz"
    )
    assert manifest["output_artifacts"]["raw_landscape_csv"] == str(
        output_dir / "data" / "raw_landscape.csv"
    )
    assert manifest["output_artifacts"]["visualization_report_json"] == str(
        output_dir / "visualizations" / "raw_default" / "visualization_report.json"
    )


def test_run_summary_contains_counts_and_canonicalization_status(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir), "--k-neighbors", "7"]
    )

    run_command(args, runner=FakeRunner())

    summary = json.loads((output_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["artifact_type"] == "run_summary"
    assert summary["input_paths"] == {"ref": "ref.cs", "mov": "mov.cs"}
    assert summary["source_types"] == {"ref": "cryosparc", "mov": "cryosparc"}
    assert summary["ref_domain"] == "ref"
    assert summary["mov_domain"] == "mov"
    assert summary["row_counts"] == {"ref": 3, "mov": 3}
    assert summary["matched_count"] == 3
    assert summary["match_key"] == "uid"
    assert summary["matched_row_count"] == 3
    assert summary["dropped_ref_only_count"] == 0
    assert summary["dropped_mov_only_count"] == 0
    assert summary["matched_rows_reordered"] is False
    assert summary["match_warnings"] == []
    assert summary["landscape_row_count"] == 3
    assert summary["k_neighbors"] == 7
    assert summary["canonicalization_performed"] is False
    assert summary["selection_performed"] is False
    assert summary["output_artifacts"]["raw_landscape_npz"] == str(
        output_dir / "data" / "raw_landscape.npz"
    )
    assert summary["output_artifacts"]["raw_landscape_csv"] == str(
        output_dir / "data" / "raw_landscape.csv"
    )
    assert summary["run_backend_requested"] == "auto"
    assert summary["run_backend_resolved"] == "dataframe_compat"
    assert summary["raw_csv_performed"] is True
    assert summary["raw_csv_backend"] == "dataframe_compat_full_table"
    assert summary["raw_csv_chunk_size"] == 100000
    assert summary["raw_visualization_performed"] is True
    assert summary["raw_visualization_style"] == "legacy_rainbow"
    assert summary["raw_visualization_preview_groups"] == [
        "all_particles",
        "filter_particles_by_sld_gt_1p5",
    ]
    assert summary["raw_visualization_sld_preview_cutoff"] == 1.5
    assert summary["raw_visualization_display_vmax_cap"] == 100.0
    assert summary["density_backend"] == "current_dataframe_compat"
    assert summary["density_query_batch_size"] == 100000
    assert summary["timing_profile_performed"] is False
    assert summary["memory_profile_performed"] is False
    assert summary["euler_convention"] == "extrinsic_zyx"
    assert summary["scipy_euler_sequence"] == "zyx"
    assert summary["euler_angle_columns"] == [
        "raw_ea_zyx_alpha_deg",
        "raw_ea_zyx_beta_deg",
        "raw_ea_zyx_gamma_deg",
    ]
    assert summary["euler_convention_source"] == "cli_default"


def test_run_euler_convention_intrinsic_override_is_recorded(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--euler-convention",
            "intrinsic_zyx",
            "--no-visualize",
        ]
    )
    runner = FakeRunner()

    run_command(args, runner=runner)

    summary = json.loads((output_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["euler_convention"] == "intrinsic_zyx"
    assert summary["scipy_euler_sequence"] == "ZYX"
    assert summary["euler_convention_source"] == "cli_override"
    assert runner.calls[1][0] == "compute_ro_for_phase1_result"
    assert runner.calls[1][2] == {"euler_sequence": "ZYX"}


def test_run_records_explicit_domain_labels(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--ref-domain",
            "rigid_a",
            "--mov-domain",
            "rigid_b",
            "--output-dir",
            str(output_dir),
        ]
    )

    run_command(args, runner=FakeRunner())

    summary = json.loads((output_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["ref_domain"] == "rigid_a"
    assert summary["mov_domain"] == "rigid_b"


def test_run_no_visualize_skips_raw_visualization_artifacts(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--no-visualize",
        ]
    )

    run_command(args, runner=FakeRunner())

    assert not (output_dir / "visualizations" / "raw_default" / "visualization_report.json").exists()
    assert not (output_dir / "visualizations" / "raw_default" / "all_particles").exists()
    assert not (output_dir / "visualizations" / "raw_default" / "filter_particles_by_sld_gt_1p5").exists()
    summary = json.loads((output_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["raw_visualization_performed"] is False
    assert "visualization_report_json" not in summary["output_artifacts"]


def test_run_profile_time_writes_timing_profile(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--profile-time",
            "--no-visualize",
        ]
    )

    run_command(args, runner=FakeRunner())

    profile_path = output_dir / "reports" / "run_timing_profile.json"
    assert profile_path.exists()
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    assert profile["artifact_type"] == "run_timing_profile"
    assert profile["schema_version"] == "1"
    assert profile["total_elapsed_sec"] >= 0.0
    assert {stage["stage"] for stage in profile["stages"]} >= {
        "read_normalize_match",
        "compute_ro",
        "compute_sld",
        "write_npz",
        "write_raw_csv",
        "manifest",
    }
    summary = json.loads((output_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["timing_profile_performed"] is True
    assert summary["timing_profile"] == str(profile_path)


def test_run_profile_memory_writes_memory_profile(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--profile-memory",
            "--no-visualize",
        ]
    )

    run_command(args, runner=FakeRunner())

    profile_path = output_dir / "reports" / "run_memory_profile.json"
    assert profile_path.exists()
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    assert profile["artifact_type"] == "run_memory_profile"
    assert profile["schema_version"] == "1"
    assert profile["metric"] == "rss"
    assert profile["units"] == "bytes"
    assert profile["peak_rss_bytes"] >= 0
    assert {sample["stage"] for sample in profile["samples"]} >= {
        "start",
        "source_loaded",
        "ro_computed",
        "density_computed",
        "raw_npz_written",
        "raw_csv_written",
    }
    summary = json.loads((output_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["memory_profile_performed"] is True
    assert summary["memory_profile"] == str(profile_path)


def test_run_progress_goes_to_stderr_and_stdout_is_final_path(tmp_path, capsys) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir)])

    run_command(args, runner=FakeRunner())

    captured = capsys.readouterr()
    assert captured.out == f"{output_dir}\n"
    assert "[run] stage" in captured.err
    assert str(output_dir) not in captured.err.splitlines()[0]


def test_run_quiet_suppresses_routine_progress(tmp_path, capsys) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--quiet",
        ]
    )

    run_command(args, runner=FakeRunner())

    captured = capsys.readouterr()
    assert captured.out == f"{output_dir}\n"
    assert "[run] stage" not in captured.err


@pytest.mark.parametrize(
    "option",
    ["--raw-csv-chunk-size", "--density-query-batch-size"],
)
def test_run_positive_integer_controls_reject_zero(option) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "run",
                "--ref",
                "ref.cs",
                "--mov",
                "mov.cs",
                option,
                "0",
            ]
        )


def test_run_array_native_backend_fails_clearly(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--run-backend",
            "array_native",
        ]
    )

    with pytest.raises(ValueError, match="array_native is not implemented"):
        run_command(args, runner=FakeRunner())


def test_visualization_report_is_ok_and_records_display_color_field(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir)])

    run_command(args, runner=FakeRunner())

    report = json.loads(
        (
            output_dir
            / "visualizations"
            / "raw_default"
            / "visualization_report.json"
        ).read_text(encoding="utf-8")
    )
    assert report["artifact_type"] == "run_default_visualization_report"
    assert report["preview_only"] is True
    assert report["not_a_selection"] is True
    assert report["visual_style"] == "legacy_rainbow"
    assert report["color_field"] == "sld_raw"
    assert report["sld_preview_cutoff"] == 1.5
    assert report["display_vmax_cap"] == 100.0
    assert report["preview_groups"] == [
        "all_particles",
        "filter_particles_by_sld_gt_1p5",
    ]

    all_report = json.loads(
        (
            output_dir
            / "visualizations"
            / "raw_default"
            / "all_particles"
            / "visualization_report.json"
        ).read_text(encoding="utf-8")
    )
    assert all_report["backend"] == "matplotlib_static"
    assert all_report["color_field"] == "sld_raw"
    assert all_report["display_filter_mode"] == "none"
    assert all_report["display_vmax_cap"] == 100.0
    assert all_report["max_displayed_sld"] == 1.0
    assert all_report["resolved_vmax"] == 1.0
    assert all_report["display_color_vmax"] == 1.0
    assert all_report["full_landscape_table"] == "data/raw_landscape.csv"
    assert all_report["display_table_filename"] == "display_table.csv"

    threshold_report = json.loads(
        (
            output_dir
            / "visualizations"
            / "raw_default"
            / "filter_particles_by_sld_gt_1p5"
            / "visualization_report.json"
        ).read_text(encoding="utf-8")
    )
    assert threshold_report["display_filter_mode"] == "threshold"
    assert threshold_report["display_sld_threshold"] == 1.5
    assert threshold_report["display_vmax_cap"] == 100.0
    assert threshold_report["max_displayed_sld"] is None
    assert threshold_report["resolved_vmax"] is None


def test_run_preview_resolved_vmax_caps_only_above_100(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir)])

    run_command(args, runner=HighSldRunner())

    all_report = json.loads(
        (
            output_dir
            / "visualizations"
            / "raw_default"
            / "all_particles"
            / "visualization_report.json"
        ).read_text(encoding="utf-8")
    )
    threshold_report = json.loads(
        (
            output_dir
            / "visualizations"
            / "raw_default"
            / "filter_particles_by_sld_gt_1p5"
            / "visualization_report.json"
        ).read_text(encoding="utf-8")
    )

    assert all_report["max_displayed_sld"] == 150.0
    assert all_report["resolved_vmax"] == 100.0
    assert all_report["display_color_vmax"] == 100.0
    assert threshold_report["max_displayed_sld"] == 150.0
    assert threshold_report["resolved_vmax"] == 100.0


def test_run_writes_expected_figure_artifacts_and_no_projection_csvs_by_default(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir)])

    run_command(args, runner=FakeRunner())

    expected = [
        "display_table.csv",
        "euler_alpha_beta.png",
        "euler_beta_gamma.png",
        "euler_alpha_gamma.png",
        "euler_3view_projection.png",
        "rotvec_xy.png",
        "rotvec_yz.png",
        "rotvec_xz.png",
        "rotvec_3view_projection.png",
        "landscape_3d_euler.png",
        "landscape_3d_rotvec.png",
    ]
    for group in ("all_particles", "filter_particles_by_sld_gt_1p5"):
        group_dir = output_dir / "visualizations" / "raw_default" / group
        for filename in expected:
            assert (group_dir / filename).exists()
        assert not (group_dir / "analysis_projection_xy.csv").exists()
        assert not (group_dir / "euler_alpha_beta.svg").exists()
        assert not (group_dir / "euler_alpha_beta.pdf").exists()


def test_run_output_formats_alias_writes_requested_formats(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--output-formats",
            "png,pdf,svg",
        ]
    )

    run_command(args, runner=FakeRunner())

    visualization_dir = output_dir / "visualizations" / "raw_default" / "all_particles"
    assert (visualization_dir / "euler_alpha_beta.png").exists()
    assert (visualization_dir / "euler_alpha_beta.pdf").exists()
    assert (visualization_dir / "euler_alpha_beta.svg").exists()


def test_run_write_projection_csvs_writes_debug_projection_csvs(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--write-projection-csvs",
        ]
    )

    run_command(args, runner=FakeRunner())

    assert (
        output_dir
        / "visualizations"
        / "raw_default"
        / "all_particles"
        / "analysis_projection_xy.csv"
    ).exists()
    assert (
        output_dir
        / "visualizations"
        / "raw_default"
        / "all_particles"
        / "analysis_euler_projection_alpha_beta.csv"
    ).exists()


def test_run_display_filter_does_not_change_complete_landscape_outputs(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--display-top-fraction",
            "0.34",
        ]
    )

    run_command(args, runner=FakeRunner())

    raw = pd.read_csv(output_dir / "data" / "raw_landscape.csv")
    display = pd.read_csv(
        output_dir
        / "visualizations"
        / "raw_default"
        / "filter_particles_by_sld_gt_1p5"
        / "display_table.csv"
    )
    assert len(raw) == 3
    assert len(display) == 0


def test_run_command_does_not_select_or_canonicalize_by_default(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir)])
    runner = FakeRunner()

    run_command(args, runner=runner)

    call_names = [call[0] for call in runner.calls]
    assert "canonicalize_density_result" not in call_names
    assert not (output_dir / "canonical" / "default" / "canonical_landscape.npz").exists()
    assert not (output_dir / "canonical" / "default" / "canonicalization_report.json").exists()
    assert not (output_dir / "selection.json").exists()
    summary = json.loads((output_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["selection_performed"] is False


def test_run_command_canonicalizes_only_when_requested(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir), "--canonicalize"]
    )
    runner = FakeRunner()

    run_command(args, runner=runner)

    call_names = [call[0] for call in runner.calls]
    assert "canonicalize_density_result" in call_names
    canonicalize_call = next(call for call in runner.calls if call[0] == "canonicalize_density_result")
    assert canonicalize_call[2].sign_rule == "density_weighted_skewness"
    assert canonicalize_call[2].positive_side == "high_density_skew"
    canonical_dir = output_dir / "canonical" / "default"
    assert (canonical_dir / "canonical_landscape.npz").exists()
    assert (canonical_dir / "canonical_landscape.csv").exists()
    assert (canonical_dir / "canonicalization_report.json").exists()
    assert not (canonical_dir / "visualizations").exists()
    report = json.loads(
        (
            output_dir
            / "visualizations"
            / "raw_default"
            / "all_particles"
            / "visualization_report.json"
        ).read_text(encoding="utf-8")
    )
    assert report["coordinate_source_resolved"] == ["analysis"]
    summary = json.loads((output_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["canonicalization_performed"] is True
    assert "canonical_visualization_report_json" not in summary["output_artifacts"]
    assert "canonical_euler_projection_2d_png" not in summary["output_artifacts"]
    manifest = json.loads((output_dir / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["active_policies"]["canonicalization_policy"]["method"] == "motion_axis_pca"
    assert "canonical_visualization_report_json" not in manifest["output_artifacts"]


def test_visualize_resolves_landscape_from_run_dir(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    visualize_args = parser.parse_args(["visualize", "--run-dir", str(run_dir), "--overwrite"])

    assert visualize_command(visualize_args) == 0

    output_dir = run_dir / "visualizations" / "raw" / "default"
    report = json.loads((output_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert report["status"] == "ok"
    assert report["coordinate_source_resolved"] == ["analysis"]
    assert report["color_field"] == "sld_display"
    assert report["color_field_usage"] == "display_only"
    assert report["representations"] == ["euler", "rotvec"]
    assert report["visual_style"] == "legacy"
    assert report["color_map"] == "rainbow_r"
    assert report["visual_id"] == "default"
    assert report["space"] == "raw"
    assert report["display_density_field"] == "sld_display"
    assert (output_dir / "display_table.csv").exists()
    assert (output_dir / "euler_3view_projection.png").exists()
    assert (output_dir / "landscape_3d_euler.png").exists()
    assert (output_dir / "rotvec_3view_projection.png").exists()
    assert (output_dir / "landscape_3d_rotvec.png").exists()
    assert (output_dir / "distributions_1d" / "distribution_1d_report.json").exists()
    assert (output_dir / "distributions_1d" / "distribution_1d_stats.csv").exists()
    assert (output_dir / "distributions_1d" / "euler_alpha_distribution.png").exists()
    assert (output_dir / "distributions_1d" / "rotvec_x_distribution.png").exists()


def test_run_canonicalize_accepts_low_density_skew_positive_side(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--output-dir",
            str(output_dir),
            "--canonicalize",
            "--skewness-positive-side",
            "low_density_skew",
        ]
    )
    runner = FakeRunner()

    run_command(args, runner=runner)

    canonicalize_call = next(call for call in runner.calls if call[0] == "canonicalize_density_result")
    assert canonicalize_call[2].positive_side == "low_density_skew"


def test_visualize_help_shows_compact_public_surface(capsys) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["visualize", "--help"])

    help_text = capsys.readouterr().out
    for public_option in (
        "--run-dir",
        "--space",
        "--canonical-id",
        "--selection-id",
        "--use-selected-landscape",
        "--visual-id",
        "--representation",
        "--colormap",
        "--range",
        "--top-fraction",
        "--threshold",
        "--formats",
        "--vmin",
        "--vmax",
        "--bins",
        "--hist-mode",
        "--kde-bandwidth",
        "--xlim",
        "--ylim",
        "--overwrite",
    ):
        assert public_option in help_text
    assert "Default: 72" in help_text
    assert "Default: scott" in help_text
    for removed_option in (
        "--output-dir",
        "--landscape",
        "--coordinate-source",
        "--euler-convention",
        "--euler-sequence",
        "--euler-radians",
        "--color-field",
        "--display-density-field",
        "--display-top-fraction",
        "--display-sld-threshold",
        "--display-filter-mode",
        "--display-max-divisor",
        "--visual-style",
        "--color-map",
        "--color-vmin",
        "--color-vmax",
        "--point-size",
        "--point-alpha",
        "--figure-width",
        "--figure-height",
        "--max-points-2d",
        "--max-points-3d",
        "--write-projection-csvs",
        "--generate-histograms",
        "--generate-axis-direction-map",
    ):
        assert removed_option not in help_text


@pytest.mark.parametrize(
    "removed_arg",
    ["--output-dir", "--landscape", "--euler-convention", "--display-top-fraction", "--color-field"],
)
def test_visualize_rejects_removed_public_arguments(removed_arg) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["visualize", "--run-dir", "run", removed_arg, "value"])


def test_visualize_requires_run_dir() -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["visualize"])


def test_visualize_visual_id_controls_raw_output_path(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    args = parser.parse_args(["visualize", "--run-dir", str(run_dir), "--visual-id", "abc"])

    visualize_command(args)

    output_dir = run_dir / "visualizations" / "raw" / "abc"
    assert (output_dir / "display_table.csv").exists()
    report = json.loads((output_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert report["visual_id"] == "abc"


def test_visualize_canonical_writes_compact_path(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    canonicalize_command(
        parser.parse_args(["canonicalize", "--run-dir", str(run_dir), "--canonical-id", "default"])
    )
    args = parser.parse_args(
        ["visualize", "--run-dir", str(run_dir), "--space", "canonical", "--canonical-id", "default"]
    )

    visualize_command(args)

    output_dir = run_dir / "visualizations" / "canonical" / "default" / "default"
    assert (output_dir / "display_table.csv").exists()
    assert (output_dir / "euler_3view_projection.png").exists()
    assert not (run_dir / "visualizations" / "canonical" / "default" / "views").exists()
    assert not (run_dir / "canonical" / "default" / "visualizations" / "default").exists()


def test_visualize_representation_controls_2d_3d_and_1d_outputs(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    args = parser.parse_args(
        ["visualize", "--run-dir", str(run_dir), "--representation", "rotvec", "--visual-id", "rot"]
    )

    visualize_command(args)

    output_dir = run_dir / "visualizations" / "raw" / "rot"
    display = pd.read_csv(output_dir / "display_table.csv")
    assert {"rotvec_x", "rotvec_y", "rotvec_z"}.issubset(display.columns)
    assert "euler_alpha" not in display.columns
    assert (output_dir / "rotvec_3view_projection.png").exists()
    assert not (output_dir / "euler_3view_projection.png").exists()
    assert (output_dir / "distributions_1d" / "rotvec_x_distribution.png").exists()
    assert not (output_dir / "distributions_1d" / "euler_alpha_distribution.png").exists()


def test_visualize_formats_apply_to_2d_3d_and_1d_outputs(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    args = parser.parse_args(
        [
            "visualize",
            "--run-dir",
            str(run_dir),
            "--representation",
            "euler",
            "--visual-id",
            "formats",
            "--formats",
            "png,svg,pdf",
        ]
    )

    visualize_command(args)

    output_dir = run_dir / "visualizations" / "raw" / "formats"
    assert (output_dir / "euler_3view_projection.png").exists()
    assert (output_dir / "euler_3view_projection.svg").exists()
    assert (output_dir / "euler_3view_projection.pdf").exists()
    assert (output_dir / "distributions_1d" / "euler_alpha_distribution.svg").exists()
    assert (output_dir / "distributions_1d" / "euler_alpha_distribution.pdf").exists()


def test_visualize_range_filters_display_rows_without_modifying_parent(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    parent_before = (run_dir / "data" / "raw_landscape.npz").read_bytes()
    args = parser.parse_args(
        [
            "visualize",
            "--run-dir",
            str(run_dir),
            "--representation",
            "rotvec",
            "--range",
            "x:1:",
            "--visual-id",
            "range",
        ]
    )

    visualize_command(args)

    output_dir = run_dir / "visualizations" / "raw" / "range"
    display = pd.read_csv(output_dir / "display_table.csv")
    assert display["particle_key"].tolist() == ["p2"]
    report = json.loads((output_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert report["display_range_filter_applied"] is True
    assert (run_dir / "data" / "raw_landscape.npz").read_bytes() == parent_before


def test_visualize_euler_range_sets_matching_axis_viewport(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    args = parser.parse_args(
        [
            "visualize",
            "--run-dir",
            str(run_dir),
            "--representation",
            "euler",
            "--range",
            "alpha:-20:20",
            "--range",
            "beta:-20:20",
            "--range",
            "gamma:-30:30",
            "--visual-id",
            "euler_range",
        ]
    )

    visualize_command(args)

    output_dir = run_dir / "visualizations" / "raw" / "euler_range"
    report = json.loads((output_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert report["display_range_filter_applied"] is True
    assert report["axis_limits"]["alpha"] == [-20.0, 20.0]
    assert report["axis_limits"]["beta"] == [-20.0, 20.0]
    assert report["axis_limits"]["gamma"] == [-30.0, 30.0]


def test_visualize_density_filters_and_colormap_are_recorded(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=HighSldRunner(),
    )
    args = parser.parse_args(
        [
            "visualize",
            "--run-dir",
            str(run_dir),
            "--visual-id",
            "dense",
            "--top-fraction",
            "0.5",
            "--colormap",
            "viridis",
            "--vmin",
            "0",
            "--vmax",
            "100",
        ]
    )

    visualize_command(args)

    output_dir = run_dir / "visualizations" / "raw" / "dense"
    report = json.loads((output_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert report["display_filter_mode"] == "top_fraction"
    assert report["top_fraction"] == 0.5
    assert report["color_field"] == "sld_display"
    assert report["color_map"] == "viridis"
    assert report["colormap"] == "viridis"
    assert report["vmin"] == 0.0
    assert report["vmax"] == 100.0


def test_visualize_rejects_filter_conflict(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    conflict = parser.parse_args(
        ["visualize", "--run-dir", str(run_dir), "--top-fraction", "0.5", "--threshold", "1.5"]
    )

    with pytest.raises(ValueError, match="Use only one display density filter"):
        visualize_command(conflict)


def test_visualize_1d_distribution_options_are_recorded(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    args = parser.parse_args(
        [
            "visualize",
            "--run-dir",
            str(run_dir),
            "--visual-id",
            "one_d",
            "--bins",
            "12",
            "--hist-mode",
            "percent",
            "--kde-bandwidth",
            "silverman",
            "--xlim=-180:180",
            "--ylim",
            "0:100",
        ]
    )

    visualize_command(args)

    dist_dir = run_dir / "visualizations" / "raw" / "one_d" / "distributions_1d"
    report = json.loads((dist_dir / "distribution_1d_report.json").read_text(encoding="utf-8"))
    stats = pd.read_csv(dist_dir / "distribution_1d_stats.csv")
    assert report["bins"] == 12
    assert report["hist_mode"] == "percent"
    assert report["kde_bandwidth"] == "silverman"
    assert report["xlim"] == [-180.0, 180.0]
    assert report["ylim"] == [0.0, 100.0]
    assert {"axis", "n", "min", "max", "median", "q05", "q25", "q75", "q95", "kde_status"}.issubset(stats.columns)
    assert set(report["generated_axes"]) == {
        "euler_alpha",
        "euler_beta",
        "euler_gamma",
        "rotvec_x",
        "rotvec_y",
        "rotvec_z",
    }


def test_visualize_rejects_invalid_1d_options(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    with pytest.raises(SystemExit):
        parser.parse_args(["visualize", "--run-dir", str(run_dir), "--hist-mode", "density"])
    args = parser.parse_args(["visualize", "--run-dir", str(run_dir), "--kde-bandwidth", "bad"])
    with pytest.raises(ValueError, match="--kde-bandwidth"):
        visualize_command(args)


def test_visualize_selection_id_filters_raw_landscape_from_selected_keys_csv(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    selection_dir = run_dir / "selections" / "subset"
    export_selection(
        _make_subset_selection(selection_id="subset"),
        policy=SelectionExportPolicy(output_dir=selection_dir),
    )
    selection_json = selection_dir / "selection.json"
    selection_payload = json.loads(selection_json.read_text(encoding="utf-8"))
    selection_payload.pop("selected_particle_keys")
    selection_json.write_text(json.dumps(selection_payload), encoding="utf-8")
    args = parser.parse_args(
        [
            "visualize",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "subset",
            "--formats",
            "png",
        ]
    )

    visualize_command(args)

    output_dir = run_dir / "visualizations" / "selections" / "subset" / "parent_raw" / "default"
    display = pd.read_csv(output_dir / "display_table.csv")
    assert display["particle_key"].tolist() == ["p1", "p3"]
    report = json.loads((output_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert report["selection_id"] == "subset"
    assert report["selected_count"] == 2
    assert report["total_count"] == 3
    assert report["selection_filter_applied"] is True
    assert report["landscape_count_before_selection_filter"] == 3
    assert report["landscape_count_after_selection_filter"] == 2
    assert report["missing_selected_key_count"] == 0
    assert report["selection_source_path"].endswith("selected_particle_keys.csv")
    assert report["n_points_input"] == 2


def test_visualize_can_use_selected_derived_landscape(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    select_command(
        parser.parse_args(
            [
                "select",
                "--run-dir",
                str(run_dir),
                "--selection-id",
                "random-subset",
                    "--mode",
                    "random",
                    "--fraction",
                    "0.67",
                    "--seed",
                    "11",
                "--write-selected-landscape",
                "--recompute-sld",
            ]
        )
    )
    args = parser.parse_args(
        [
            "visualize",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "random-subset",
            "--use-selected-landscape",
            "--formats",
            "png",
        ]
    )

    visualize_command(args)

    output_dir = run_dir / "visualizations" / "selections" / "random-subset" / "selected_landscape" / "default"
    report = json.loads((output_dir / "visualization_report.json").read_text(encoding="utf-8"))
    display = pd.read_csv(output_dir / "display_table.csv")
    assert report["selection_id"] == "random-subset"
    assert report["selected_landscape_used"] is True
    assert report["selection_filter_applied"] is False
    assert report["selected_landscape_density_source"] == "recomputed_on_selection"
    assert report["selected_landscape_sld_recomputed"] is True
    selected_landscape_path = Path(report["selected_landscape_path"])
    assert selected_landscape_path.name == "landscape.npz"
    assert selected_landscape_path.parent.name == "selected_landscape"
    assert report["n_points_input"] == len(display)


def test_visualize_selection_id_falls_back_to_selection_json(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    selection_dir = run_dir / "selections" / "json-subset"
    export_selection(
        _make_subset_selection(selection_id="json-subset", selected_particle_keys=("p2",)),
        policy=SelectionExportPolicy(output_dir=selection_dir),
    )
    (selection_dir / "selected_particle_keys.csv").unlink()
    args = parser.parse_args(
        [
            "visualize",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "json-subset",
            "--formats",
            "png",
        ]
    )

    visualize_command(args)

    output_dir = run_dir / "visualizations" / "selections" / "json-subset" / "parent_raw" / "default"
    display = pd.read_csv(output_dir / "display_table.csv")
    assert display["particle_key"].tolist() == ["p2"]
    report = json.loads((output_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert report["selection_source_path"].endswith("selection.json")
    assert report["selected_count"] == 1
    assert report["total_count"] == 3


def test_visualize_selection_id_filters_canonical_landscape(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    canonicalize_command(
        parser.parse_args(["canonicalize", "--run-dir", str(run_dir), "--canonical-id", "default"])
    )
    selection_dir = run_dir / "selections" / "canonical-subset"
    export_selection(
        _make_subset_selection(selection_id="canonical-subset"),
        policy=SelectionExportPolicy(output_dir=selection_dir),
    )
    args = parser.parse_args(
        [
            "visualize",
            "--run-dir",
            str(run_dir),
            "--space",
            "canonical",
            "--canonical-id",
            "default",
            "--selection-id",
            "canonical-subset",
            "--formats",
            "png",
        ]
    )

    visualize_command(args)

    output_dir = (
        run_dir
        / "visualizations"
        / "selections"
        / "canonical-subset"
        / "parent_canonical"
        / "default"
        / "default"
    )
    display = pd.read_csv(output_dir / "display_table.csv")
    assert display["particle_key"].tolist() == ["p1", "p3"]
    report = json.loads((output_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert report["selection_id"] == "canonical-subset"
    assert report["coordinate_source_resolved"] == ["canonical"]


def test_visualize_selection_id_requires_run_dir() -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["visualize", "--selection-id", "subset"])


def test_visualize_selection_id_missing_artifact_fails_clearly(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    args = parser.parse_args(
        ["visualize", "--run-dir", str(run_dir), "--selection-id", "missing"]
    )

    with pytest.raises(ValueError, match="Selection directory does not exist"):
        visualize_command(args)


def test_visualize_selection_id_missing_keys_fail_fast(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    selection_dir = run_dir / "selections" / "bad-subset"
    export_selection(
        _make_subset_selection(selection_id="bad-subset", selected_particle_keys=("p9",)),
        policy=SelectionExportPolicy(output_dir=selection_dir),
    )
    args = parser.parse_args(
        ["visualize", "--run-dir", str(run_dir), "--selection-id", "bad-subset"]
    )

    with pytest.raises(ValueError, match="absent from the requested landscape"):
        visualize_command(args)


def test_visualize_overwrite_behavior(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_command(
        parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]),
        runner=FakeRunner(),
    )
    args = parser.parse_args(
        ["visualize", "--run-dir", str(run_dir), "--visual-id", "same"]
    )
    overwrite_args = parser.parse_args(
        ["visualize", "--run-dir", str(run_dir), "--visual-id", "same", "--overwrite"]
    )

    assert visualize_command(args) == 0
    with pytest.raises(FileExistsError, match="Output path already exists"):
        visualize_command(args)
    assert visualize_command(overwrite_args) == 0
    output_dir = run_dir / "visualizations" / "raw" / "same"
    payload = json.loads((output_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert payload["status"] == "ok"


def test_canonicalize_loads_landscape_and_writes_outputs(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
        ]
    )

    assert canonicalize_command(args) == 0

    assert (output_dir / "canonical_landscape.npz").exists()
    assert (output_dir / "canonical_landscape.csv").exists()
    assert (output_dir / "canonicalization_report.json").exists()
    assert (output_dir / "canonicalize_summary.json").exists()
    assert (output_dir / "canonical_frame.json").exists()
    assert (output_dir / "canonical_frame.npz").exists()
    assert (output_dir / "visualizations" / "visualization_report.json").exists()
    assert (output_dir / "visualizations" / "all_particles" / "visualization_report.json").exists()
    assert not (output_dir / "selection.json").exists()


def test_canonicalize_summary_records_source_and_artifacts(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
        ]
    )

    canonicalize_command(args)

    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    assert summary["artifact_type"] == "canonicalize_summary"
    assert summary["schema_version"] == "1"
    assert "timestamp" in summary
    assert summary["source_landscape_path"] == str(landscape_path)
    assert summary["source_landscape_artifact_type"] == "landscape"
    assert summary["source_landscape_schema_version"] == "1"
    assert summary["source_landscape_row_count"] == 3
    assert summary["canonicalization_policy"]["method"] == "motion_axis_pca"
    assert summary["canonicalization_policy"]["axis_assignment"] == "pc123_to_alpha_beta_gamma"
    assert summary["canonicalization_policy"]["sign_rule"] == "density_weighted_skewness"
    assert summary["canonicalization_policy"]["positive_side"] == "low_density_skew"
    assert summary["canonicalization_policy"]["sign_weight_field"] == "sld_raw"
    assert summary["canonicalization_policy"]["fit_top_fraction"] == 0.40
    assert summary["fit_top_fraction"] == 0.40
    assert summary["canonicalization_performed"] is True
    assert summary["frame_applied"] is False
    assert summary["frame_source_path"] is None
    assert summary["canonicalization_backend"] == "dataframe_compat"
    assert summary["csv_performed"] is True
    assert summary["csv_backend"] == "dataframe_compat"
    assert summary["csv_chunk_size"] == 100000
    assert summary["visualization_performed"] is True
    assert summary["visualization_report"] == str(
        output_dir / "visualizations" / "visualization_report.json"
    )
    assert "canonical_visualization_hint" not in summary
    assert summary["canonical_frame_json"] == str(output_dir / "canonical_frame.json")
    assert summary["canonical_frame_npz"] == str(output_dir / "canonical_frame.npz")
    assert summary["memory_profile_performed"] is False
    assert summary["memory_profile"] is None
    assert summary["selection_performed"] is False
    assert summary["output_artifacts"]["canonical_landscape_npz"] == str(
        output_dir / "canonical_landscape.npz"
    )
    assert summary["output_artifacts"]["canonical_landscape_csv"] == str(
        output_dir / "canonical_landscape.csv"
    )
    assert summary["output_artifacts"]["canonical_frame_json"] == str(
        output_dir / "canonical_frame.json"
    )
    assert summary["output_artifacts"]["canonical_visualization_report_json"] == str(
        output_dir / "visualizations" / "visualization_report.json"
    )
    assert summary["updated_star_cs_export_performed"] is False
    assert summary["composite_map_export_performed"] is False
    assert "memory_profile_json" not in summary["output_artifacts"]


def test_canonicalize_accepts_low_density_skew_positive_side(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
            "--skewness-positive-side",
            "low_density_skew",
        ]
    )

    canonicalize_command(args)

    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    report = json.loads((output_dir / "canonicalization_report.json").read_text(encoding="utf-8"))
    assert summary["canonicalization_policy"]["positive_side"] == "low_density_skew"
    assert report["report"]["positive_side"] == "low_density_skew"
    assert report["report"]["sign_rule"] == "density_weighted_skewness"


def test_canonicalize_positive_side_high_maps_to_high_density_skew(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
            "--positive-side",
            "high",
        ]
    )

    canonicalize_command(args)

    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    assert summary["canonicalization_policy"]["positive_side"] == "high_density_skew"


def test_canonicalize_fit_top_fraction_is_recorded(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
            "--fit-top-fraction",
            "0.50",
        ]
    )

    canonicalize_command(args)

    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    report = json.loads((output_dir / "canonicalization_report.json").read_text(encoding="utf-8"))
    assert summary["fit_top_fraction"] == 0.50
    assert summary["canonicalization_policy"]["fit_top_fraction"] == 0.50
    assert report["report"]["fit_top_fraction"] == 0.50


def test_canonicalize_fit_top_alias_is_recorded(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
            "--fit-top",
            "0.50",
        ]
    )

    canonicalize_command(args)

    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    assert summary["fit_top_fraction"] == 0.50


@pytest.mark.parametrize("value", ["0", "-0.1", "1.1"])
def test_canonicalize_invalid_fit_top_fraction_fails_clearly(tmp_path, value) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
            "--fit-top-fraction",
            value,
        ]
    )

    with pytest.raises(ValueError, match="--fit-top-fraction must be > 0 and <= 1"):
        canonicalize_command(args)


def test_canonicalize_help_omits_visualization_options(capsys) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["canonicalize", "--help"])

    help_text = capsys.readouterr().out
    assert "--run-dir" in help_text
    assert "--canonical-id" in help_text
    assert "--fit-top" in help_text
    assert "--fit-top-fraction" in help_text
    assert "--positive-side" in help_text
    assert "--use-frame" in help_text
    assert "--no-visualize" in help_text
    assert "--landscape" not in help_text
    assert "--output-dir" not in help_text
    assert "--output " not in help_text
    assert "--sign-rule" not in help_text
    assert "--euler-convention" not in help_text
    assert "--visualize" not in help_text
    assert "--display-top-fraction" not in help_text
    assert "--display-sld-threshold" not in help_text
    assert "--display-density-field" not in help_text
    assert "--max-points-2d" not in help_text
    assert "--max-points-3d" not in help_text
    assert "--write-projection-csvs" not in help_text
    assert "--visual-style" not in help_text
    assert "--formats" not in help_text
    assert "--output-formats" not in help_text


def test_canonicalize_requires_run_dir_for_public_entrypoint() -> None:
    parser = build_parser()
    args = parser.parse_args(["canonicalize"])

    with pytest.raises(ValueError, match="canonicalize requires --run-dir"):
        canonicalize_command(args)


def test_canonicalize_profile_memory_writes_rss_profile(tmp_path, monkeypatch) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    rss_state = {"value": 90_000_000}

    def next_rss():
        rss_state["value"] += 10_000_000
        return rss_state["value"], "test_rss"

    monkeypatch.setattr(
        "cryorole.cli.main._current_rss_bytes",
        next_rss,
    )
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
            "--profile-memory",
        ]
    )

    canonicalize_command(args)

    profile_path = output_dir / "canonical_memory_profile.json"
    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    profile = json.loads(profile_path.read_text(encoding="utf-8"))
    assert profile_path.exists()
    assert summary["memory_profile_performed"] is True
    assert summary["memory_profile"] == str(profile_path)
    assert summary["output_artifacts"]["memory_profile_json"] == str(profile_path)
    assert profile["artifact_type"] == "canonical_memory_profile"
    assert profile["metric"] == "rss"
    assert profile["backend"] == "test_rss"
    assert profile["peak_rss_bytes"] == 180_000_000
    assert profile["samples"][0]["stage"] == "start"
    assert profile["samples"][-1]["stage"] == "canonicalize_summary_written"
    assert profile["samples"][-1]["delta_from_start_bytes"] == 80_000_000


def test_canonicalize_does_not_write_memory_profile_by_default(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
        ]
    )

    canonicalize_command(args)

    assert not (output_dir / "canonical_memory_profile.json").exists()


def test_canonicalize_output_landscape_contains_canonical_coordinates(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
        ]
    )

    canonicalize_command(args)

    with np.load(output_dir / "canonical_landscape.npz", allow_pickle=False) as payload:
        assert str(payload["artifact_type"].item()) == "canonical_landscape"
        assert "coordinates_canonical" in payload
        np.testing.assert_allclose(payload["coordinates_analysis"][0], [0.0, 0.0, 0.0])
    canonical_csv = pd.read_csv(output_dir / "canonical_landscape.csv")
    assert {"raw_rotvec_x", "canonical_rotvec_x", "canonical_euler_alpha"}.issubset(
        canonical_csv.columns
    )


def test_canonicalize_writes_default_visualization_artifacts(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
        ]
    )

    canonicalize_command(args)

    canonical_csv = pd.read_csv(output_dir / "canonical_landscape.csv")
    assert len(canonical_csv) == 3
    assert not (output_dir / "visualization_report.json").exists()
    assert not (output_dir / "display_coordinates.csv").exists()
    visualization_dir = output_dir / "visualizations"
    assert not (visualization_dir / "default").exists()
    aggregate = json.loads((visualization_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert aggregate["artifact_type"] == "canonical_default_visualization_report"
    assert aggregate["preview_groups"] == [
        "all_particles",
        "filter_particles_by_sld_gt_1p5",
        "filter_particles_by_top_sld_40pct",
    ]
    assert aggregate["fit_top_fraction"] == 0.40
    assert aggregate["top_sld_preview_group"] == "filter_particles_by_top_sld_40pct"
    for group in aggregate["preview_groups"]:
        group_report = _assert_default_visualization_group(visualization_dir / group)
        assert group_report["preview_only"] is True
        assert group_report["not_a_selection"] is True
        assert group_report["display_vmax_cap"] == 100.0
        assert "max_displayed_sld" in group_report
        assert "resolved_vmax" in group_report
    top_report = json.loads(
        (visualization_dir / "filter_particles_by_top_sld_40pct" / "visualization_report.json").read_text(
            encoding="utf-8"
        )
    )
    assert top_report["top_sld_fraction"] == 0.40
    assert not (output_dir / "selections").exists()


def test_canonicalize_fit_top_controls_top_sld_preview_name(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(sld_values=[1.0, 2.0, 3.0]), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
            "--fit-top",
            "0.50",
        ]
    )

    canonicalize_command(args)

    visualization_dir = output_dir / "visualizations"
    aggregate = json.loads((visualization_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert aggregate["fit_top_fraction"] == 0.50
    assert aggregate["top_sld_preview_group"] == "filter_particles_by_top_sld_50pct"
    top_report = _assert_default_visualization_group(
        visualization_dir / "filter_particles_by_top_sld_50pct"
    )
    assert top_report["display_filter_mode"] == "top_fraction"
    assert top_report["top_sld_fraction"] == 0.50


def test_canonicalize_no_visualize_skips_default_visualization(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
            "--no-visualize",
        ]
    )

    canonicalize_command(args)

    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    assert summary["visualization_performed"] is False
    assert summary["visualization_report"] is None
    assert not (output_dir / "visualizations").exists()


def test_canonicalize_rejects_removed_visualize_option() -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["canonicalize", "--landscape", "l.json", "--visualize"])


def test_canonicalize_does_not_modify_input_landscape_artifact(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    before = landscape_path.read_text(encoding="utf-8")
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
        ]
    )

    canonicalize_command(args)

    assert landscape_path.read_text(encoding="utf-8") == before


def test_canonicalize_existing_output_directory_fails_by_default(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    output_dir.mkdir()
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
        ]
    )

    with pytest.raises(FileExistsError, match="Canonicalize output directory already exists"):
        canonicalize_command(args)


def test_canonicalize_overwrite_permits_existing_output_directory(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    output_dir.mkdir()
    (output_dir / "canonical_landscape.npz").write_text("old\n", encoding="utf-8")
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(landscape_path),
            "--output-dir",
            str(output_dir),
            "--overwrite",
        ]
    )

    assert canonicalize_command(args) == 0
    with np.load(output_dir / "canonical_landscape.npz", allow_pickle=False) as payload:
        assert str(payload["artifact_type"].item()) == "canonical_landscape"


def test_canonicalize_legacy_output_is_directory_alias(tmp_path) -> None:
    landscape_path = tmp_path / "landscape.json"
    output_dir = tmp_path / "canonicalize"
    write_landscape_json(_make_landscape(), landscape_path)
    parser = build_parser()
    args = parser.parse_args(
        ["canonicalize", "--landscape", str(landscape_path), "--output", str(output_dir)]
    )

    canonicalize_command(args)

    assert (output_dir / "canonical_landscape.npz").exists()


def test_canonicalize_resolves_raw_landscape_from_run_dir(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    canonicalize_args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--canonical-id",
            "default",
        ]
    )

    canonicalize_command(canonicalize_args)

    output_dir = run_dir / "canonical" / "default"
    assert (output_dir / "canonical_landscape.npz").exists()
    assert (output_dir / "canonical_landscape.csv").exists()
    assert not (output_dir / "visualizations" / "display_table.csv").exists()
    assert (
        run_dir
        / "visualizations"
        / "canonical"
        / "default"
        / "all_particles"
        / "display_table.csv"
    ).exists()
    assert not (run_dir / "visualizations" / "canonical" / "default" / "default").exists()
    assert (
        run_dir / "visualizations" / "canonical" / "default" / "visualization_report.json"
    ).exists()
    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    assert summary["canonicalization_backend"] == "array_native"
    assert summary["csv_backend"] == "array_native_chunked"
    assert summary["visualization_performed"] is True


def test_canonicalize_canonical_id_controls_visualization_path(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    canonicalize_args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--canonical-id",
            "abc",
        ]
    )

    canonicalize_command(canonicalize_args)

    visualization_dir = run_dir / "visualizations" / "canonical" / "abc"
    assert (visualization_dir / "visualization_report.json").exists()
    assert (visualization_dir / "all_particles" / "display_table.csv").exists()
    assert not (run_dir / "visualizations" / "canonical" / "abc" / "default").exists()


def test_canonicalize_run_dir_array_backend_does_not_call_dataframe_reader(
    tmp_path,
    monkeypatch,
) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    monkeypatch.setattr(
        "cryorole.cli.main.read_landscape",
        lambda *_args, **_kwargs: pytest.fail("NPZ canonicalize must not read DataFrame Landscape"),
    )
    canonicalize_args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--canonical-id",
            "default",
        ]
    )

    canonicalize_command(canonicalize_args)

    summary = json.loads(
        (run_dir / "canonical" / "default" / "canonicalize_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert summary["canonicalization_backend"] == "array_native"


def test_canonicalize_explicit_npz_uses_array_native_backend(tmp_path, monkeypatch) -> None:
    run_dir = tmp_path / "run"
    output_dir = tmp_path / "canonicalize"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    monkeypatch.setattr(
        "cryorole.cli.main.read_landscape",
        lambda *_args, **_kwargs: pytest.fail("Explicit NPZ canonicalize must not read DataFrame Landscape"),
    )
    canonicalize_args = parser.parse_args(
        [
            "canonicalize",
            "--landscape",
            str(run_dir / "data" / "raw_landscape.npz"),
            "--output-dir",
            str(output_dir),
        ]
    )

    canonicalize_command(canonicalize_args)

    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    assert summary["canonicalization_backend"] == "array_native"


def test_canonicalize_writes_npz_before_chunked_csv_export(tmp_path, monkeypatch) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())

    def checked_writer(npz_path, csv_path, **kwargs):
        assert npz_path.exists()
        return real_write_canonical_landscape_csv_from_npz(npz_path, csv_path, **kwargs)

    monkeypatch.setattr(
        "cryorole.cli.main.write_canonical_landscape_csv_from_npz",
        checked_writer,
    )

    canonicalize_args = parser.parse_args(
        ["canonicalize", "--run-dir", str(run_dir), "--canonical-id", "default"]
    )

    canonicalize_command(canonicalize_args)


def test_canonicalize_no_csv_skips_csv_and_records_status(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    canonicalize_args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--canonical-id",
            "default",
            "--no-csv",
        ]
    )

    canonicalize_command(canonicalize_args)

    output_dir = run_dir / "canonical" / "default"
    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    assert (output_dir / "canonical_landscape.npz").exists()
    assert not (output_dir / "canonical_landscape.csv").exists()
    assert summary["csv_performed"] is False
    assert summary["csv_backend"] == "skipped"
    assert summary["output_artifacts"]["canonical_landscape_csv"] is None


def test_canonicalize_use_frame_applies_existing_frame_and_skips_fitting(tmp_path, monkeypatch) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    base_args = parser.parse_args(
        ["canonicalize", "--run-dir", str(run_dir), "--canonical-id", "base", "--no-visualize"]
    )
    canonicalize_command(base_args)
    frame_path = run_dir / "canonical" / "base" / "canonical_frame.json"

    monkeypatch.setattr(
        "cryorole.cli.main.canonicalize_landscape_arrays",
        lambda *_args, **_kwargs: pytest.fail("--use-frame must skip PCA fitting"),
    )
    aligned_args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--canonical-id",
            "aligned",
            "--use-frame",
            str(frame_path),
            "--no-visualize",
        ]
    )

    canonicalize_command(aligned_args)

    output_dir = run_dir / "canonical" / "aligned"
    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    report = json.loads((output_dir / "canonicalization_report.json").read_text(encoding="utf-8"))
    assert (output_dir / "canonical_landscape.npz").exists()
    assert (output_dir / "canonical_frame.json").exists()
    assert summary["frame_applied"] is True
    assert summary["frame_source_path"] == str(frame_path)
    assert summary["canonicalization_performed"] is False
    assert report["report"]["frame_applied"] is True
    assert report["report"]["frame_source_path"] == str(frame_path)


def test_canonicalize_use_frame_uses_frame_fit_fraction_for_preview(tmp_path, monkeypatch) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    base_args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--canonical-id",
            "base",
            "--fit-top",
            "0.50",
            "--no-visualize",
        ]
    )
    canonicalize_command(base_args)
    frame_path = run_dir / "canonical" / "base" / "canonical_frame.json"

    monkeypatch.setattr(
        "cryorole.cli.main.canonicalize_landscape_arrays",
        lambda *_args, **_kwargs: pytest.fail("--use-frame must skip PCA fitting"),
    )
    aligned_args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--canonical-id",
            "aligned",
            "--use-frame",
            str(frame_path),
        ]
    )

    canonicalize_command(aligned_args)

    visualization_dir = run_dir / "visualizations" / "canonical" / "aligned"
    aggregate = json.loads((visualization_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert aggregate["fit_top_fraction"] == 0.50
    assert aggregate["top_sld_preview_group"] == "filter_particles_by_top_sld_50pct"
    assert (visualization_dir / "filter_particles_by_top_sld_50pct" / "display_table.csv").exists()


def test_canonicalize_use_frame_missing_fit_fraction_warns_and_uses_default_preview(
    tmp_path,
    monkeypatch,
) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    frame_path = tmp_path / "frame_without_fraction.json"
    frame_path.write_text(
        json.dumps(
            {
                "canonical_transform": np.eye(3).tolist(),
                "transform_direction": "canonical_rv = raw_rv @ canonical_transform",
                "coordinate_space": "rotvec_ro_radians",
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "cryorole.cli.main.canonicalize_landscape_arrays",
        lambda *_args, **_kwargs: pytest.fail("--use-frame must skip PCA fitting"),
    )
    args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--canonical-id",
            "fallback",
            "--use-frame",
            str(frame_path),
        ]
    )

    canonicalize_command(args)

    visualization_dir = run_dir / "visualizations" / "canonical" / "fallback"
    aggregate = json.loads((visualization_dir / "visualization_report.json").read_text(encoding="utf-8"))
    assert aggregate["fit_top_fraction"] == 0.40
    assert aggregate["top_sld_preview_group"] == "filter_particles_by_top_sld_40pct"
    assert aggregate["warnings"]
    assert (visualization_dir / "filter_particles_by_top_sld_40pct" / "display_table.csv").exists()


def test_canonicalize_use_frame_rejects_invalid_frame_shape(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    frame_path = tmp_path / "bad_frame.json"
    frame_path.write_text(
        json.dumps(
            {
                "canonical_transform": [[1, 0], [0, 1]],
                "transform_direction": "canonical_rv = raw_rv @ canonical_transform",
                "coordinate_space": "rotvec_ro_radians",
            }
        ),
        encoding="utf-8",
    )
    args = parser.parse_args(
        ["canonicalize", "--run-dir", str(run_dir), "--use-frame", str(frame_path)]
    )

    with pytest.raises(ValueError, match="canonical_transform must be shape"):
        canonicalize_command(args)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("transform_direction", "wrong", "transform_direction"),
        ("coordinate_space", "wrong", "coordinate_space"),
    ],
)
def test_canonicalize_use_frame_rejects_invalid_frame_metadata(
    tmp_path,
    field,
    value,
    message,
) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    payload = {
        "canonical_transform": np.eye(3).tolist(),
        "transform_direction": "canonical_rv = raw_rv @ canonical_transform",
        "coordinate_space": "rotvec_ro_radians",
    }
    payload[field] = value
    frame_path = tmp_path / "bad_frame.json"
    frame_path.write_text(json.dumps(payload), encoding="utf-8")
    args = parser.parse_args(
        ["canonicalize", "--run-dir", str(run_dir), "--use-frame", str(frame_path)]
    )

    with pytest.raises(ValueError, match=message):
        canonicalize_command(args)


@pytest.mark.parametrize(
    "extra_args",
    [
        ["--fit-top", "0.5"],
        ["--fit-top-fraction", "0.5"],
        ["--positive-side", "high"],
    ],
)
def test_canonicalize_use_frame_conflicts_with_explicit_fitting_controls(
    tmp_path,
    extra_args,
) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    frame_path = tmp_path / "frame.json"
    parser = build_parser()
    args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--use-frame",
            str(frame_path),
            *extra_args,
        ]
    )

    with pytest.raises(ValueError, match="--use-frame conflicts"):
        canonicalize_command(args)


def test_canonicalize_csv_chunk_size_is_recorded(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    canonicalize_args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--canonical-id",
            "default",
            "--csv-chunk-size",
            "2",
        ]
    )

    canonicalize_command(canonicalize_args)

    output_dir = run_dir / "canonical" / "default"
    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    assert summary["csv_chunk_size"] == 2
    assert len(pd.read_csv(output_dir / "canonical_landscape.csv")) == 3


def test_canonicalize_no_visualize_does_not_write_run_dir_visualization(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    canonicalize_args = parser.parse_args(
        [
            "canonicalize",
            "--run-dir",
            str(run_dir),
            "--canonical-id",
            "default",
            "--no-visualize",
        ]
    )

    canonicalize_command(canonicalize_args)

    output_dir = run_dir / "canonical" / "default"
    assert not (output_dir / "visualizations").exists()
    assert not (run_dir / "visualizations" / "canonical" / "default").exists()
    summary = json.loads((output_dir / "canonicalize_summary.json").read_text(encoding="utf-8"))
    assert summary["visualization_performed"] is False


def test_visualize_can_read_canonical_npz_from_run_dir_after_canonicalize(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    canonicalize_command(
        parser.parse_args(["canonicalize", "--run-dir", str(run_dir), "--canonical-id", "default"])
    )
    visualize_args = parser.parse_args(
        [
            "visualize",
            "--run-dir",
            str(run_dir),
            "--space",
            "canonical",
            "--canonical-id",
            "default",
            "--overwrite",
        ]
    )

    visualize_command(visualize_args)

    assert (
        run_dir
        / "visualizations"
        / "canonical"
        / "default"
        / "default"
        / "euler_3view_projection.png"
    ).exists()
    assert not (
        run_dir
        / "canonical"
        / "default"
        / "visualizations"
        / "default"
    ).exists()


def test_select_can_read_canonical_npz_from_run_dir_after_canonicalize(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    canonicalize_command(
        parser.parse_args(["canonicalize", "--run-dir", str(run_dir), "--canonical-id", "default"])
    )
    select_args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--space",
            "canonical",
            "--canonical-id",
            "default",
            "--selection-id",
            "canonical-core",
            "--mode",
            "threshold",
            "--sld-min",
            "0",
        ]
    )

    select_command(select_args)

    assert (
        run_dir / "selections" / "canonical-core" / "selected_particle_keys.csv"
    ).exists()


def test_run_command_builds_explicit_row_aligned_identity_policy(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.star",
            "--mov",
            "mov.star",
            "--row-aligned",
            "--output-dir",
            str(output_dir),
        ]
    )
    runner = FakeRunner()

    run_command(args, runner=runner)

    kwargs = runner.calls[0][2]
    assert kwargs["identity_policy_a"].identity_mode == "row_aligned"
    assert kwargs["identity_policy_b"].identity_mode == "row_aligned"


def test_run_command_defaults_star_identity_to_image_name(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.star",
            "--mov",
            "mov.star",
            "--output-dir",
            str(output_dir),
        ]
    )
    runner = FakeRunner()

    run_command(args, runner=runner)

    kwargs = runner.calls[0][2]
    assert kwargs["identity_policy_a"].identity_mode == "relion_image_name"
    assert kwargs["identity_policy_b"].identity_mode == "relion_image_name"


def test_run_command_row_aligned_applies_to_cryosparc_inputs(tmp_path) -> None:
    output_dir = tmp_path / "run"
    parser = build_parser()
    args = parser.parse_args(
        [
            "run",
            "--ref",
            "ref.cs",
            "--mov",
            "mov.cs",
            "--row-aligned",
            "--output-dir",
            str(output_dir),
        ]
    )
    runner = FakeRunner()

    run_command(args, runner=runner)

    kwargs = runner.calls[0][2]
    assert kwargs["identity_policy_a"].identity_mode == "row_aligned"
    assert kwargs["identity_policy_b"].identity_mode == "row_aligned"


def test_run_existing_output_directory_fails_by_default(tmp_path) -> None:
    output_dir = tmp_path / "run"
    output_dir.mkdir()
    parser = build_parser()
    args = parser.parse_args(["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir)])

    with pytest.raises(FileExistsError, match="Run output directory already exists"):
        run_command(args, runner=FakeRunner())


def test_run_overwrite_permits_existing_output_directory(tmp_path) -> None:
    output_dir = tmp_path / "run"
    output_dir.mkdir()
    (output_dir / "data").mkdir()
    (output_dir / "data" / "raw_landscape.npz").write_text("old\n", encoding="utf-8")
    parser = build_parser()
    args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(output_dir), "--overwrite"]
    )

    assert run_command(args, runner=FakeRunner()) == 0
    with np.load(output_dir / "data" / "raw_landscape.npz", allow_pickle=False) as payload:
        assert str(payload["artifact_type"].item()) == "raw_landscape"


def test_main_prints_clean_value_error(capsys) -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(["run", "--ref", "ref.star", "--mov", "mov.star", "--identity-mode", "relion_user_columns"])

    captured = capsys.readouterr()
    assert excinfo.value.code == 2
    assert "cryorole: error: relion_user_columns identity requires --identity-column" in captured.err
    assert "Traceback" not in captured.err


def test_select_command_rejects_sld_display_support_field() -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "select",
                "--run-dir",
                "run",
                "--selection-id",
                "sel",
                "--density-support-field",
                "sld_display",
            ]
        )


def test_select_loads_landscape_and_writes_selection_outputs(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--center",
            "0",
            "0",
            "0",
            "--radius",
            "30",
        ]
    )

    assert select_command(args) == 0

    assert (output_dir / "selection.json").exists()
    assert (output_dir / "selected_particle_keys.csv").exists()
    assert (output_dir / "export_report.json").exists()
    assert (output_dir / "selection_summary.json").exists()


def test_select_resolves_raw_landscape_from_run_dir(tmp_path) -> None:
    run_dir = tmp_path / "run"
    parser = build_parser()
    run_args = parser.parse_args(
        ["run", "--ref", "ref.cs", "--mov", "mov.cs", "--output-dir", str(run_dir)]
    )
    run_command(run_args, runner=FakeRunner())
    select_args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "core",
            "--center",
            "0",
            "0",
            "0",
            "--radius",
            "30",
        ]
    )

    select_command(select_args)

    output_dir = run_dir / "selections" / "core"
    assert (output_dir / "selection.json").exists()
    summary = json.loads((output_dir / "selection_summary.json").read_text(encoding="utf-8"))
    assert summary["landscape_path"] == str(run_dir / "data" / "raw_landscape.npz")


def test_select_summary_includes_landscape_metadata_policy_and_export_report(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--mode",
            "threshold",
            "--sld-min",
            "1.0",
        ]
    )

    select_command(args)

    summary = json.loads((output_dir / "selection_summary.json").read_text(encoding="utf-8"))
    assert summary["artifact_type"] == "selection_summary"
    assert summary["schema_version"] == "1"
    assert "timestamp" in summary
    assert summary["overwrite"] is False
    assert summary["landscape_path"] == str(run_dir / "data" / "raw_landscape.npz")
    assert summary["parent_space"] == "raw"
    assert summary["candidate_count"] == 3
    assert summary["landscape_artifact_type"] == "raw_landscape"
    assert summary["landscape_schema_version"] == "1"
    assert summary["landscape_row_count"] == 3
    assert summary["selection_policy"]["selection_mode"] == "threshold_by_density"
    assert summary["selection_policy"]["sld_min"] == 1.0
    assert summary["selection_policy"]["density_support_field"] == "sld_raw"
    assert summary["export_report_output_paths"]["selection_json"] == str(
        output_dir / "selection.json"
    )
    assert summary["export_report_row_counts"]["selected_particle_keys_csv"] == 3


def test_select_default_radius_preserves_selected_order(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--center",
            "0",
            "0",
            "0",
            "--radius",
            "30",
        ]
    )

    select_command(args)

    rows = (output_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    assert rows == ["selection_rank,particle_key", "1,p1"]


def test_select_threshold_includes_display_outlier_rows_by_default(tmp_path) -> None:
    landscape = _make_landscape()
    landscape.data["sld_raw"] = [100.0, 9.0, 8.0]
    landscape.data["sld_display_is_outlier"] = [True, False, False]
    run_dir = _make_select_run_bundle(tmp_path, landscape)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--mode",
            "threshold",
            "--sld-min",
            "50",
        ]
    )

    select_command(args)

    rows = (output_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    summary = json.loads((output_dir / "selection_summary.json").read_text(encoding="utf-8"))
    selection_payload = json.loads((output_dir / "selection.json").read_text(encoding="utf-8"))
    assert rows == ["selection_rank,particle_key", "1,p1"]
    assert summary["density_artifact_policy"] == "include_all"
    assert summary["density_artifact_excluded_count"] == 0
    assert selection_payload["density_artifact_candidate_count_after"] == 3


def test_select_threshold_by_density_works(tmp_path) -> None:
    landscape = _make_landscape()
    landscape.data["sld_raw"] = [0.5, 2.0, 3.0]
    landscape.data["sld_display"] = [100.0, 100.0, 0.0]
    run_dir = _make_select_run_bundle(tmp_path, landscape)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--mode",
            "threshold",
            "--sld-min",
            "2.0",
            "--sld-max",
            "2.5",
        ]
    )

    select_command(args)

    rows = (output_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    payload = json.loads((output_dir / "selection.json").read_text(encoding="utf-8"))
    assert rows == ["selection_rank,particle_key", "1,p2"]
    assert payload["density_support_field"] == "sld_raw"
    assert payload["sld_min"] == 2.0
    assert payload["sld_max"] == 2.5


def test_select_threshold_accepts_sld_max_only_and_rejects_bad_window(tmp_path) -> None:
    landscape = _make_landscape()
    landscape.data["sld_raw"] = [0.5, 2.0, 3.0]
    run_dir = _make_select_run_bundle(tmp_path, landscape)
    parser = build_parser()

    select_command(
        parser.parse_args(
            [
                "select",
                "--run-dir",
                str(run_dir),
                "--selection-id",
                "max_only",
                "--mode",
                "threshold",
                "--sld-max",
                "2.0",
            ]
        )
    )
    rows = (run_dir / "selections" / "max_only" / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    assert rows == ["selection_rank,particle_key", "1,p1", "2,p2"]

    bad_args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "bad_window",
            "--mode",
            "threshold",
            "--sld-min",
            "3.0",
            "--sld-max",
            "1.0",
        ]
    )
    with pytest.raises(ValueError, match="sld_min"):
        select_command(bad_args)


def test_select_random_mode_records_policy_and_is_reproducible(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    first_dir = run_dir / "selections" / "selection_first"
    second_dir = run_dir / "selections" / "selection_second"
    parser = build_parser()
    base_args = [
        "select",
        "--run-dir",
        str(run_dir),
        "--mode",
        "random",
        "--fraction",
        "0.67",
        "--seed",
        "7",
    ]

    select_command(parser.parse_args([*base_args, "--selection-id", "selection_first"]))
    select_command(parser.parse_args([*base_args, "--selection-id", "selection_second"]))

    first_rows = (first_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    second_rows = (second_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    summary = json.loads((first_dir / "selection_summary.json").read_text(encoding="utf-8"))
    payload = json.loads((first_dir / "selection.json").read_text(encoding="utf-8"))
    assert first_rows == second_rows
    assert summary["selection_mode"] == "random"
    assert summary["random_fraction"] == 0.67
    assert summary["random_seed"] == 7
    assert summary["random_candidate_count"] == 3
    assert payload["random_fraction"] == 0.67
    assert payload["random_seed"] == 7
    assert payload["random_candidate_count"] == 3


def test_select_random_mode_rejects_invalid_fraction(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--mode",
            "random",
            "--fraction",
            "0",
        ]
    )

    with pytest.raises(ValueError, match="random_fraction"):
        select_command(args)


def test_select_metadata_value_from_run_time_ref_star(tmp_path) -> None:
    run_dir, ref_path, _mov_path = _make_metadata_selection_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "ref_class1",
            "--mode",
            "metadata",
            "--metadata-domain",
            "ref",
            "--metadata-column",
            "rlnClassNumber",
            "--metadata-value",
            "1",
        ]
    )

    select_command(args)

    selection_dir = run_dir / "selections" / "ref_class1"
    rows = (selection_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    selected_rows = pd.read_csv(selection_dir / "selected_landscape_rows.csv")
    summary = json.loads((selection_dir / "selection_summary.json").read_text(encoding="utf-8"))
    payload = json.loads((selection_dir / "selection.json").read_text(encoding="utf-8"))
    assert rows == ["selection_rank,particle_key", "1,p1", "2,p3"]
    assert selected_rows["ref_source_row_id"].tolist() == [0, 2]
    assert summary["metadata_domain"] == "ref"
    assert summary["metadata_source_file"] == str(ref_path)
    assert summary["metadata_column"] == "rlnClassNumber"
    assert summary["metadata_values"] == ["1"]
    assert summary["metadata_candidate_count"] == 3
    assert summary["metadata_missing_count"] == 0
    assert payload["metric"] == "metadata_value_match"
    assert payload["metadata_source_row_id_field"] == "ref_source_row_id"


def test_select_metadata_multi_value_selects_union(tmp_path) -> None:
    run_dir, _ref_path, _mov_path = _make_metadata_selection_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "ref_classes_1_2",
            "--mode",
            "metadata",
            "--metadata-domain",
            "ref",
            "--metadata-column",
            "rlnClassNumber",
            "--metadata-value",
            "1,2",
        ]
    )

    select_command(args)

    rows = (
        run_dir / "selections" / "ref_classes_1_2" / "selected_particle_keys.csv"
    ).read_text(encoding="utf-8").splitlines()
    assert rows == ["selection_rank,particle_key", "1,p1", "2,p2", "3,p3"]


def test_select_metadata_group_splits_unique_matched_values(tmp_path) -> None:
    run_dir, _ref_path, _mov_path = _make_metadata_selection_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "ref_rlnClassNumber",
            "--mode",
            "metadata",
            "--metadata-domain",
            "ref",
            "--metadata-column",
            "rlnClassNumber",
            "--split-by-value",
        ]
    )

    select_command(args)

    class1_dir = run_dir / "selections" / "ref_rlnClassNumber_1"
    class2_dir = run_dir / "selections" / "ref_rlnClassNumber_2"
    assert (class1_dir / "selection.json").exists()
    assert (class1_dir / "selected_landscape_rows.csv").exists()
    assert (class2_dir / "selection.json").exists()
    class1_rows = (class1_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    class2_rows = (class2_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    assert class1_rows == ["selection_rank,particle_key", "1,p1", "2,p3"]
    assert class2_rows == ["selection_rank,particle_key", "1,p2"]


def test_select_metadata_split_by_value_conflicts_with_metadata_value(tmp_path) -> None:
    run_dir, _ref_path, _mov_path = _make_metadata_selection_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "bad_split",
            "--mode",
            "metadata",
            "--metadata-domain",
            "ref",
            "--metadata-column",
            "rlnClassNumber",
            "--metadata-value",
            "1",
            "--split-by-value",
        ]
    )

    with pytest.raises(ValueError, match="split-by-value"):
        select_command(args)


def test_select_metadata_requires_explicit_domain(tmp_path) -> None:
    run_dir, _ref_path, _mov_path = _make_metadata_selection_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "bad_metadata",
            "--mode",
            "metadata",
            "--metadata-column",
            "rlnClassNumber",
            "--metadata-value",
            "1",
        ]
    )

    with pytest.raises(ValueError, match="metadata-domain"):
        select_command(args)


def test_select_metadata_missing_column_fails_clearly(tmp_path) -> None:
    run_dir, _ref_path, _mov_path = _make_metadata_selection_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "missing_column",
            "--mode",
            "metadata",
            "--metadata-domain",
            "ref",
            "--metadata-column",
            "rlnNoSuchColumn",
            "--metadata-value",
            "1",
        ]
    )

    with pytest.raises(ValueError, match="missing column"):
        select_command(args)


def test_metadata_selection_export_still_uses_selected_row_provenance(tmp_path) -> None:
    run_dir, _ref_path, _mov_path = _make_metadata_selection_run_bundle(tmp_path)
    parser = build_parser()
    select_command(
        parser.parse_args(
            [
                "select",
                "--run-dir",
                str(run_dir),
                "--selection-id",
                "ref_class1",
                "--mode",
                "metadata",
                "--metadata-domain",
                "ref",
                "--metadata-column",
                "rlnClassNumber",
                "--metadata-value",
                "1",
            ]
        )
    )

    export_args = parser.parse_args(
        [
            "export",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "ref_class1",
            "--domain",
            "ref",
            "--format",
            "keys",
        ]
    )
    export_metadata_command(export_args)

    exported_rows = pd.read_csv(run_dir / "exports" / "ref_class1" / "selected_landscape_rows.csv")
    assert exported_rows["particle_key"].tolist() == ["p1", "p3"]
    assert exported_rows["ref_source_row_id"].tolist() == [0, 2]


def test_select_writes_inherited_selected_landscape(tmp_path) -> None:
    landscape = _make_landscape()
    landscape.data["sld_raw"] = [11.0, 22.0, 33.0]
    run_dir = _make_select_run_bundle(tmp_path, landscape)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--mode",
            "random",
            "--fraction",
            "0.67",
            "--seed",
            "7",
            "--write-selected-landscape",
        ]
    )

    select_command(args)

    selected_dir = output_dir / "selected_landscape"
    report = json.loads((selected_dir / "landscape_report.json").read_text(encoding="utf-8"))
    selected = read_landscape(selected_dir / "landscape.npz")
    assert (selected_dir / "landscape.csv").exists()
    assert (output_dir / "selected_landscape_rows.csv").exists()
    assert report["density_source"] == "parent_landscape"
    assert report["sld_recomputed"] is False
    assert len(selected.data) == report["selected_count"]
    assert "parent_sld_raw" not in selected.data.columns


def test_select_recomputed_selected_landscape_preserves_parent_sld(tmp_path) -> None:
    landscape = _make_landscape()
    landscape.data["sld_raw"] = [11.0, 22.0, 33.0]
    landscape.data["sld_display"] = [1.1, 2.2, 3.3]
    landscape.data["sld_display_is_outlier"] = [False, True, False]
    run_dir = _make_select_run_bundle(tmp_path, landscape)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--mode",
            "random",
            "--fraction",
            "0.67",
            "--seed",
            "7",
            "--write-selected-landscape",
            "--recompute-sld",
        ]
    )

    select_command(args)

    selected_dir = output_dir / "selected_landscape"
    report = json.loads((selected_dir / "landscape_report.json").read_text(encoding="utf-8"))
    selected = read_landscape(selected_dir / "landscape.npz")
    csv_table = pd.read_csv(selected_dir / "landscape.csv")
    assert report["density_source"] == "recomputed_on_selection"
    assert report["sld_recomputed"] is True
    assert report["effective_k_neighbors"] == len(selected.data) - 1
    assert "parent_sld_raw" in selected.data.columns
    assert "parent_sld_display" in selected.data.columns
    assert "parent_sld_display_is_outlier" in selected.data.columns
    assert "parent_sld_raw" in csv_table.columns
    assert not np.allclose(selected.data["sld_raw"], selected.data["parent_sld_raw"])


def test_select_range_by_coordinates_with_euler_alpha_bounds(tmp_path) -> None:
    landscape = _make_landscape()
    landscape.data["coordinates_analysis"] = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.0, 0.0, np.pi / 2]),
        np.array([0.0, 0.0, np.pi]),
    ]
    landscape.data["coordinates_display"] = landscape.data["coordinates_analysis"]
    run_dir = _make_select_run_bundle(tmp_path, landscape)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--mode",
            "range",
            "--range-bound",
            "alpha:60:120",
        ]
    )

    select_command(args)

    rows = (output_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    assert rows == ["selection_rank,particle_key", "1,p2"]


def test_select_radius_around_center_with_rotvec_center(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--center",
            "0",
            "0",
            "0",
            "--center-representation",
            "rotvec",
            "--radius",
            "30",
        ]
    )

    select_command(args)

    rows = (output_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    assert rows == ["selection_rank,particle_key", "1,p1"]
    selection_payload = json.loads((output_dir / "selection.json").read_text(encoding="utf-8"))
    assert selection_payload["metric"] == "so3_geodesic"
    assert selection_payload["radius_unit"] == "degrees"
    assert selection_payload["active_policy"]["metric"] == "so3_geodesic"
    assert selection_payload["active_policy"]["radius_unit"] == "degrees"
    summary = json.loads((output_dir / "selection_summary.json").read_text(encoding="utf-8"))
    assert summary["metric"] == "so3_geodesic"
    assert summary["radius_unit"] == "degrees"


def test_select_radius_deg_records_radius_unit_and_metric(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--center",
            "0",
            "0",
            "0",
            "--radius",
            "30",
        ]
    )

    select_command(args)

    selection_payload = json.loads((output_dir / "selection.json").read_text(encoding="utf-8"))
    assert selection_payload["metric"] == "so3_geodesic"
    assert selection_payload["radius_unit"] == "degrees"
    assert selection_payload["active_policy"]["metric"] == "so3_geodesic"
    assert selection_payload["active_policy"]["radius_unit"] == "degrees"
    assert selection_payload["radius"] == pytest.approx(30.0)


def test_select_radius_rad_records_radius_unit(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    output_dir = run_dir / "selections" / "selection"
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--center",
            "0",
            "0",
            "0",
            "--radius-rad",
            "0.5",
        ]
    )

    select_command(args)

    selection_payload = json.loads((output_dir / "selection.json").read_text(encoding="utf-8"))
    assert selection_payload["radius_unit"] == "radians"
    assert selection_payload["radius"] == pytest.approx(0.5)


def test_select_rejects_multiple_radius_unit_arguments(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--center",
            "0",
            "0",
            "0",
            "--radius",
            "30",
            "--radius-rad",
            "0.5",
        ]
    )

    with pytest.raises(ValueError, match="Use only one radius argument"):
        select_command(args)


def test_select_help_documents_compact_public_surface(capsys) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["select", "--help"])

    help_text = capsys.readouterr().out
    for option in (
        "--run-dir",
        "--selection-id",
        "--space",
        "--canonical-id",
        "--mode",
        "--center",
        "--radius",
        "--radius-rad",
        "--center-representation",
        "--metric",
        "--sld-min",
        "--sld-max",
        "--range-bound",
        "--fraction",
        "--seed",
        "--metadata-domain",
        "--metadata-column",
        "--metadata-value",
        "--split-by-value",
        "--write-selected-landscape",
        "--recompute-sld",
        "--overwrite",
    ):
        assert option in help_text
    assert "selected_landscape" in help_text
    assert "use-selected-landscape" in help_text
    assert "Recompute SLD only" in help_text
    assert "raw/canonical landscapes are unchanged" in help_text
    for removed in (
        "--output",
        "--landscape",
        "--radius-deg",
        "--center-input-space",
        "--center-euler-sequence",
        "--center-radians",
        "--euler-convention",
        "--range-coordinate-source",
        "--range-representation",
        "--range-euler-sequence",
        "--range-radians",
        "--density-support-field",
        "--density-artifact-policy",
        "--evaluation-space",
        "--random-fraction",
        "--random-seed",
        "--metadata-values",
        "--split-by-metadata",
        "--no-recompute-sld",
    ):
        assert removed not in help_text


@pytest.mark.parametrize(
    "removed_arg",
    [
        "--output",
        "--landscape",
        "--radius-deg",
        "--center-input-space",
        "--center-euler-sequence",
        "--center-radians",
        "--euler-convention",
        "--range-coordinate-source",
        "--range-representation",
        "--range-euler-sequence",
        "--range-radians",
        "--density-support-field",
        "--density-artifact-policy",
        "--evaluation-space",
        "--random-fraction",
        "--random-seed",
        "--metadata-values",
        "--split-by-metadata",
        "--no-recompute-sld",
    ],
)
def test_select_rejects_removed_public_arguments(removed_arg) -> None:
    parser = build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            ["select", "--run-dir", "run", "--selection-id", "selection", removed_arg, "value"]
        )


def test_select_canonical_space_fails_without_canonical_coordinates(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--space",
            "canonical",
            "--center",
            "0",
            "0",
            "0",
            "--radius",
            "1.0",
        ]
    )

    with pytest.raises(ValueError, match="canonical_landscape"):
        select_command(args)


def test_select_existing_output_directory_fails_by_default(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    output_dir = run_dir / "selections" / "selection"
    output_dir.mkdir(parents=True)
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--center",
            "0",
            "0",
            "0",
            "--radius",
            "30",
        ]
    )

    with pytest.raises(FileExistsError, match="Selection output directory already exists"):
        select_command(args)


def test_select_overwrite_permits_existing_output_directory(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    raw_landscape_path = run_dir / "data" / "raw_landscape.npz"
    raw_landscape_before = raw_landscape_path.read_bytes()
    other_selection_dir = run_dir / "selections" / "other_id"
    other_selection_dir.mkdir(parents=True)
    other_marker = other_selection_dir / "marker.txt"
    other_marker.write_text("untouched", encoding="utf-8")
    parser = build_parser()
    initial_args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--center",
            "0",
            "0",
            "0",
            "--radius",
            "30",
        ]
    )
    overwrite_args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--center",
            "0",
            "0",
            "0",
            "--radius",
            "30",
            "--overwrite",
        ]
    )

    assert select_command(initial_args) == 0
    assert select_command(overwrite_args) == 0

    output_dir = run_dir / "selections" / "selection"
    for artifact_name in (
        "selection.json",
        "selection.csv",
        "selected_particle_keys.csv",
        "selected_landscape_rows.csv",
        "selection_summary.json",
    ):
        assert (output_dir / artifact_name).exists()
    summary = json.loads((output_dir / "selection_summary.json").read_text(encoding="utf-8"))
    assert summary["overwrite"] is True
    assert raw_landscape_path.read_bytes() == raw_landscape_before
    assert other_marker.read_text(encoding="utf-8") == "untouched"


def test_select_does_not_modify_landscape_artifact(tmp_path) -> None:
    run_dir = _make_select_run_bundle(tmp_path)
    landscape_path = run_dir / "data" / "raw_landscape.npz"
    before = landscape_path.read_bytes()
    parser = build_parser()
    args = parser.parse_args(
        [
            "select",
            "--run-dir",
            str(run_dir),
            "--selection-id",
            "selection",
            "--center",
            "0",
            "0",
            "0",
            "--radius",
            "30",
        ]
    )

    select_command(args)

    assert landscape_path.read_bytes() == before


def test_export_selection_loads_existing_selection_json_and_writes_artifacts(tmp_path) -> None:
    source_dir = tmp_path / "source"
    output_dir = tmp_path / "reexport"
    source_report = export_selection(
        _make_selection(),
        policy=SelectionExportPolicy(output_dir=source_dir),
    )
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "selection",
            "--selection",
            source_report.output_paths["selection_json"],
            "--output-dir",
            str(output_dir),
        ]
    )

    assert export_selection_command(args) == 0

    assert (output_dir / "selection.json").exists()
    assert (output_dir / "selected_particle_keys.csv").exists()
    assert (output_dir / "export_report.json").exists()
    rows = (output_dir / "selected_particle_keys.csv").read_text(
        encoding="utf-8"
    ).splitlines()
    assert rows == ["selection_rank,particle_key", "1,p3", "2,p1", "3,p2"]


def test_export_selection_existing_output_directory_fails_by_default(tmp_path) -> None:
    source_dir = tmp_path / "source"
    output_dir = tmp_path / "reexport"
    output_dir.mkdir()
    source_report = export_selection(
        _make_selection(),
        policy=SelectionExportPolicy(output_dir=source_dir),
    )
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "selection",
            "--selection",
            source_report.output_paths["selection_json"],
            "--output-dir",
            str(output_dir),
        ]
    )

    with pytest.raises(FileExistsError, match="Selection export output directory already exists"):
        export_selection_command(args)


def test_export_selection_overwrite_permits_existing_output_directory(tmp_path) -> None:
    source_dir = tmp_path / "source"
    output_dir = tmp_path / "reexport"
    output_dir.mkdir()
    (output_dir / "selection.json").write_text("old\n", encoding="utf-8")
    source_report = export_selection(
        _make_selection(),
        policy=SelectionExportPolicy(output_dir=source_dir),
    )
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "selection",
            "--selection",
            source_report.output_paths["selection_json"],
            "--output-dir",
            str(output_dir),
            "--overwrite",
        ]
    )

    assert export_selection_command(args) == 0
    payload = json.loads((output_dir / "selection.json").read_text(encoding="utf-8"))
    assert payload["selection_id"] == "sel-1"


def test_export_selection_does_not_mutate_input_selection_json(tmp_path) -> None:
    source_dir = tmp_path / "source"
    output_dir = tmp_path / "reexport"
    source_report = export_selection(
        _make_selection(),
        policy=SelectionExportPolicy(output_dir=source_dir),
    )
    selection_path = source_dir / "selection.json"
    before = selection_path.read_text(encoding="utf-8")
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "selection",
            "--selection",
            source_report.output_paths["selection_json"],
            "--output-dir",
            str(output_dir),
        ]
    )

    export_selection_command(args)

    assert selection_path.read_text(encoding="utf-8") == before


def test_export_selection_does_not_load_landscape_or_recompute_selection(tmp_path, monkeypatch) -> None:
    source_dir = tmp_path / "source"
    output_dir = tmp_path / "reexport"
    source_report = export_selection(
        _make_selection(),
        policy=SelectionExportPolicy(output_dir=source_dir),
    )
    monkeypatch.setattr(
        "cryorole.cli.main.read_landscape_json",
        lambda *_args, **_kwargs: pytest.fail("export selection must not load a landscape"),
    )
    monkeypatch.setattr(
        "cryorole.cli.main.select_particles",
        lambda *_args, **_kwargs: pytest.fail("export selection must not recompute selection"),
    )
    parser = build_parser()
    args = parser.parse_args(
        [
            "export",
            "selection",
            "--selection",
            source_report.output_paths["selection_json"],
            "--output-dir",
            str(output_dir),
        ]
    )

    export_selection_command(args)

    assert (output_dir / "export_report.json").exists()


def test_manifest_command_writes_manifest(tmp_path) -> None:
    output_path = tmp_path / "manifest.json"
    parser = build_parser()
    args = parser.parse_args(["manifest", "--output", str(output_path)])

    assert manifest_command(args) == 0
    assert output_path.exists()
