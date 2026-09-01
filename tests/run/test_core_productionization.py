from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.transform import Rotation

from cryorole.cli.main import build_parser, run_command
from cryorole.io.writers.landscape_store import read_landscape_npz_arrays
import cryorole.workflows.run_service as run_service


def _write_cs(path: Path, uids: list[int], rotvecs: np.ndarray) -> None:
    dtype = np.dtype(
        [
            ("uid", "<u8"),
            ("alignments3D/pose", "<f8", (3,)),
            ("alignments3D/shift", "<f8", (2,)),
            ("blob/path", "S32"),
        ]
    )
    values = np.zeros(len(uids), dtype=dtype)
    values["uid"] = np.asarray(uids, dtype=np.uint64)
    values["alignments3D/pose"] = np.asarray(rotvecs, dtype=float)
    values["blob/path"] = b"particles.mrc"
    with path.open("wb") as handle:
        np.save(handle, values)


def _run_args(ref: Path, mov: Path, output: Path, *extra: str):
    return build_parser().parse_args(
        [
            "run",
            "--ref",
            str(ref),
            "--mov",
            str(mov),
            "--output-dir",
            str(output),
            "--k-neighbors",
            "2",
            "--no-visualize",
            *extra,
        ]
    )


def _inputs(tmp_path: Path, *, mov_order: list[int] | None = None) -> tuple[Path, Path]:
    ref = tmp_path / "ref.cs"
    mov = tmp_path / "mov.cs"
    uids = [100, 101, 102, 103, 104, 105]
    ref_rv = np.array(
        [
            [0.00, 0.00, 0.00],
            [0.02, 0.01, 0.00],
            [0.04, 0.00, 0.01],
            [0.01, 0.05, 0.00],
            [0.00, 0.03, 0.04],
            [0.05, 0.02, 0.03],
        ]
    )
    mov_rv = ref_rv + np.array(
        [
            [0.10, 0.00, 0.00],
            [0.00, 0.12, 0.00],
            [0.00, 0.00, 0.14],
            [0.08, 0.06, 0.00],
            [0.00, 0.09, 0.07],
            [0.06, 0.00, 0.11],
        ]
    )
    order = mov_order or list(range(len(uids)))
    _write_cs(ref, uids, ref_rv)
    _write_cs(mov, [uids[index] for index in order], mov_rv[order])
    return ref, mov


def test_auto_run_is_array_native_transactional_and_records_source_identity(tmp_path) -> None:
    ref, mov = _inputs(tmp_path, mov_order=[2, 0, 5, 1, 4, 3])
    output = tmp_path / "run"

    assert run_command(_run_args(ref, mov, output, "--raw-csv-chunk-size", "2")) == 0

    summary = json.loads((output / "run_summary.json").read_text(encoding="utf-8"))
    manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    assert summary["run_backend_resolved"] == "array_native"
    assert summary["raw_csv_backend"] == "array_native_chunked"
    assert summary["raw_csv_chunk_size"] == 2
    assert summary["density_backend"] == "array_native_batched"
    assert summary["bundle_state"] == "completed"
    assert manifest["bundle_transaction"]["state"] == "completed"
    assert manifest["run_id"] == summary["run_id"]
    assert (output / ".cryorole_bundle_complete").is_file()
    for domain, source in (("ref", ref), ("mov", mov)):
        identity = summary["source_identities"][domain]
        assert identity["original_path"] == str(source)
        assert identity["resolved_path"] == str(source.resolve())
        assert identity["size_bytes"] == source.stat().st_size
        assert len(identity["sha256"]) == 64
        assert identity["row_count"] == 6

    arrays = read_landscape_npz_arrays(output / "data/raw_landscape.npz")
    assert arrays.n_points == 6
    assert arrays.ref_source_row_id.tolist() == [0, 1, 2, 3, 4, 5]
    assert arrays.mov_source_row_id.tolist() == [1, 3, 0, 5, 4, 2]
    with (output / "data/raw_landscape.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == arrays.n_points


def test_array_native_run_writes_complete_quicklook_and_report(tmp_path) -> None:
    ref, mov = _inputs(tmp_path)
    output = tmp_path / "run"
    args = build_parser().parse_args(
        [
            "run",
            "--ref",
            str(ref),
            "--mov",
            str(mov),
            "--output-dir",
            str(output),
            "--k-neighbors",
            "2",
        ]
    )

    assert run_command(args) == 0

    quicklook = output / "visualizations" / "quicklook"
    assert {path.name for path in quicklook.iterdir()} == {
        "all_euler_3view_projection.png",
        "all_rotvec_3view_projection.png",
        "sld_ge_1_euler_3view_projection.png",
        "sld_ge_1_rotvec_3view_projection.png",
        "top_40pct_euler_3view_projection.png",
        "top_40pct_rotvec_3view_projection.png",
        "sld_log_distribution.png",
    }
    summary = json.loads((output / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["raw_visualization_sld_distribution"]["input_count"] == 6
    assert summary["raw_visualization_subset_counts"]["all"] == 6
    assert (output / "run_report.md").is_file()


def test_array_native_and_dataframe_reference_backends_are_equivalent(tmp_path) -> None:
    ref, mov = _inputs(tmp_path)
    array_dir = tmp_path / "array"
    compat_dir = tmp_path / "compat"

    run_command(_run_args(ref, mov, array_dir, "--run-backend", "array_native"))
    run_command(_run_args(ref, mov, compat_dir, "--run-backend", "dataframe_compat"))
    array = read_landscape_npz_arrays(array_dir / "data/raw_landscape.npz")
    compat = read_landscape_npz_arrays(compat_dir / "data/raw_landscape.npz")

    assert array.particle_key.tolist() == compat.particle_key.tolist()
    assert np.array_equal(array.ref_source_row_id, compat.ref_source_row_id)
    assert np.array_equal(array.mov_source_row_id, compat.mov_source_row_id)
    assert np.allclose(array.coordinates_analysis, compat.coordinates_analysis, atol=1e-12)
    assert np.allclose(
        Rotation.from_rotvec(array.coordinates_analysis).as_matrix(),
        Rotation.from_rotvec(compat.coordinates_analysis).as_matrix(),
        atol=1e-12,
    )
    assert np.allclose(array.sld_raw, compat.sld_raw, atol=1e-12)
    assert np.array_equal(array.sld_was_floored, compat.sld_was_floored)
    array_csv = pd.read_csv(array_dir / "data/raw_landscape.csv")
    compat_csv = pd.read_csv(compat_dir / "data/raw_landscape.csv")
    euler_columns = [
        "raw_ea_zyx_alpha_deg",
        "raw_ea_zyx_beta_deg",
        "raw_ea_zyx_gamma_deg",
    ]
    assert np.allclose(array_csv[euler_columns], compat_csv[euler_columns], atol=1e-10)
    assert np.allclose(array_csv["raw_angle_deg"], compat_csv["raw_angle_deg"], atol=1e-10)

    array_summary = json.loads((array_dir / "run_summary.json").read_text(encoding="utf-8"))
    compat_summary = json.loads((compat_dir / "run_summary.json").read_text(encoding="utf-8"))
    for key in (
        "matched_count",
        "match_key",
        "ref_coverage",
        "mov_coverage",
        "overlap_smaller_input",
        "low_overlap_allowed",
        "requested_sld_metric",
        "resolved_sld_metric",
        "euler_convention",
        "scipy_euler_sequence",
    ):
        assert array_summary[key] == compat_summary[key]

    array_manifest = json.loads((array_dir / "run_manifest.json").read_text(encoding="utf-8"))
    compat_manifest = json.loads((compat_dir / "run_manifest.json").read_text(encoding="utf-8"))
    assert array_manifest["schema_version"] == compat_manifest["schema_version"] == "3.0"
    for policy_name in ("match_policy", "density_policy", "representation_policy"):
        assert (
            array_manifest["active_policies"][policy_name]
            == compat_manifest["active_policies"][policy_name]
        )
    assert (
        array_manifest["reports"]["match_report"]
        == compat_manifest["reports"]["match_report"]
    )


def test_low_overlap_requires_explicit_allow_and_is_recorded(tmp_path) -> None:
    ref, _ = _inputs(tmp_path)
    mov = tmp_path / "low_mov.cs"
    rotvecs = np.array(
        [
            [0.15, 0.00, 0.00],
            [0.00, 0.18, 0.00],
            [0.01, 0.02, 0.03],
            [0.02, 0.03, 0.04],
            [0.03, 0.04, 0.05],
            [0.04, 0.05, 0.06],
        ]
    )
    _write_cs(mov, [100, 101, 200, 201, 202, 203], rotvecs)

    with pytest.raises(ValueError, match="below MatchPolicy threshold"):
        run_command(_run_args(ref, mov, tmp_path / "failed"))
    assert not (tmp_path / "failed").exists()

    output = tmp_path / "allowed"
    run_command(_run_args(ref, mov, output, "--allow-low-overlap"))
    report = json.loads((output / "reports/match_report.json").read_text(encoding="utf-8"))["report"]
    summary = json.loads((output / "run_summary.json").read_text(encoding="utf-8"))
    assert report["matched_count"] == 2
    assert report["ref_coverage"] == pytest.approx(2 / 6)
    assert report["mov_coverage"] == pytest.approx(2 / 6)
    assert report["overlap_smaller_input"] == pytest.approx(2 / 6)
    assert report["low_overlap_allowed"] is True
    assert summary["low_overlap_allowed"] is True


def test_zero_matches_and_duplicate_uid_fail_before_pose_normalization(tmp_path) -> None:
    ref, _ = _inputs(tmp_path)
    zero = tmp_path / "zero.cs"
    malformed = np.zeros((6, 3), dtype=float)
    _write_cs(zero, [200, 201, 202, 203, 204, 205], malformed)
    with pytest.raises(ValueError, match="zero matches"):
        run_command(_run_args(ref, zero, tmp_path / "zero-run"))

    duplicate = tmp_path / "duplicate.cs"
    _write_cs(duplicate, [100, 100, 102, 103, 104, 105], malformed)
    with pytest.raises(ValueError, match="failed_duplicates|duplicate"):
        run_command(_run_args(ref, duplicate, tmp_path / "duplicate-run"))


@pytest.mark.parametrize(
    ("symbol", "expected_stage"),
    [
        ("preflight_and_normalize_matched_arrays", "matching"),
        ("compute_relative_orientation_arrays", "computing_ro"),
        ("compute_landscape_density_arrays", "computing_sld"),
        ("write_landscape_npz_arrays", "writing"),
        ("write_raw_landscape_csv_from_arrays", "writing"),
        ("write_report_json", "writing"),
    ],
)
def test_run_fault_injection_never_publishes_partial_bundle(
    tmp_path,
    monkeypatch,
    symbol,
    expected_stage,
) -> None:
    ref, mov = _inputs(tmp_path)
    output = tmp_path / f"run-{symbol}"

    def injected_failure(*args, **kwargs):
        raise RuntimeError(f"injected {symbol} failure")

    monkeypatch.setattr(run_service, symbol, injected_failure)
    with pytest.raises(RuntimeError, match=f"injected {symbol} failure"):
        run_command(_run_args(ref, mov, output))

    assert not output.exists()
    failed = list(tmp_path.glob(f".{output.name}.failed-*"))
    assert len(failed) == 1
    report = json.loads((failed[0] / "failure_report.json").read_text(encoding="utf-8"))
    assert report["failed_stage"] == expected_stage
