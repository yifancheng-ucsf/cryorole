from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.transform import Rotation

from cryorole.cli.main import main as cryorole_main
from cryorole.export import (
    read_landscape,
    resolve_landscape_path,
    write_canonical_landscape_csv,
    write_raw_landscape_csv,
)
from cryorole.models.landscape import Landscape
from cryorole.workflows.rotate_landscape import main, rotate_landscape_bundle


def _landscape(*, canonical_transform: np.ndarray | None = None) -> Landscape:
    raw = np.array(
        [
            [0.10, 0.02, -0.03],
            [-0.15, 0.08, 0.04],
            [0.05, -0.12, 0.09],
            [0.18, 0.11, -0.07],
            [-0.06, -0.04, 0.16],
            [0.03, 0.17, 0.12],
            [-0.11, 0.06, -0.14],
            [0.14, -0.09, 0.05],
        ],
        dtype=float,
    )
    data = pd.DataFrame(
        {
            "particle_key": [f"p{index}" for index in range(len(raw))],
            "ref_source_row_id": np.arange(10, 10 + len(raw)),
            "mov_source_row_id": np.arange(30, 30 + len(raw)),
            "coordinates_analysis": list(raw),
            "coordinates_display": list(raw.copy()),
            "sld_unfloored": np.linspace(0.8, 1.5, len(raw)),
            "sld_raw": np.linspace(0.9, 1.6, len(raw)),
            "sld_display": np.linspace(0.9, 1.6, len(raw)),
            "sld_display_is_outlier": [False] * len(raw),
            "sld_was_floored": [False, True] * (len(raw) // 2),
            "sld_local_k_mean": np.linspace(0.10, 0.17, len(raw)),
            "sld_effective_local_k_mean": np.linspace(0.11, 0.18, len(raw)),
            "sld_distance_floor": np.full(len(raw), 1e-4),
        }
    )
    if canonical_transform is not None:
        data["coordinates_canonical"] = list(raw @ canonical_transform)
    return Landscape(data=data, canonical_transform=canonical_transform)


def _rotations_from(landscape: Landscape, column: str) -> np.ndarray:
    coordinates = np.vstack(landscape.data[column])
    return Rotation.from_rotvec(coordinates).as_matrix()


def test_raw_rotation_composes_on_right_and_writes_visualize_bundle(tmp_path) -> None:
    source = tmp_path / "source.csv"
    original = _landscape()
    write_raw_landscape_csv(original, source)
    source_before = source.read_bytes()
    output_dir = tmp_path / "rotated"

    report = rotate_landscape_bundle(
        source,
        output_dir,
        rotation_euler_deg=(20.0, 10.0, 0.0),
        space="raw",
    )

    assert source.read_bytes() == source_before
    source_landscape = read_landscape(source)
    assert resolve_landscape_path(output_dir, space="raw") == (
        output_dir / "data" / "raw_landscape.npz"
    )
    rotated = read_landscape(output_dir, space="raw")
    increment = Rotation.from_euler("zyx", [20.0, 10.0, 0.0], degrees=True).as_matrix()
    expected = _rotations_from(original, "coordinates_analysis") @ increment
    np.testing.assert_allclose(
        _rotations_from(rotated, "coordinates_analysis"),
        expected,
        atol=1e-12,
    )
    for field in (
        "sld_unfloored",
        "sld_raw",
        "sld_display",
        "sld_local_k_mean",
        "sld_effective_local_k_mean",
        "sld_distance_floor",
    ):
        np.testing.assert_array_equal(rotated.data[field], source_landscape.data[field])
    np.testing.assert_array_equal(
        rotated.data["sld_was_floored"],
        source_landscape.data["sld_was_floored"],
    )
    assert rotated.data["particle_key"].tolist() == original.data["particle_key"].tolist()
    assert rotated.data["ref_source_row_id"].tolist() == list(range(10, 18))
    assert report["operation"] == "moving_domain_local_right_rotation"
    assert report["ro_transform"] == "R_ro_shifted = R_ro @ G_raw"
    assert report["density_source"] == "parent_landscape"
    assert report["sld_recomputed"] is False
    assert report["rotation_euler_deg"] == {
        "alpha": 20.0,
        "beta": 10.0,
        "gamma": 0.0,
    }
    assert not (output_dir / "canonical" / "default" / "canonical_landscape.npz").exists()

    summary = json.loads((output_dir / "run_summary.json").read_text(encoding="utf-8"))
    assert summary["euler_convention"] == "extrinsic_zyx"
    assert summary["scipy_euler_sequence"] == "zyx"
    manifest = json.loads((output_dir / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["derived_bundle"] is True
    assert manifest["source_landscape_path"] == str(source.resolve())


def test_canonical_csv_rotation_is_mapped_back_to_consistent_raw_coordinates(
    tmp_path,
) -> None:
    canonical_transform = Rotation.from_euler(
        "xyz",
        [17.0, -23.0, 31.0],
        degrees=True,
    ).as_matrix()
    original = _landscape(canonical_transform=canonical_transform)
    source = tmp_path / "canonical_landscape.csv"
    write_canonical_landscape_csv(original, source)
    output_dir = tmp_path / "rotated"

    report = rotate_landscape_bundle(
        source,
        output_dir,
        rotation_euler_deg=(20.0, 10.0, 0.0),
        space="canonical",
    )

    raw_rotated = read_landscape(output_dir, space="raw")
    canonical_rotated = read_landscape(
        output_dir,
        space="canonical",
        canonical_id="default",
    )
    increment_canonical = Rotation.from_euler(
        "zyx",
        [20.0, 10.0, 0.0],
        degrees=True,
    ).as_matrix()
    increment_raw = (
        canonical_transform @ increment_canonical @ canonical_transform.T
    )
    expected_raw = (
        _rotations_from(original, "coordinates_analysis") @ increment_raw
    )
    expected_canonical = (
        _rotations_from(original, "coordinates_canonical") @ increment_canonical
    )
    np.testing.assert_allclose(
        _rotations_from(raw_rotated, "coordinates_analysis"),
        expected_raw,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        _rotations_from(canonical_rotated, "coordinates_canonical"),
        expected_canonical,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        canonical_rotated.canonical_transform,
        canonical_transform,
        atol=1e-12,
    )
    assert report["canonical_transform_source"] == "inferred_from_landscape_coordinates"
    assert resolve_landscape_path(
        output_dir,
        space="canonical",
        canonical_id="default",
    ) == output_dir / "canonical" / "default" / "canonical_landscape.npz"
    assert (output_dir / "canonical" / "default" / "canonical_frame.json").exists()


def test_raw_run_bundle_reuses_existing_canonical_frame_and_writes_both_spaces(
    tmp_path,
) -> None:
    canonical_transform = Rotation.from_euler(
        "xyz",
        [5.0, 15.0, -12.0],
        degrees=True,
    ).as_matrix()
    source_dir = tmp_path / "parent"
    (source_dir / "data").mkdir(parents=True)
    (source_dir / "canonical" / "default").mkdir(parents=True)
    landscape = _landscape(canonical_transform=canonical_transform)
    write_raw_landscape_csv(landscape, source_dir / "data" / "raw_landscape.csv")
    np.savez(
        source_dir / "canonical" / "default" / "canonical_frame.npz",
        canonical_transform=canonical_transform,
    )
    (source_dir / "run_summary.json").write_text(
        json.dumps(
            {
                "requested_sld_metric": "so3_geodesic",
                "resolved_sld_metric": "so3_geodesic",
            }
        ),
        encoding="utf-8",
    )

    output_dir = tmp_path / "rotated"
    report = rotate_landscape_bundle(
        source_dir,
        output_dir,
        rotation_euler_deg=(5.0, 0.0, 0.0),
        space="raw",
    )

    assert resolve_landscape_path(output_dir, space="raw").exists()
    assert resolve_landscape_path(
        output_dir,
        space="canonical",
        canonical_id="default",
    ).exists()
    assert report["canonical_transform_source"] == "parent_run_canonical_frame"
    assert report["resolved_sld_metric"] == "so3_geodesic"
    assert "geodesic distances are invariant" in report["warnings"][0]


def test_canonical_mode_requires_canonical_coordinates(tmp_path) -> None:
    source = tmp_path / "raw.csv"
    write_raw_landscape_csv(_landscape(), source)

    with pytest.raises(ValueError, match="canonical coordinates"):
        rotate_landscape_bundle(
            source,
            tmp_path / "rotated",
            rotation_euler_deg=(20.0, 10.0, 0.0),
            space="canonical",
        )


def test_existing_output_requires_overwrite(tmp_path) -> None:
    source = tmp_path / "source.csv"
    write_raw_landscape_csv(_landscape(), source)
    output_dir = tmp_path / "rotated"
    output_dir.mkdir()
    (output_dir / "keep.txt").write_text("user data", encoding="utf-8")
    stale_canonical = (
        output_dir / "canonical" / "default" / "canonical_landscape.npz"
    )
    stale_canonical.parent.mkdir(parents=True)
    stale_canonical.write_bytes(b"stale")

    with pytest.raises(FileExistsError, match="already exists"):
        rotate_landscape_bundle(
            source,
            output_dir,
            rotation_euler_deg=(1.0, 2.0, 3.0),
        )

    rotate_landscape_bundle(
        source,
        output_dir,
        rotation_euler_deg=(1.0, 2.0, 3.0),
        overwrite=True,
    )
    assert (output_dir / "keep.txt").read_text(encoding="utf-8") == "user data"
    assert not stale_canonical.exists()


def test_script_main_accepts_csv_and_prints_output_dir(tmp_path, capsys) -> None:
    source = tmp_path / "source.csv"
    write_raw_landscape_csv(_landscape(), source)
    output_dir = tmp_path / "rotated"

    result = main(
        [
            "--input",
            str(source),
            "--space",
            "raw",
            "--rotation-euler",
            "20",
            "10",
            "0",
            "--output-dir",
            str(output_dir),
        ]
    )

    assert result == 0
    assert capsys.readouterr().out.strip() == str(output_dir.resolve())
    assert (output_dir / "data" / "raw_landscape.csv").exists()


def test_output_bundle_runs_public_visualize_directly(tmp_path, capsys) -> None:
    source = tmp_path / "source.csv"
    write_raw_landscape_csv(_landscape(), source)
    output_dir = tmp_path / "rotated"
    rotate_landscape_bundle(
        source,
        output_dir,
        rotation_euler_deg=(20.0, 10.0, 0.0),
    )

    result = cryorole_main(
        [
            "visualize",
            "--run-dir",
            str(output_dir),
            "--space",
            "raw",
            "--representation",
            "rotvec",
            "--visual-id",
            "rotation_test",
        ]
    )

    assert result == 0
    visualization_dir = output_dir / "visualizations" / "raw" / "rotation_test"
    assert (visualization_dir / "visualization_report.json").exists()
    assert (visualization_dir / "rotvec_3view_projection.png").exists()
    assert str(visualization_dir) in capsys.readouterr().out
