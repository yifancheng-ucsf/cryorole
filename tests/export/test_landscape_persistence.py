from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from cryorole.export import (
    read_landscape,
    read_landscape_csv,
    read_landscape_json,
    read_landscape_npz,
    read_landscape_npz_arrays,
    write_canonical_landscape_csv_from_arrays,
    write_landscape_json,
    write_landscape_npz,
    write_landscape_npz_arrays,
    write_raw_landscape_csv,
)
from cryorole.models.landscape import Landscape
from cryorole.models.landscape_arrays import LandscapeArrays


def _make_landscape() -> Landscape:
    return Landscape(
        data=pd.DataFrame(
            {
                "particle_key": ["p1", "p2"],
                "coordinates_analysis": [
                    np.array([0.0, 0.0, 0.0]),
                    np.array([1.0, 0.0, 0.0]),
                ],
                "coordinates_display": [
                    np.array([0.0, 0.0, 0.0]),
                    np.array([1.0, 0.0, 0.0]),
                ],
                "sld_unfloored": [1.0, 1.0],
                "sld_raw": [1.0, 1.0],
                "sld_display": [1.0, 1.0],
                "sld_display_is_outlier": [False, False],
                "sld_was_floored": [False, False],
                "sld_local_k_mean": [1.0, 1.0],
                "sld_effective_local_k_mean": [1.0, 1.0],
                "sld_distance_floor": [0.0001, 0.0001],
            }
        )
    )


def test_read_landscape_json_rejects_unsupported_schema_version(tmp_path) -> None:
    path = tmp_path / "landscape.json"
    write_landscape_json(_make_landscape(), path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["schema_version"] = "2"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="schema_version must be '1'"):
        read_landscape_json(path)


def test_read_landscape_json_rejects_mismatched_row_count(tmp_path) -> None:
    path = tmp_path / "landscape.json"
    write_landscape_json(_make_landscape(), path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["row_count"] = 99
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="row_count does not match data length"):
        read_landscape_json(path)


def test_read_landscape_json_uses_landscape_validation_for_missing_columns(tmp_path) -> None:
    path = tmp_path / "landscape.json"
    write_landscape_json(_make_landscape(), path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    for row in payload["data"]:
        row.pop("sld_raw")
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="Landscape missing required columns"):
        read_landscape_json(path)


def test_landscape_npz_roundtrip_preserves_core_arrays(tmp_path) -> None:
    path = tmp_path / "raw_landscape.npz"
    landscape = _make_landscape()

    write_landscape_npz(landscape, path)
    loaded = read_landscape_npz(path)

    assert loaded.data["particle_key"].tolist() == ["p1", "p2"]
    np.testing.assert_allclose(
        np.vstack(loaded.data["coordinates_analysis"]),
        np.vstack(landscape.data["coordinates_analysis"]),
    )
    assert loaded.data["sld_raw"].tolist() == [1.0, 1.0]
    assert loaded.data["sld_display_is_outlier"].tolist() == [False, False]


def test_landscape_npz_arrays_roundtrip_preserves_ndarray_coordinates(tmp_path) -> None:
    path = tmp_path / "raw_landscape.npz"
    arrays = LandscapeArrays(
        particle_key=np.array(["p1", "p2"]),
        coordinates_analysis=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        coordinates_display=np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        sld_unfloored=np.ones(2),
        sld_raw=np.ones(2),
        sld_display=np.ones(2),
        sld_display_is_outlier=np.array([False, True]),
        sld_was_floored=np.array([False, False]),
        sld_local_k_mean=np.ones(2),
        sld_effective_local_k_mean=np.ones(2),
        sld_distance_floor=np.full(2, 1e-4),
        ref_source_row_id=np.array([10, 11]),
        mov_source_row_id=np.array([20, 21]),
    )

    write_landscape_npz_arrays(arrays, path)
    loaded = read_landscape_npz_arrays(path)

    assert loaded.coordinates_analysis.shape == (2, 3)
    assert isinstance(loaded.coordinates_analysis, np.ndarray)
    assert loaded.particle_key.tolist() == ["p1", "p2"]
    np.testing.assert_allclose(loaded.coordinates_analysis, arrays.coordinates_analysis)
    np.testing.assert_array_equal(loaded.ref_source_row_id, np.array([10, 11]))
    np.testing.assert_array_equal(loaded.sld_display_is_outlier, np.array([False, True]))


def test_chunked_canonical_csv_writer_writes_expected_columns_and_rows(tmp_path) -> None:
    path = tmp_path / "canonical_landscape.csv"
    arrays = LandscapeArrays(
        particle_key=np.array(["p1", "p2", "p3"]),
        coordinates_analysis=np.array(
            [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [0.0, 0.3, 0.0]]
        ),
        coordinates_display=np.array(
            [[0.0, 0.0, 0.0], [0.2, 0.0, 0.0], [0.0, 0.3, 0.0]]
        ),
        coordinates_canonical=np.array(
            [[0.0, 0.0, 0.0], [0.0, 0.2, 0.0], [0.3, 0.0, 0.0]]
        ),
        sld_unfloored=np.ones(3),
        sld_raw=np.ones(3),
        sld_display=np.ones(3),
        sld_display_is_outlier=np.array([False, False, False]),
        sld_was_floored=np.array([False, False, False]),
        sld_local_k_mean=np.ones(3),
        sld_effective_local_k_mean=np.ones(3),
        sld_distance_floor=np.full(3, 1e-4),
        ref_source_row_id=np.array([1, 2, 3]),
        mov_source_row_id=np.array([4, 5, 6]),
        canonical_transform=np.eye(3),
    )

    write_canonical_landscape_csv_from_arrays(arrays, path, chunk_size=2)

    table = pd.read_csv(path)
    assert len(table) == 3
    assert {
        "particle_key",
        "raw_rv_x_rad",
        "canonical_rv_x_rad",
        "canonical_ea_zyx_alpha_deg",
        "sld_raw",
        "coordinate_source",
    }.issubset(table.columns)
    assert table["particle_key"].tolist() == ["p1", "p2", "p3"]
    assert table["coordinate_source"].unique().tolist() == ["canonical"]


def test_landscape_csv_roundtrip_supports_unit_bearing_columns(tmp_path) -> None:
    path = tmp_path / "raw_landscape.csv"
    landscape = _make_landscape()

    write_raw_landscape_csv(landscape, path)
    loaded = read_landscape_csv(path)

    table = pd.read_csv(path)
    assert {"raw_rv_x_rad", "raw_ea_zyx_alpha_deg", "raw_angle_deg"}.issubset(
        table.columns
    )
    assert "sld_display_is_outlier" in table.columns
    assert loaded.data["particle_key"].tolist() == ["p1", "p2"]
    np.testing.assert_allclose(
        np.vstack(loaded.data["coordinates_analysis"]),
        np.vstack(landscape.data["coordinates_analysis"]),
    )


def test_raw_csv_keeps_zyx_column_names_for_intrinsic_compatibility(tmp_path) -> None:
    path = tmp_path / "raw_landscape.csv"

    write_raw_landscape_csv(_make_landscape(), path, euler_sequence="ZYX")

    table = pd.read_csv(path)
    assert {
        "raw_ea_zyx_alpha_deg",
        "raw_ea_zyx_beta_deg",
        "raw_ea_zyx_gamma_deg",
    }.issubset(table.columns)
    assert not any("intrinsic" in column or "extrinsic" in column for column in table.columns)


def test_read_legacy_landscape_csv_without_display_outlier_flag_defaults_false(tmp_path) -> None:
    path = tmp_path / "legacy_raw_landscape.csv"
    write_raw_landscape_csv(_make_landscape(), path)
    table = pd.read_csv(path)
    table = table.drop(columns=["sld_display_is_outlier"])
    table.to_csv(path, index=False)

    loaded = read_landscape_csv(path)

    assert loaded.data["sld_display_is_outlier"].tolist() == [False, False]


def test_read_legacy_landscape_csv_maps_display_clip_flag_to_outlier_flag(tmp_path) -> None:
    path = tmp_path / "legacy_raw_landscape.csv"
    write_raw_landscape_csv(_make_landscape(), path)
    table = pd.read_csv(path)
    table = table.drop(columns=["sld_display_is_outlier"])
    table["sld_display_was_clipped"] = [False, True]
    table.to_csv(path, index=False)

    loaded = read_landscape_csv(path)

    assert loaded.data["sld_display_is_outlier"].tolist() == [False, True]


def test_read_legacy_landscape_npz_without_display_outlier_flag_defaults_false(tmp_path) -> None:
    path = tmp_path / "legacy_raw_landscape.npz"
    write_landscape_npz(_make_landscape(), path)
    with np.load(path, allow_pickle=False) as payload:
        arrays = {
            key: np.asarray(payload[key])
            for key in payload.files
            if key != "sld_display_is_outlier"
        }
    np.savez_compressed(path, **arrays)

    loaded = read_landscape_npz(path)

    assert loaded.data["sld_display_is_outlier"].tolist() == [False, False]


def test_read_landscape_resolves_raw_npz_from_run_dir(tmp_path) -> None:
    run_dir = tmp_path / "run"
    data_dir = run_dir / "data"
    data_dir.mkdir(parents=True)
    write_landscape_npz(_make_landscape(), data_dir / "raw_landscape.npz")

    loaded = read_landscape(run_dir, space="raw")

    assert loaded.data["particle_key"].tolist() == ["p1", "p2"]
