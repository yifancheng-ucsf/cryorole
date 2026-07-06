"""Production landscape persistence for run-bundle artifacts."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from cryorole.core.euler_conventions import resolve_euler_convention
from cryorole.export.landscape import (
    read_landscape_json,
    read_landscape_json_metadata,
    write_canonical_landscape_csv,
    write_raw_landscape_csv,
)
from cryorole.models.landscape_arrays import LandscapeArrays
from cryorole.models.landscape import Landscape


LANDSCAPE_NPZ_SCHEMA_VERSION = "1"


def write_landscape_npz(
    landscape: Landscape,
    path: str | Path,
    *,
    overwrite: bool = False,
    artifact_type: str = "raw_landscape",
) -> Path:
    """Write compact machine-readable landscape arrays."""

    output_path = Path(path)
    _prepare_output_path(output_path, overwrite=overwrite)
    data = landscape.data
    arrays: dict[str, Any] = {
        "artifact_type": np.asarray(artifact_type),
        "schema_version": np.asarray(LANDSCAPE_NPZ_SCHEMA_VERSION),
        "particle_key": np.asarray(data["particle_key"].map(str), dtype=str),
        "coordinates_analysis": _stack_coordinates(
            data["coordinates_analysis"],
            "coordinates_analysis",
        ),
        "coordinates_display": _stack_coordinates(
            data["coordinates_display"],
            "coordinates_display",
        ),
        "sld_unfloored": np.asarray(data["sld_unfloored"], dtype=float),
        "sld_raw": np.asarray(data["sld_raw"], dtype=float),
        "sld_display": np.asarray(data["sld_display"], dtype=float),
        "sld_display_is_outlier": np.asarray(
            data["sld_display_is_outlier"],
            dtype=bool,
        ),
        "sld_was_floored": np.asarray(data["sld_was_floored"], dtype=bool),
        "sld_local_k_mean": np.asarray(data["sld_local_k_mean"], dtype=float),
        "sld_effective_local_k_mean": np.asarray(
            data["sld_effective_local_k_mean"],
            dtype=float,
        ),
        "sld_distance_floor": np.asarray(data["sld_distance_floor"], dtype=float),
        "ref_source_row_id": _integer_column_or_default(data, "ref_source_row_id"),
        "mov_source_row_id": _integer_column_or_default(data, "mov_source_row_id"),
    }
    if "coordinates_canonical" in data.columns:
        arrays["coordinates_canonical"] = _stack_coordinates(
            data["coordinates_canonical"],
            "coordinates_canonical",
        )
    for column in (
        "parent_sld_unfloored",
        "parent_sld_raw",
        "parent_sld_display",
        "parent_sld_display_is_outlier",
        "parent_sld_was_floored",
        "parent_sld_local_k_mean",
        "parent_sld_effective_local_k_mean",
        "parent_sld_distance_floor",
    ):
        if column in data.columns:
            dtype = bool if column.endswith("_is_outlier") or column.endswith("_was_floored") else float
            arrays[column] = np.asarray(data[column], dtype=dtype)
    if landscape.canonical_transform is not None:
        arrays["canonical_transform"] = np.asarray(
            landscape.canonical_transform,
            dtype=float,
        )
    np.savez_compressed(output_path, **arrays)
    return output_path


def write_landscape_npz_arrays(
    arrays: LandscapeArrays,
    path: str | Path,
    *,
    overwrite: bool = False,
    artifact_type: str = "raw_landscape",
) -> Path:
    """Write compact machine-readable landscape arrays without DataFrame inflation."""

    output_path = Path(path)
    _prepare_output_path(output_path, overwrite=overwrite)
    payload: dict[str, Any] = {
        "artifact_type": np.asarray(artifact_type),
        "schema_version": np.asarray(LANDSCAPE_NPZ_SCHEMA_VERSION),
        "particle_key": np.asarray(arrays.particle_key, dtype=str),
        "coordinates_analysis": np.asarray(arrays.coordinates_analysis, dtype=float),
        "sld_unfloored": np.asarray(arrays.sld_unfloored, dtype=float),
        "sld_raw": np.asarray(arrays.sld_raw, dtype=float),
        "sld_display": np.asarray(arrays.sld_display, dtype=float),
        "sld_display_is_outlier": np.asarray(arrays.sld_display_is_outlier, dtype=bool),
        "sld_was_floored": np.asarray(arrays.sld_was_floored, dtype=bool),
        "sld_local_k_mean": np.asarray(arrays.sld_local_k_mean, dtype=float),
        "sld_effective_local_k_mean": np.asarray(
            arrays.sld_effective_local_k_mean,
            dtype=float,
        ),
        "sld_distance_floor": np.asarray(arrays.sld_distance_floor, dtype=float),
    }
    if arrays.coordinates_display is not None:
        payload["coordinates_display"] = np.asarray(arrays.coordinates_display, dtype=float)
    if arrays.ref_source_row_id is not None:
        payload["ref_source_row_id"] = np.asarray(arrays.ref_source_row_id, dtype=np.int64)
    if arrays.mov_source_row_id is not None:
        payload["mov_source_row_id"] = np.asarray(arrays.mov_source_row_id, dtype=np.int64)
    if arrays.coordinates_canonical is not None:
        payload["coordinates_canonical"] = np.asarray(arrays.coordinates_canonical, dtype=float)
    if arrays.canonical_transform is not None:
        payload["canonical_transform"] = np.asarray(arrays.canonical_transform, dtype=float)
    np.savez(output_path, **payload)
    return output_path


def read_landscape_npz(path: str | Path) -> Landscape:
    """Read a compact NPZ landscape artifact into a Landscape."""

    input_path = Path(path)
    with np.load(input_path, allow_pickle=False) as payload:
        schema_version = str(payload["schema_version"].item())
        if schema_version != LANDSCAPE_NPZ_SCHEMA_VERSION:
            raise ValueError(
                "Landscape NPZ schema_version must be "
                f"{LANDSCAPE_NPZ_SCHEMA_VERSION!r}, got {schema_version!r}"
            )
        data = pd.DataFrame(
            {
                "particle_key": payload["particle_key"].astype(str),
                "coordinates_analysis": list(
                    np.asarray(payload["coordinates_analysis"], dtype=float)
                ),
                "coordinates_display": list(
                    np.asarray(payload["coordinates_display"], dtype=float)
                ),
                "sld_unfloored": np.asarray(payload["sld_unfloored"], dtype=float),
                "sld_raw": np.asarray(payload["sld_raw"], dtype=float),
                "sld_display": np.asarray(payload["sld_display"], dtype=float),
                "sld_display_is_outlier": _display_outlier_from_npz_payload(payload),
                "sld_was_floored": np.asarray(payload["sld_was_floored"], dtype=bool),
                "sld_local_k_mean": np.asarray(payload["sld_local_k_mean"], dtype=float),
                "sld_effective_local_k_mean": np.asarray(
                    payload["sld_effective_local_k_mean"],
                    dtype=float,
                ),
                "sld_distance_floor": np.asarray(payload["sld_distance_floor"], dtype=float),
            }
        )
        if "ref_source_row_id" in payload:
            data["ref_source_row_id"] = np.asarray(payload["ref_source_row_id"], dtype=np.int64)
        if "mov_source_row_id" in payload:
            data["mov_source_row_id"] = np.asarray(payload["mov_source_row_id"], dtype=np.int64)
        if "coordinates_canonical" in payload:
            data["coordinates_canonical"] = list(
                np.asarray(payload["coordinates_canonical"], dtype=float)
            )
        for column in (
            "parent_sld_unfloored",
            "parent_sld_raw",
            "parent_sld_display",
            "parent_sld_display_is_outlier",
            "parent_sld_was_floored",
            "parent_sld_local_k_mean",
            "parent_sld_effective_local_k_mean",
            "parent_sld_distance_floor",
        ):
            if column in payload:
                dtype = bool if column.endswith("_is_outlier") or column.endswith("_was_floored") else float
                data[column] = np.asarray(payload[column], dtype=dtype)
        canonical_transform = (
            np.asarray(payload["canonical_transform"], dtype=float)
            if "canonical_transform" in payload
            else None
        )
    return Landscape(data=data, canonical_transform=canonical_transform)


def read_landscape_npz_arrays(path: str | Path) -> LandscapeArrays:
    """Read an NPZ landscape artifact as compact arrays."""

    input_path = Path(path)
    with np.load(input_path, allow_pickle=False) as payload:
        schema_version = str(payload["schema_version"].item())
        if schema_version != LANDSCAPE_NPZ_SCHEMA_VERSION:
            raise ValueError(
                "Landscape NPZ schema_version must be "
                f"{LANDSCAPE_NPZ_SCHEMA_VERSION!r}, got {schema_version!r}"
            )
        return LandscapeArrays(
            particle_key=payload["particle_key"].astype(str),
            coordinates_analysis=np.asarray(payload["coordinates_analysis"], dtype=float),
            coordinates_display=(
                np.asarray(payload["coordinates_display"], dtype=float)
                if "coordinates_display" in payload
                else None
            ),
            sld_unfloored=np.asarray(payload["sld_unfloored"], dtype=float),
            sld_raw=np.asarray(payload["sld_raw"], dtype=float),
            sld_display=np.asarray(payload["sld_display"], dtype=float),
            sld_display_is_outlier=_display_outlier_from_npz_payload(payload),
            sld_was_floored=np.asarray(payload["sld_was_floored"], dtype=bool),
            sld_local_k_mean=np.asarray(payload["sld_local_k_mean"], dtype=float),
            sld_effective_local_k_mean=np.asarray(
                payload["sld_effective_local_k_mean"],
                dtype=float,
            ),
            sld_distance_floor=np.asarray(payload["sld_distance_floor"], dtype=float),
            ref_source_row_id=(
                np.asarray(payload["ref_source_row_id"], dtype=np.int64)
                if "ref_source_row_id" in payload
                else None
            ),
            mov_source_row_id=(
                np.asarray(payload["mov_source_row_id"], dtype=np.int64)
                if "mov_source_row_id" in payload
                else None
            ),
            coordinates_canonical=(
                np.asarray(payload["coordinates_canonical"], dtype=float)
                if "coordinates_canonical" in payload
                else None
            ),
            canonical_transform=(
                np.asarray(payload["canonical_transform"], dtype=float)
                if "canonical_transform" in payload
                else None
            ),
        )


def write_canonical_landscape_csv_from_npz(
    npz_path: str | Path,
    csv_path: str | Path,
    *,
    overwrite: bool = False,
    chunk_size: int = 100000,
    euler_sequence: str | None = None,
    euler_degrees: bool = True,
) -> Path:
    """Write canonical CSV from NPZ arrays in chunks."""

    arrays = read_landscape_npz_arrays(npz_path)
    return write_canonical_landscape_csv_from_arrays(
        arrays,
        csv_path,
        overwrite=overwrite,
        chunk_size=chunk_size,
        euler_sequence=euler_sequence,
        euler_degrees=euler_degrees,
    )


def write_canonical_landscape_csv_from_arrays(
    arrays: LandscapeArrays,
    path: str | Path,
    *,
    overwrite: bool = False,
    chunk_size: int = 100000,
    euler_sequence: str | None = None,
    euler_degrees: bool = True,
) -> Path:
    """Write a complete canonical landscape CSV without full object-column DataFrames."""

    if arrays.coordinates_canonical is None:
        raise ValueError("canonical_landscape.csv requires coordinates_canonical")
    if chunk_size <= 0:
        raise ValueError("csv chunk_size must be positive")
    output_path = Path(path)
    _prepare_output_path(output_path, overwrite=overwrite)
    header = _canonical_csv_header(euler_sequence)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for start in range(0, arrays.n_points, chunk_size):
            stop = min(start + chunk_size, arrays.n_points)
            _write_canonical_csv_chunk(
                writer,
                arrays,
                start=start,
                stop=stop,
                euler_sequence=euler_sequence,
                euler_degrees=euler_degrees,
            )
    return output_path


def read_landscape_csv(path: str | Path) -> Landscape:
    """Read a user-facing raw or canonical landscape CSV."""

    table = pd.read_csv(path)
    raw_coordinates = _coordinates_from_columns(
        table,
        preferred=("raw_rv_x_rad", "raw_rv_y_rad", "raw_rv_z_rad"),
        aliases=("rotvec_x", "rotvec_y", "rotvec_z"),
        canonical_aliases=("raw_rotvec_x", "raw_rotvec_y", "raw_rotvec_z"),
    )
    data = pd.DataFrame(
        {
            "particle_key": table["particle_key"].astype(str),
            "coordinates_analysis": list(raw_coordinates),
            "coordinates_display": list(raw_coordinates),
            "sld_unfloored": table["sld_unfloored"],
            "sld_raw": table["sld_raw"],
            "sld_display": table["sld_display"],
            "sld_display_is_outlier": _display_outlier_from_csv_table(table),
            "sld_was_floored": table["sld_was_floored"].map(_parse_bool),
            "sld_local_k_mean": table["sld_local_k_mean"],
            "sld_effective_local_k_mean": table["sld_effective_local_k_mean"],
            "sld_distance_floor": table["sld_distance_floor"],
        }
    )
    for column in ("ref_source_row_id", "mov_source_row_id"):
        if column in table.columns:
            data[column] = table[column]
    if {"canonical_rv_x_rad", "canonical_rv_y_rad", "canonical_rv_z_rad"}.issubset(
        table.columns
    ):
        data["coordinates_canonical"] = list(
            table[["canonical_rv_x_rad", "canonical_rv_y_rad", "canonical_rv_z_rad"]].to_numpy(
                dtype=float
            )
        )
    elif {"canonical_rotvec_x", "canonical_rotvec_y", "canonical_rotvec_z"}.issubset(
        table.columns
    ):
        data["coordinates_canonical"] = list(
            table[["canonical_rotvec_x", "canonical_rotvec_y", "canonical_rotvec_z"]].to_numpy(
                dtype=float
            )
        )
    for column in (
        "parent_sld_unfloored",
        "parent_sld_raw",
        "parent_sld_display",
        "parent_sld_display_is_outlier",
        "parent_sld_was_floored",
        "parent_sld_local_k_mean",
        "parent_sld_effective_local_k_mean",
        "parent_sld_distance_floor",
    ):
        if column in table.columns:
            if column.endswith("_is_outlier") or column.endswith("_was_floored"):
                data[column] = table[column].map(_parse_bool).to_numpy(dtype=bool)
            else:
                data[column] = table[column].to_numpy(dtype=float)
    return Landscape(data=data)


def read_landscape(
    path_or_run_dir: str | Path,
    *,
    space: str = "raw",
    canonical_id: str | None = None,
) -> Landscape:
    """Resolve and read a landscape from a run directory or explicit artifact."""

    path = resolve_landscape_path(
        path_or_run_dir,
        space=space,
        canonical_id=canonical_id,
    )
    suffix = path.suffix.lower()
    if suffix == ".npz":
        return read_landscape_npz(path)
    if suffix == ".csv":
        return read_landscape_csv(path)
    if suffix == ".json":
        return read_landscape_json(path)
    raise ValueError(f"Unsupported landscape artifact suffix: {suffix}")


def read_landscape_metadata(
    path_or_run_dir: str | Path,
    *,
    space: str = "raw",
    canonical_id: str | None = None,
) -> dict[str, Any]:
    """Read lightweight metadata for any supported landscape artifact."""

    path = resolve_landscape_path(
        path_or_run_dir,
        space=space,
        canonical_id=canonical_id,
    )
    if path.suffix.lower() == ".json":
        metadata = read_landscape_json_metadata(path)
        metadata["path"] = str(path)
        return metadata
    if path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as payload:
            return {
                "artifact_type": str(payload["artifact_type"].item()),
                "schema_version": str(payload["schema_version"].item()),
                "row_count": int(len(payload["particle_key"])),
                "path": str(path),
            }
    if path.suffix.lower() == ".csv":
        row_count = max(sum(1 for _ in path.open("r", encoding="utf-8")) - 1, 0)
        return {
            "artifact_type": "landscape_csv",
            "schema_version": "1",
            "row_count": row_count,
            "path": str(path),
        }
    raise ValueError(f"Unsupported landscape artifact suffix: {path.suffix}")


def resolve_landscape_path(
    path_or_run_dir: str | Path,
    *,
    space: str = "raw",
    canonical_id: str | None = None,
) -> Path:
    """Resolve run-dir landscape defaults or return an explicit artifact path."""

    path = Path(path_or_run_dir)
    if path.is_file():
        return path
    if not path.is_dir():
        raise ValueError(f"Landscape path or run directory does not exist: {path}")
    if space not in {"raw", "canonical"}:
        raise ValueError("space must be 'raw' or 'canonical'")
    if space == "raw":
        return _first_existing(
            (
                path / "data" / "raw_landscape.npz",
                path / "data" / "raw_landscape.csv",
                path / "debug" / "landscape_debug.json",
                path / "landscape.json",
            )
        )
    canonical = canonical_id or "default"
    return _first_existing(
        (
            path / "canonical" / canonical / "canonical_landscape.npz",
            path / "canonical" / canonical / "canonical_landscape.csv",
            path / "canonical_landscape.npz",
            path / "canonical_landscape.csv",
            path / "debug" / f"canonical_{canonical}_landscape_debug.json",
            path / "canonical_landscape.json",
        )
    )


def _first_existing(candidates: tuple[Path, ...]) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    formatted = ", ".join(str(candidate) for candidate in candidates)
    raise ValueError(f"No supported landscape artifact found; checked: {formatted}")


def _prepare_output_path(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)


def _stack_coordinates(column, column_name: str) -> np.ndarray:
    coordinates = np.vstack([np.asarray(value, dtype=float) for value in column])
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(f"{column_name} must contain length-3 coordinates")
    return coordinates


def _integer_column_or_default(data: pd.DataFrame, column: str) -> np.ndarray:
    if column not in data.columns:
        return np.full(len(data), -1, dtype=np.int64)
    return np.asarray(data[column], dtype=np.int64)


def _display_outlier_from_npz_payload(payload) -> np.ndarray:
    if "sld_display_is_outlier" in payload:
        return np.asarray(payload["sld_display_is_outlier"], dtype=bool)
    if "sld_display_was_clipped" in payload:
        return np.asarray(payload["sld_display_was_clipped"], dtype=bool)
    return np.zeros(len(payload["particle_key"]), dtype=bool)


def _display_outlier_from_csv_table(table: pd.DataFrame) -> np.ndarray:
    if "sld_display_is_outlier" in table.columns:
        return table["sld_display_is_outlier"].map(_parse_bool).to_numpy(dtype=bool)
    if "sld_display_was_clipped" in table.columns:
        return table["sld_display_was_clipped"].map(_parse_bool).to_numpy(dtype=bool)
    return np.zeros(len(table), dtype=bool)


def _coordinates_from_columns(
    table: pd.DataFrame,
    *,
    preferred: tuple[str, str, str],
    aliases: tuple[str, str, str],
    canonical_aliases: tuple[str, str, str],
) -> np.ndarray:
    for columns in (preferred, aliases, canonical_aliases):
        if set(columns).issubset(table.columns):
            return table[list(columns)].to_numpy(dtype=float)
    raise ValueError(
        "Landscape CSV missing raw coordinate columns; expected raw_rv_*_rad "
        "or compatible rotvec aliases"
    )


def _parse_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().casefold()
        if normalized in {"true", "1", "yes"}:
            return True
        if normalized in {"false", "0", "no"}:
            return False
    return bool(value)


def _canonical_csv_header(euler_sequence: str) -> list[str]:
    return [
        "particle_key",
        "ref_source_row_id",
        "mov_source_row_id",
        "raw_rv_x_rad",
        "raw_rv_y_rad",
        "raw_rv_z_rad",
        "raw_ea_zyx_alpha_deg",
        "raw_ea_zyx_beta_deg",
        "raw_ea_zyx_gamma_deg",
        "canonical_rv_x_rad",
        "canonical_rv_y_rad",
        "canonical_rv_z_rad",
        "canonical_ea_zyx_alpha_deg",
        "canonical_ea_zyx_beta_deg",
        "canonical_ea_zyx_gamma_deg",
        "raw_rotvec_x",
        "raw_rotvec_y",
        "raw_rotvec_z",
        "canonical_rotvec_x",
        "canonical_rotvec_y",
        "canonical_rotvec_z",
        "canonical_euler_alpha",
        "canonical_euler_beta",
        "canonical_euler_gamma",
        "sld_unfloored",
        "sld_raw",
        "sld_display",
        "sld_display_is_outlier",
        "sld_was_floored",
        "sld_local_k_mean",
        "sld_effective_local_k_mean",
        "sld_distance_floor",
        "coordinate_source",
    ]


def _write_canonical_csv_chunk(
    writer: csv.writer,
    arrays: LandscapeArrays,
    *,
    start: int,
    stop: int,
    euler_sequence: str | None,
    euler_degrees: bool,
) -> None:
    raw = arrays.coordinates_analysis[start:stop]
    canonical = arrays.coordinates_canonical[start:stop]
    resolved_euler = resolve_euler_convention(scipy_euler_sequence=euler_sequence)
    raw_euler = Rotation.from_rotvec(raw).as_euler(
        resolved_euler.scipy_euler_sequence,
        degrees=euler_degrees,
    )
    canonical_euler = Rotation.from_rotvec(canonical).as_euler(
        resolved_euler.scipy_euler_sequence,
        degrees=euler_degrees,
    )
    ref_source_row_id = _source_rows_for_chunk(arrays.ref_source_row_id, start, stop)
    mov_source_row_id = _source_rows_for_chunk(arrays.mov_source_row_id, start, stop)
    for offset, row_index in enumerate(range(start, stop)):
        writer.writerow(
            [
                arrays.particle_key[row_index],
                ref_source_row_id[offset],
                mov_source_row_id[offset],
                raw[offset, 0],
                raw[offset, 1],
                raw[offset, 2],
                raw_euler[offset, 0],
                raw_euler[offset, 1],
                raw_euler[offset, 2],
                canonical[offset, 0],
                canonical[offset, 1],
                canonical[offset, 2],
                canonical_euler[offset, 0],
                canonical_euler[offset, 1],
                canonical_euler[offset, 2],
                raw[offset, 0],
                raw[offset, 1],
                raw[offset, 2],
                canonical[offset, 0],
                canonical[offset, 1],
                canonical[offset, 2],
                canonical_euler[offset, 0],
                canonical_euler[offset, 1],
                canonical_euler[offset, 2],
                arrays.sld_unfloored[row_index],
                arrays.sld_raw[row_index],
                arrays.sld_display[row_index],
                bool(arrays.sld_display_is_outlier[row_index]),
                bool(arrays.sld_was_floored[row_index]),
                arrays.sld_local_k_mean[row_index],
                arrays.sld_effective_local_k_mean[row_index],
                arrays.sld_distance_floor[row_index],
                "canonical",
            ]
        )


def _source_rows_for_chunk(
    values: np.ndarray | None,
    start: int,
    stop: int,
) -> np.ndarray:
    if values is None:
        return np.full(stop - start, -1, dtype=np.int64)
    return values[start:stop]
