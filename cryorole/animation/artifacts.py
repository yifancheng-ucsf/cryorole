"""Strict, compact artifact resolution for offline animation."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from cryorole.canonicalize.transforms import validate_canonical_transform
from cryorole.core.display_policy import (
    ColorScale,
    EULER_AXIS_NAMES,
    resolve_color_scale,
    resolve_display_indices,
)
from cryorole.core.euler_conventions import (
    CANONICAL_EULER_ANGLE_COLUMNS,
    RAW_EULER_ANGLE_COLUMNS,
    EulerConventionResolution,
    resolve_euler_convention,
)
CANONICAL_FRAME_TRANSFORM_DIRECTION = (
    "canonical_rv = raw_rv @ canonical_transform"
)
CANONICAL_FRAME_COORDINATE_SPACE = "rotvec_ro_radians"


@dataclass(frozen=True)
class AnimationLandscape:
    """Compact plotting fields loaded once for a selected landscape."""

    source_path: Path
    euler_degrees: np.ndarray
    sld_display: np.ndarray
    sld_display_is_outlier: np.ndarray
    total_rows: int
    displayed_rows: int
    source_row_indices: np.ndarray
    color_scale: ColorScale
    sld_field: str = "sld_display"


@dataclass(frozen=True)
class CanonicalFrame:
    """Validated canonical frame plus explicit vector-convention metadata."""

    source_path: Path
    transform: np.ndarray
    canonical_id: str
    physical_change_of_basis: bool
    transform_direction: str = CANONICAL_FRAME_TRANSFORM_DIRECTION
    coordinate_space: str = CANONICAL_FRAME_COORDINATE_SPACE
    row_vector_convention: str = "canonical_rv = raw_rv @ C"
    raw_to_canonical_column_convention: str = "x_canonical = C.T @ x_raw"


def resolve_recorded_euler_convention(
    run_dir: str | Path,
    *,
    coordinate_set: str,
    canonical_id: str | None = None,
    requested: str = "auto",
) -> EulerConventionResolution:
    """Resolve Euler policy strictly from selected-landscape provenance."""

    root = Path(run_dir)
    payloads: list[tuple[Path, dict[str, object]]] = []
    if coordinate_set == "canonical":
        canonical = canonical_id or "default"
        for name in ("canonicalize_summary.json", "canonicalization_report.json"):
            path = root / "canonical" / canonical / name
            payload = _read_json_object(path)
            if payload is not None:
                payloads.append((path, payload))
    elif coordinate_set != "raw":
        raise ValueError("coordinate_set must be 'raw' or 'canonical'")
    for name in ("run_summary.json", "run_manifest.json"):
        path = root / name
        payload = _read_json_object(path)
        if payload is not None:
            payloads.append((path, payload))

    found: list[tuple[Path, EulerConventionResolution]] = []
    for path, payload in payloads:
        metadata_candidates = [payload]
        results = payload.get("results")
        if isinstance(results, dict):
            euler_metadata = results.get("euler_metadata")
            if isinstance(euler_metadata, dict):
                metadata_candidates.append(euler_metadata)
        for metadata in metadata_candidates:
            convention = metadata.get("euler_convention")
            sequence = metadata.get("scipy_euler_sequence")
            if isinstance(convention, str) or isinstance(sequence, str):
                found.append(
                    (
                        path,
                        resolve_euler_convention(
                            convention if isinstance(convention, str) else None,
                            scipy_euler_sequence=sequence if isinstance(sequence, str) else None,
                            source=f"inherited:{path}",
                        ),
                    )
                )
                break
    if not found:
        raise ValueError(
            "Selected landscape has no recorded Euler convention provenance; "
            "animation refuses to infer intrinsic/extrinsic semantics"
        )
    resolved = found[0][1]
    for path, candidate in found[1:]:
        if candidate.euler_convention != resolved.euler_convention:
            raise ValueError(
                "Conflicting Euler convention provenance: "
                f"{resolved.euler_convention!r} vs {candidate.euler_convention!r} in {path}"
            )
    if requested != "auto" and requested != resolved.euler_convention:
        raise ValueError(
            f"Explicit Euler convention {requested!r} conflicts with selected "
            f"landscape provenance {resolved.euler_convention!r}"
        )
    return resolved


def load_canonical_frame(run_dir: str | Path, canonical_id: str) -> CanonicalFrame:
    """Load and audit the canonical frame belonging to ``canonical_id``."""

    root = Path(run_dir)
    path = root / "canonical" / canonical_id / "canonical_frame.json"
    payload = _read_json_object(path)
    if payload is None:
        raise ValueError(f"Canonical frame does not exist or is not a JSON object: {path}")
    if payload.get("transform_direction") != CANONICAL_FRAME_TRANSFORM_DIRECTION:
        raise ValueError(
            "canonical_frame transform_direction must be "
            f"{CANONICAL_FRAME_TRANSFORM_DIRECTION!r}"
        )
    if payload.get("coordinate_space") != CANONICAL_FRAME_COORDINATE_SPACE:
        raise ValueError(
            f"canonical_frame coordinate_space must be {CANONICAL_FRAME_COORDINATE_SPACE!r}"
        )
    source_id = payload.get("source_canonical_id")
    if source_id not in {None, canonical_id}:
        raise ValueError(
            f"canonical_frame source_canonical_id {source_id!r} does not match {canonical_id!r}"
        )
    if "canonical_transform" in payload:
        transform = np.asarray(payload["canonical_transform"], dtype=float)
    else:
        npz_value = payload.get("canonical_frame_npz")
        if not isinstance(npz_value, str):
            raise ValueError("canonical_frame requires canonical_transform or canonical_frame_npz")
        npz_path = Path(npz_value)
        if not npz_path.is_absolute():
            npz_path = path.parent / npz_path
        with np.load(npz_path, allow_pickle=False) as archive:
            if "canonical_transform" not in archive:
                raise ValueError("canonical_frame.npz missing canonical_transform")
            transform = np.asarray(archive["canonical_transform"], dtype=float)
    validate_canonical_transform(transform)
    return CanonicalFrame(
        source_path=path,
        transform=transform,
        canonical_id=canonical_id,
        physical_change_of_basis=payload.get("physical_change_of_basis") is True,
    )


def load_animation_landscape(
    run_dir: str | Path,
    *,
    coordinate_set: str,
    canonical_id: str | None,
    scipy_euler_sequence: str,
    sld_threshold: float | None,
    top_fraction: float | None,
    range_bounds: dict[str, tuple[float | None, float | None]] | None = None,
    color_vmin: float | None = None,
    color_vmax: float | None = None,
) -> AnimationLandscape:
    """Load only projection coordinates and fixed ``sld_display`` values."""

    from cryorole.io.landscape_resolver import resolve_landscape_path

    source = resolve_landscape_path(
        run_dir,
        space=coordinate_set,
        canonical_id=canonical_id,
    )
    euler, sld, outlier = _read_projection_fields(
        source,
        coordinate_set=coordinate_set,
        scipy_euler_sequence=scipy_euler_sequence,
    )
    total = len(sld)
    color_scale = resolve_color_scale(
        sld,
        display_outlier_mask=outlier,
        visual_style="legacy",
        color_vmin=color_vmin,
        color_vmax=color_vmax,
        display_threshold=sld_threshold,
    )
    unknown_ranges = sorted(set(range_bounds or {}) - set(EULER_AXIS_NAMES))
    if unknown_ranges:
        raise ValueError(f"Animation Euler ranges use unsupported axes: {unknown_ranges}")
    selected = resolve_display_indices(
        sld,
        threshold=sld_threshold,
        top_fraction=top_fraction,
        coordinates=euler,
        axis_names=EULER_AXIS_NAMES,
        range_bounds=range_bounds,
        wraparound=True,
    )
    if not selected.size:
        raise ValueError("Display filter removed every landscape row")
    return AnimationLandscape(
        source_path=source,
        euler_degrees=np.asarray(euler[selected], dtype=np.float32),
        sld_display=np.asarray(sld[selected], dtype=np.float32),
        sld_display_is_outlier=np.asarray(outlier[selected], dtype=bool),
        total_rows=total,
        displayed_rows=int(selected.size),
        source_row_indices=selected,
        color_scale=color_scale,
    )


def _read_projection_fields(
    path: Path,
    *,
    coordinate_set: str,
    scipy_euler_sequence: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if path.suffix.lower() == ".npz":
        with np.load(path, allow_pickle=False) as payload:
            coordinate_key = (
                "coordinates_canonical" if coordinate_set == "canonical" else "coordinates_analysis"
            )
            if coordinate_key not in payload:
                raise ValueError(f"Landscape NPZ missing {coordinate_key!r}")
            rv = np.asarray(payload[coordinate_key], dtype=np.float64)
            sld = np.asarray(payload["sld_display"], dtype=np.float64)
            outlier = (
                np.asarray(payload["sld_display_is_outlier"], dtype=bool)
                if "sld_display_is_outlier" in payload
                else np.zeros(len(sld), dtype=bool)
            )
        euler = Rotation.from_rotvec(rv).as_euler(scipy_euler_sequence, degrees=True)
        return euler, sld, outlier
    if path.suffix.lower() != ".csv":
        raise ValueError("Animation supports production NPZ/CSV landscapes only")
    euler_columns = (
        CANONICAL_EULER_ANGLE_COLUMNS
        if coordinate_set == "canonical"
        else RAW_EULER_ANGLE_COLUMNS
    )
    header = pd.read_csv(path, nrows=0).columns
    required = list(euler_columns) + ["sld_display"]
    if not set(required).issubset(header):
        raise ValueError(
            "Landscape CSV lacks recorded-space Euler projection columns: "
            + ", ".join(euler_columns)
        )
    usecols = required + (
        ["sld_display_is_outlier"] if "sld_display_is_outlier" in header else []
    )
    table = pd.read_csv(path, usecols=usecols, dtype={column: "float32" for column in required})
    euler = table[list(euler_columns)].to_numpy(dtype=np.float64)
    sld = table["sld_display"].to_numpy(dtype=np.float64)
    outlier = (
        table["sld_display_is_outlier"].astype(bool).to_numpy()
        if "sld_display_is_outlier" in table
        else np.zeros(len(table), dtype=bool)
    )
    return euler, sld, outlier


def _read_json_object(path: Path) -> dict[str, object] | None:
    if not path.is_file():
        return None
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, dict) else None
