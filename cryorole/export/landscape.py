"""JSON persistence helpers for landscape artifacts and reports."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from cryorole.core.euler_conventions import resolve_euler_convention
from cryorole.export.serialization import to_json_safe
from cryorole.models.canonicalization_report import CanonicalizationReport
from cryorole.models.density_report import DensityReport
from cryorole.models.landscape import Landscape


def write_landscape_json(
    landscape: Landscape,
    path: str | Path,
    *,
    overwrite: bool = False,
    artifact_type: str = "landscape",
) -> Path:
    """Write a portable v1 Landscape JSON artifact.

    This intentionally favors auditability and dependency-light portability over
    storage efficiency. It is not optimized for million-particle production
    persistence; future production writers should add chunked/binary formats.
    """

    output_path = Path(path)
    _prepare_output_path(output_path, overwrite=overwrite)
    payload = {
        "artifact_type": artifact_type,
        "schema_version": "1",
        "persistence_note": (
            "Portable JSON v1 format for auditability; not optimized for "
            "million-particle production persistence."
        ),
        "row_count": len(landscape.data),
        "data": landscape.data.to_dict(orient="records"),
        "canonical_transform": landscape.canonical_transform,
        "active_policies": landscape.active_policies or {},
        "density_report": landscape.density_report,
        "canonicalization_report": landscape.canonicalization_report,
    }
    _write_json(output_path, payload)
    return output_path


def write_json_artifact(
    payload: Any,
    path: str | Path,
    *,
    overwrite: bool = False,
) -> Path:
    """Write a JSON-safe artifact payload."""

    output_path = Path(path)
    _prepare_output_path(output_path, overwrite=overwrite)
    _write_json(output_path, payload)
    return output_path


def write_raw_landscape_csv(
    landscape: Landscape,
    path: str | Path,
    *,
    overwrite: bool = False,
    euler_sequence: str | None = None,
    euler_degrees: bool = True,
) -> Path:
    """Write a complete user-facing raw analysis landscape table."""

    output_path = Path(path)
    _prepare_output_path(output_path, overwrite=overwrite)
    data = landscape.data
    coordinates = _stack_coordinate_column(data["coordinates_analysis"], "coordinates_analysis")
    resolved_euler = resolve_euler_convention(scipy_euler_sequence=euler_sequence)
    euler = Rotation.from_rotvec(coordinates).as_euler(
        resolved_euler.scipy_euler_sequence,
        degrees=euler_degrees,
    )
    angles_deg = np.degrees(np.linalg.norm(coordinates, axis=1))
    output = pd.DataFrame(
        {
            "particle_key": data["particle_key"],
            "ref_source_row_id": data.get("ref_source_row_id", -1),
            "mov_source_row_id": data.get("mov_source_row_id", -1),
            "raw_rv_x_rad": coordinates[:, 0],
            "raw_rv_y_rad": coordinates[:, 1],
            "raw_rv_z_rad": coordinates[:, 2],
            "raw_angle_deg": angles_deg,
            "raw_ea_zyx_alpha_deg": euler[:, 0],
            "raw_ea_zyx_beta_deg": euler[:, 1],
            "raw_ea_zyx_gamma_deg": euler[:, 2],
            # Backward-compatible aliases retained for existing scripts/tests.
            "rotvec_x": coordinates[:, 0],
            "rotvec_y": coordinates[:, 1],
            "rotvec_z": coordinates[:, 2],
            "euler_alpha": euler[:, 0],
            "euler_beta": euler[:, 1],
            "euler_gamma": euler[:, 2],
            "sld_unfloored": data["sld_unfloored"],
            "sld_raw": data["sld_raw"],
            "sld_display": data["sld_display"],
            "sld_display_is_outlier": data["sld_display_is_outlier"],
            "sld_was_floored": data["sld_was_floored"],
            "sld_local_k_mean": data["sld_local_k_mean"],
            "sld_effective_local_k_mean": data["sld_effective_local_k_mean"],
            "sld_distance_floor": data["sld_distance_floor"],
            "coordinate_source": "analysis",
        }
    )
    output.to_csv(output_path, index=False)
    return output_path


def write_canonical_landscape_csv(
    landscape: Landscape,
    path: str | Path,
    *,
    overwrite: bool = False,
    euler_sequence: str | None = None,
    euler_degrees: bool = True,
) -> Path:
    """Write a complete user-facing canonical landscape table."""

    if "coordinates_canonical" not in landscape.data.columns:
        raise ValueError("canonical_landscape.csv requires coordinates_canonical")
    output_path = Path(path)
    _prepare_output_path(output_path, overwrite=overwrite)
    data = landscape.data
    raw_coordinates = _stack_coordinate_column(
        data["coordinates_analysis"],
        "coordinates_analysis",
    )
    canonical_coordinates = _stack_coordinate_column(
        data["coordinates_canonical"],
        "coordinates_canonical",
    )
    resolved_euler = resolve_euler_convention(scipy_euler_sequence=euler_sequence)
    raw_euler = Rotation.from_rotvec(raw_coordinates).as_euler(
        resolved_euler.scipy_euler_sequence,
        degrees=euler_degrees,
    )
    canonical_euler = Rotation.from_rotvec(canonical_coordinates).as_euler(
        resolved_euler.scipy_euler_sequence,
        degrees=euler_degrees,
    )
    output = pd.DataFrame(
        {
            "particle_key": data["particle_key"],
            "ref_source_row_id": data.get("ref_source_row_id", -1),
            "mov_source_row_id": data.get("mov_source_row_id", -1),
            "raw_rv_x_rad": raw_coordinates[:, 0],
            "raw_rv_y_rad": raw_coordinates[:, 1],
            "raw_rv_z_rad": raw_coordinates[:, 2],
            "raw_ea_zyx_alpha_deg": raw_euler[:, 0],
            "raw_ea_zyx_beta_deg": raw_euler[:, 1],
            "raw_ea_zyx_gamma_deg": raw_euler[:, 2],
            "canonical_rv_x_rad": canonical_coordinates[:, 0],
            "canonical_rv_y_rad": canonical_coordinates[:, 1],
            "canonical_rv_z_rad": canonical_coordinates[:, 2],
            "canonical_ea_zyx_alpha_deg": canonical_euler[:, 0],
            "canonical_ea_zyx_beta_deg": canonical_euler[:, 1],
            "canonical_ea_zyx_gamma_deg": canonical_euler[:, 2],
            # Backward-compatible aliases retained for existing scripts/tests.
            "raw_rotvec_x": raw_coordinates[:, 0],
            "raw_rotvec_y": raw_coordinates[:, 1],
            "raw_rotvec_z": raw_coordinates[:, 2],
            "canonical_rotvec_x": canonical_coordinates[:, 0],
            "canonical_rotvec_y": canonical_coordinates[:, 1],
            "canonical_rotvec_z": canonical_coordinates[:, 2],
            "canonical_euler_alpha": canonical_euler[:, 0],
            "canonical_euler_beta": canonical_euler[:, 1],
            "canonical_euler_gamma": canonical_euler[:, 2],
            "sld_unfloored": data["sld_unfloored"],
            "sld_raw": data["sld_raw"],
            "sld_display": data["sld_display"],
            "sld_display_is_outlier": data["sld_display_is_outlier"],
            "sld_was_floored": data["sld_was_floored"],
            "sld_local_k_mean": data["sld_local_k_mean"],
            "sld_effective_local_k_mean": data["sld_effective_local_k_mean"],
            "sld_distance_floor": data["sld_distance_floor"],
            "coordinate_source": "canonical",
        }
    )
    output.to_csv(output_path, index=False)
    return output_path


def read_landscape_json(path: str | Path) -> Landscape:
    """Read a portable v1 Landscape JSON artifact."""

    payload = _read_landscape_payload(path)
    data = pd.DataFrame(payload.get("data", []))
    for column in (
        "coordinates_analysis",
        "coordinates_display",
        "coordinates_canonical",
    ):
        if column in data.columns:
            data[column] = data[column].map(lambda value: np.asarray(value, dtype=float))
    return Landscape(
        data=data,
        canonical_transform=_optional_array(payload.get("canonical_transform")),
        active_policies=payload.get("active_policies") or {},
        density_report=_density_report_from_payload(payload.get("density_report")),
        canonicalization_report=_canonicalization_report_from_payload(
            payload.get("canonicalization_report")
        ),
    )


def read_landscape_json_metadata(path: str | Path) -> dict[str, Any]:
    """Read validated top-level metadata for a portable v1 Landscape artifact."""

    payload = _read_landscape_payload(path)
    return {
        "artifact_type": payload.get("artifact_type"),
        "schema_version": payload.get("schema_version"),
        "row_count": payload.get("row_count"),
    }


def write_report_json(
    report: Any,
    path: str | Path,
    *,
    overwrite: bool = False,
    artifact_type: str,
) -> Path:
    """Write a JSON-safe report artifact."""

    output_path = Path(path)
    _prepare_output_path(output_path, overwrite=overwrite)
    payload = {
        "artifact_type": artifact_type,
        "schema_version": "1",
        "report": report,
    }
    _write_json(output_path, payload)
    return output_path


def _prepare_output_path(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)


def _write_json(path: Path, payload: Any) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(to_json_safe(payload), handle, indent=2, sort_keys=True)
        handle.write("\n")


def _read_landscape_payload(path: str | Path) -> dict[str, Any]:
    input_path = Path(path)
    with input_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("Landscape artifact JSON must be an object")
    artifact_type = payload.get("artifact_type")
    if artifact_type not in {"landscape", "canonical_landscape"}:
        raise ValueError(
            "Landscape artifact must have artifact_type 'landscape' or "
            f"'canonical_landscape', got {artifact_type!r}"
        )
    schema_version = payload.get("schema_version")
    if schema_version != "1":
        raise ValueError(
            f"Landscape artifact schema_version must be '1', got {schema_version!r}"
        )
    data = payload.get("data", [])
    row_count = payload.get("row_count")
    if row_count is not None and row_count != len(data):
        raise ValueError(
            "Landscape artifact row_count does not match data length: "
            f"row_count={row_count}, data_length={len(data)}"
        )
    return payload


def _optional_array(value: Any) -> np.ndarray | None:
    if value is None:
        return None
    return np.asarray(value, dtype=float)


def _stack_coordinate_column(column, column_name: str) -> np.ndarray:
    coordinates = np.vstack([np.asarray(value, dtype=float) for value in column])
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(f"{column_name} must contain length-3 coordinates")
    return coordinates


def _density_report_from_payload(value: Any) -> DensityReport | None:
    if value is None:
        return None
    data = dict(value)
    _normalize_density_report_payload(data)
    data["floored_particle_keys"] = tuple(data.get("floored_particle_keys", ()))
    data["floored_rows"] = tuple(data.get("floored_rows", ()))
    data["warnings"] = tuple(data.get("warnings", ()))
    return DensityReport(**data)


def _normalize_density_report_payload(data: dict[str, Any]) -> None:
    """Tolerate transitional percentile-clipping density reports."""

    clip_vmax = data.pop("sld_display_clip_vmax", None)
    clipped_count = data.pop("n_sld_display_clipped", None)
    clipped_fraction = data.pop("fraction_sld_display_clipped", None)
    max_over_clip = data.pop("max_over_display_clip_vmax", None)
    data.pop("sld_display_clip_percentile", None)
    if data.get("sld_display_mode") == "robust_percentile_clip":
        data["sld_display_mode"] = "identity"
    if "sld_display_color_vmax" not in data and clip_vmax is not None:
        data["sld_display_color_vmax"] = clip_vmax
    if "n_sld_display_outliers" not in data and clipped_count is not None:
        data["n_sld_display_outliers"] = clipped_count
    if "fraction_sld_display_outliers" not in data and clipped_fraction is not None:
        data["fraction_sld_display_outliers"] = clipped_fraction
    if "max_over_display_vmax" not in data and max_over_clip is not None:
        data["max_over_display_vmax"] = max_over_clip


def _canonicalization_report_from_payload(
    value: Any,
) -> CanonicalizationReport | None:
    if value is None:
        return None
    data = dict(value)
    data["singular_values"] = tuple(data.get("singular_values", ()))
    data["explained_variance_ratios"] = tuple(
        data.get("explained_variance_ratios", ())
    )
    if "pca_axis_order" in data:
        data["pca_axis_order"] = tuple(data["pca_axis_order"])
    if "assigned_coordinate_names" in data:
        data["assigned_coordinate_names"] = tuple(data["assigned_coordinate_names"])
    if "assigned_rotvec_columns" in data:
        data["assigned_rotvec_columns"] = tuple(data["assigned_rotvec_columns"])
    if "axis_weighted_skewness" in data:
        data["axis_weighted_skewness"] = tuple(data["axis_weighted_skewness"])
    if "flipped_axes" in data:
        data["flipped_axes"] = tuple(data["flipped_axes"])
    if "ambiguous_sign_axes" in data:
        data["ambiguous_sign_axes"] = tuple(data["ambiguous_sign_axes"])
    if "handedness_adjusted_axes" in data:
        data["handedness_adjusted_axes"] = tuple(data["handedness_adjusted_axes"])
    data["warnings"] = tuple(data.get("warnings", ()))
    return CanonicalizationReport(**data)
