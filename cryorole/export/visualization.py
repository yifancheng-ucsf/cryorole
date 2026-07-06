"""Display-only landscape visualization artifact helpers."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from cryorole.core.euler_conventions import RAW_EULER_ANGLE_COLUMNS, resolve_euler_convention
from cryorole.export.landscape import write_json_artifact
from cryorole.models.landscape import Landscape


ROT_VECTOR_PROJECTIONS = (
    ("xy", 0, 1),
    ("yz", 1, 2),
    ("xz", 0, 2),
)
ROT_VECTOR_AXIS_NAMES = ("x", "y", "z")
EULER_PROJECTIONS = (
    ("alpha_beta", 0, 1),
    ("beta_gamma", 1, 2),
    ("alpha_gamma", 0, 2),
)
EULER_AXIS_NAMES = ("alpha", "beta", "gamma")
SLD_FIELDS = (
    "sld_unfloored",
    "sld_raw",
    "sld_display",
    "sld_display_is_outlier",
    "sld_was_floored",
    "sld_local_k_mean",
    "sld_effective_local_k_mean",
    "sld_distance_floor",
)
STYLE_PRESETS: dict[str, dict[str, Any]] = {
    "modern": {
        "color_map": "viridis",
        "point_size": 8.0,
        "point_alpha": 0.85,
        "figure_size_2d": (12.0, 4.0),
        "figure_size_3d": (6.0, 5.0),
        "colorbar_position": "right",
        "aspect": "auto",
        "sort_points_by_color": "none",
        "axis_limits": {
            "alpha": (-180.0, 180.0),
            "beta": (-180.0, 180.0),
            "gamma": (-180.0, 180.0),
        },
    },
    "legacy": {
        "color_map": "rainbow_r",
        "point_size": 1.0,
        "point_alpha": 1.0,
        "figure_size_2d": (20.0, 9.0),
        "figure_size_3d": (8.0, 7.0),
        "colorbar_position": "bottom",
        "aspect": "equal",
        "sort_points_by_color": "ascending",
        "axis_limits": {
            "alpha": (-180.0, 180.0),
            "beta": (-180.0, 180.0),
            "gamma": (-180.0, 180.0),
            "x": (-3.14, 3.14),
            "y": (-3.14, 3.14),
            "z": (-3.14, 3.14),
        },
    },
    "paper": {
        "color_map": "viridis",
        "point_size": 3.0,
        "point_alpha": 0.9,
        "figure_size_2d": (14.0, 5.0),
        "figure_size_3d": (6.0, 5.0),
        "colorbar_position": "right",
        "aspect": "equal",
        "sort_points_by_color": "ascending",
        "axis_limits": {
            "alpha": (-180.0, 180.0),
            "beta": (-180.0, 180.0),
            "gamma": (-180.0, 180.0),
        },
    },
}


def write_landscape_visualizations(
    landscape: Landscape,
    output_dir: str | Path,
    *,
    overwrite: bool = False,
    coordinate_source: str = "all_available",
    representation: str = "both",
    color_field: str = "sld_display",
    euler_sequence: str | None = None,
    euler_convention: str | None = None,
    euler_convention_source: str = "cli_default",
    euler_degrees: bool = True,
    display_top_fraction: float | None = None,
    display_sld_threshold: float | None = None,
    display_density_field: str = "sld_display",
    range_bounds: Mapping[str, tuple[float | None, float | None] | None] | None = None,
    formats: Sequence[str] = ("png", "svg", "pdf"),
    max_points_2d: int | None = None,
    max_points_3d: int = 50000,
    random_seed: int = 0,
    write_projection_csvs: bool = False,
    output_prefix: str = "",
    full_landscape_table: str | None = None,
    euler_axis_limits: tuple[float, float] | None = (-180.0, 180.0),
    visual_style: str = "modern",
    color_map: str | None = None,
    color_vmin: float | None = None,
    color_vmax: float | None = None,
    point_size: float | None = None,
    point_alpha: float | None = None,
    figure_width: float | None = None,
    figure_height: float | None = None,
    colorbar_position: str | None = None,
    sort_points_by_color: str | None = None,
    axis_limits: Mapping[str, tuple[float | None, float | None] | None] | None = None,
    display_filter_mode: str | None = None,
    display_max_divisor: float = 3.0,
    generate_histograms: bool = False,
    generate_axis_direction_map: bool = False,
    display_table_filename: str | None = None,
    artifact_layout: str = "flat",
    selection_metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Write display-only table, static figures, and optional debug CSVs."""

    _ensure_matplotlib_available()
    resolved_euler = resolve_euler_convention(
        euler_convention,
        scipy_euler_sequence=euler_sequence,
        source=euler_convention_source,
    )
    _validate_visualization_inputs(
        landscape,
        representation=representation,
        color_field=color_field,
        display_top_fraction=display_top_fraction,
        display_sld_threshold=display_sld_threshold,
        display_density_field=display_density_field,
        formats=formats,
        max_points_2d=max_points_2d,
        max_points_3d=max_points_3d,
        visual_style=visual_style,
        color_vmin=color_vmin,
        color_vmax=color_vmax,
        point_size=point_size,
        point_alpha=point_alpha,
        figure_width=figure_width,
        figure_height=figure_height,
        colorbar_position=colorbar_position,
        sort_points_by_color=sort_points_by_color,
        display_filter_mode=display_filter_mode,
        display_max_divisor=display_max_divisor,
        artifact_layout=artifact_layout,
    )
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    data = landscape.data.copy(deep=True)
    coordinate_sources = _resolve_coordinate_sources(data, coordinate_source)
    representations = _resolve_representations(representation)
    normalized_formats = tuple(_normalize_formats(formats))
    range_bounds = dict(range_bounds or {})
    _validate_range_axes(range_bounds, representations)
    style_options = _resolve_style_options(
        visual_style=visual_style,
        color_map=color_map,
        point_size=point_size,
        point_alpha=point_alpha,
        figure_width=figure_width,
        figure_height=figure_height,
        colorbar_position=colorbar_position,
        sort_points_by_color=sort_points_by_color,
        axis_limits=axis_limits,
        euler_axis_limits=euler_axis_limits,
    )
    _validate_axis_limits(style_options["axis_limits"], representations)
    resolved_display_filter_mode, resolved_display_threshold = _resolve_display_filter(
        data[display_density_field],
        display_top_fraction=display_top_fraction,
        display_sld_threshold=display_sld_threshold,
        display_filter_mode=display_filter_mode,
        display_max_divisor=display_max_divisor,
    )
    (
        resolved_color_vmin,
        resolved_color_vmax,
        color_vmin_source,
        color_vmax_source,
    ) = _resolve_color_scale(
        data[color_field],
        display_outlier_mask=_display_outlier_mask(data, color_field),
        visual_style=visual_style,
        color_vmin=color_vmin,
        color_vmax=color_vmax,
        display_threshold=resolved_display_threshold,
    )
    unused_axis_limits = _unused_axis_limits(style_options["axis_limits"], representations)
    expected_paths = _expected_visualization_paths(
        output_path,
        coordinate_sources,
        representations,
        normalized_formats,
        write_projection_csvs=write_projection_csvs,
        output_prefix=output_prefix,
        generate_histograms=generate_histograms,
        generate_axis_direction_map=generate_axis_direction_map,
        artifact_layout=artifact_layout,
        display_table_filename=display_table_filename,
    )
    _ensure_visualization_paths_available(expected_paths, overwrite=overwrite)

    generated_files: dict[str, str] = {}
    display_tables = []
    n_after_filter_by_source: dict[str, int] = {}
    n_2d_by_source: dict[str, int] = {}
    n_3d_by_source: dict[str, int] = {}
    downsampled_by_source: dict[str, bool] = {}
    generated_histograms: dict[str, str] = {}
    generated_axis_direction_map: dict[str, str] = {}

    for source in coordinate_sources:
        coordinates = _stack_coordinates(data[f"coordinates_{source}"], f"coordinates_{source}")
        euler_coordinates = Rotation.from_rotvec(coordinates).as_euler(
            resolved_euler.scipy_euler_sequence,
            degrees=euler_degrees,
        )
        display_indices = _display_indices(
            data,
            coordinates,
            euler_coordinates,
            representations=representations,
            display_top_fraction=(
                display_top_fraction if resolved_display_filter_mode == "top_fraction" else None
            ),
            display_sld_threshold=resolved_display_threshold,
            display_density_field=display_density_field,
            range_bounds=range_bounds,
        )
        source_display = data.iloc[display_indices].copy(deep=True)
        source_rotvec = coordinates[display_indices]
        source_euler = euler_coordinates[display_indices]
        display_tables.append(
            _display_table(
                source_display,
                source_rotvec,
                source_euler,
                coordinate_source=source,
                representations=representations,
            )
        )
        n_after_filter_by_source[source] = int(len(display_indices))
        n_2d_by_source[source] = int(_downsample_indices(len(display_indices), max_points_2d, random_seed).size)
        indices_3d = _downsample_indices(len(display_indices), max_points_3d, random_seed)
        n_3d_by_source[source] = int(indices_3d.size)
        downsampled_by_source[source] = bool(indices_3d.size < len(display_indices))

        for rep in representations:
            rep_coordinates, projections, axis_names = _coordinates_for_representation(
                rep,
                source_rotvec,
                source_euler,
            )
            generated_files.update(
                _write_2d_figure(
                    source_display,
                    rep_coordinates,
                    output_path,
                    coordinate_source=source,
                    representation=rep,
                    projections=projections,
                    axis_names=axis_names,
                    color_field=color_field,
                    formats=normalized_formats,
                    max_points_2d=max_points_2d,
                    random_seed=random_seed,
                    output_prefix=output_prefix,
                    style_options=style_options,
                    color_vmin=resolved_color_vmin,
                    color_vmax=resolved_color_vmax,
                    artifact_layout=artifact_layout,
                )
            )
            generated_files.update(
                _write_3d_figure(
                    source_display,
                    rep_coordinates,
                    output_path,
                    coordinate_source=source,
                    representation=rep,
                    axis_names=axis_names,
                    color_field=color_field,
                    formats=normalized_formats,
                    max_points_3d=max_points_3d,
                    random_seed=random_seed,
                    output_prefix=output_prefix,
                    style_options=style_options,
                    color_vmin=resolved_color_vmin,
                    color_vmax=resolved_color_vmax,
                    artifact_layout=artifact_layout,
                )
            )
            if write_projection_csvs:
                generated_files.update(
                    _write_projection_csvs(
                        source_display,
                        rep_coordinates,
                        output_path,
                        coordinate_source=source,
                        representation=rep,
                        projections=projections,
                        axis_names=axis_names,
                        color_field=color_field,
                        overwrite=overwrite,
                        output_prefix=output_prefix,
                    )
                )
        if generate_histograms:
            histogram_files = _write_histograms(
                source_display,
                source_rotvec,
                source_euler,
                output_path,
                coordinate_source=source,
                color_field=color_field,
                formats=normalized_formats,
                output_prefix=output_prefix,
                style_options=style_options,
            )
            generated_histograms.update(histogram_files)
            generated_files.update(histogram_files)
        if generate_axis_direction_map:
            axis_map_files = _write_axis_direction_map(
                source_display,
                source_rotvec,
                output_path,
                coordinate_source=source,
                color_field=color_field,
                formats=normalized_formats,
                output_prefix=output_prefix,
                style_options=style_options,
                color_vmin=resolved_color_vmin,
                color_vmax=resolved_color_vmax,
            )
            generated_axis_direction_map.update(axis_map_files)
            generated_files.update(axis_map_files)

    display_coordinates = (
        pd.concat(display_tables, ignore_index=True)
        if display_tables
        else pd.DataFrame()
    )
    display_coordinates_name = display_table_filename or f"{output_prefix}display_coordinates.csv"
    display_coordinates_path = output_path / display_coordinates_name
    display_coordinates.to_csv(display_coordinates_path, index=False)
    generated_files["display_coordinates_csv"] = str(display_coordinates_path)

    report = {
        "artifact_type": "visualization_report",
        "schema_version": "1",
        "status": "ok",
        "backend": "matplotlib_static",
        "plot_format": "static_figure",
        "matplotlib_required": True,
        "coordinate_source_requested": coordinate_source,
        "coordinate_source_resolved": coordinate_sources,
        "coordinate_sources": coordinate_sources,
        "representation": representation,
        "representations": representations,
        "euler_sequence": resolved_euler.scipy_euler_sequence,
        "scipy_euler_sequence": resolved_euler.scipy_euler_sequence,
        "euler_convention": resolved_euler.euler_convention,
        "euler_angle_columns": list(RAW_EULER_ANGLE_COLUMNS),
        "euler_convention_source": resolved_euler.euler_convention_source,
        "euler_degrees": euler_degrees,
        "visual_style": visual_style,
        "color_map": style_options["color_map"],
        "color_field": color_field,
        "color_field_usage": "display_only",
        "color_scale_mode": _color_scale_mode(data, color_field),
        "tail_search_fraction": _density_report_field(
            landscape,
            "sld_tail_search_fraction",
        ),
        "tail_jump_factor": _density_report_field(
            landscape,
            "sld_tail_jump_factor",
        ),
        "max_display_outlier_fraction": _density_report_field(
            landscape,
            "sld_max_display_outlier_fraction",
        ),
        "display_color_vmax": resolved_color_vmax,
        "n_sld_display_outliers": int(
            np.asarray(data["sld_display_is_outlier"], dtype=bool).sum()
        ),
        "fraction_sld_display_outliers": (
            float(np.asarray(data["sld_display_is_outlier"], dtype=bool).sum() / len(data))
            if len(data)
            else 0.0
        ),
        "color_vmin": resolved_color_vmin,
        "color_vmax": resolved_color_vmax,
        "color_vmin_source": color_vmin_source,
        "color_vmax_source": color_vmax_source,
        "point_size": style_options["point_size"],
        "point_alpha": style_options["point_alpha"],
        "figure_size": {
            "projection_2d": list(style_options["figure_size_2d"]),
            "cloud_3d": list(style_options["figure_size_3d"]),
        },
        "colorbar_position": style_options["colorbar_position"],
        "aspect": style_options["aspect"],
        "sort_points_by_color": style_options["sort_points_by_color"],
        "axis_limits": _json_axis_limits(style_options["axis_limits"]),
        "unused_axis_limits": _json_axis_limits(unused_axis_limits),
        "display_filter_mode": resolved_display_filter_mode,
        "display_top_fraction": (
            display_top_fraction if resolved_display_filter_mode == "top_fraction" else None
        ),
        "display_sld_threshold": resolved_display_threshold,
        "display_density_field": display_density_field,
        "display_max_divisor": (
            display_max_divisor
            if resolved_display_filter_mode == "legacy_max_divisor"
            else None
        ),
        "display_filter_threshold": resolved_display_threshold,
        "display_range_filter_applied": bool(range_bounds),
        "range_bounds": range_bounds,
        "euler_axis_limits": (
            list(style_options["axis_limits"]["alpha"])
            if "alpha" in style_options["axis_limits"]
            else None
        ),
        "n_points_input": len(landscape.data),
        "n_points_after_display_filter": n_after_filter_by_source,
        "n_points_2d": n_2d_by_source,
        "n_points_3d": n_3d_by_source,
        "downsampled": downsampled_by_source,
        "sampling_method": "deterministic_random_without_replacement",
        "random_seed": random_seed,
        "max_points_2d": max_points_2d,
        "max_points_3d": max_points_3d,
        "write_projection_csvs": write_projection_csvs,
        "display_coordinates_contains": "display_filtered_rows_only",
        "display_table_filename": display_coordinates_name,
        "artifact_layout": artifact_layout,
        "full_landscape_table": full_landscape_table,
        "figure_point_rasterized": True,
        "generated_histograms": generated_histograms,
        "generated_axis_direction_map": generated_axis_direction_map,
        "generated_files": generated_files,
    }
    if selection_metadata is not None:
        report.update(dict(selection_metadata))
    report_path = write_json_artifact(
        report,
        output_path / f"{output_prefix}visualization_report.json",
        overwrite=overwrite,
    )
    report["report_path"] = str(report_path)
    return report


def _ensure_matplotlib_available() -> None:
    try:
        import matplotlib  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is required for cryoROLE visualization outputs. "
            "Install the plotting dependency to generate static figures."
        ) from exc


def _validate_visualization_inputs(
    landscape: Landscape,
    *,
    representation: str,
    color_field: str,
    display_top_fraction: float | None,
    display_sld_threshold: float | None,
    display_density_field: str,
    formats: Sequence[str],
    max_points_2d: int | None,
    max_points_3d: int,
    visual_style: str,
    color_vmin: float | None,
    color_vmax: float | None,
    point_size: float | None,
    point_alpha: float | None,
    figure_width: float | None,
    figure_height: float | None,
    colorbar_position: str | None,
    sort_points_by_color: str | None,
    display_filter_mode: str | None,
    display_max_divisor: float,
    artifact_layout: str,
) -> None:
    if representation not in {"euler", "rotvec", "both"}:
        raise ValueError(f"Unsupported visualization representation: {representation}")
    if visual_style not in STYLE_PRESETS:
        raise ValueError(f"Unsupported visual_style: {visual_style}")
    if artifact_layout not in {"flat", "run_bundle"}:
        raise ValueError("artifact_layout must be 'flat' or 'run_bundle'")
    if color_field not in landscape.data.columns:
        raise ValueError(f"Landscape missing visualization color field: {color_field}")
    if display_density_field not in {"sld_display", "sld_raw"}:
        raise ValueError("display_density_field must be sld_display or sld_raw")
    if display_density_field not in landscape.data.columns:
        raise ValueError(f"Landscape missing display density field: {display_density_field}")
    if display_top_fraction is not None and display_sld_threshold is not None:
        raise ValueError("Use only one display density filter at a time")
    if display_filter_mode not in {None, "none", "top_fraction", "threshold", "legacy_max_divisor"}:
        raise ValueError(f"Unsupported display_filter_mode: {display_filter_mode}")
    if display_top_fraction is not None and not 0.0 < display_top_fraction <= 1.0:
        raise ValueError("display_top_fraction must be > 0 and <= 1")
    if display_sld_threshold is not None and not np.isfinite(display_sld_threshold):
        raise ValueError("display_sld_threshold must be finite")
    if display_filter_mode == "top_fraction" and display_top_fraction is None:
        raise ValueError("display_filter_mode='top_fraction' requires display_top_fraction")
    if display_filter_mode == "threshold" and display_sld_threshold is None:
        raise ValueError("display_filter_mode='threshold' requires display_sld_threshold")
    if display_filter_mode == "legacy_max_divisor" and display_max_divisor <= 0:
        raise ValueError("display_max_divisor must be positive")
    if color_vmin is not None and not np.isfinite(color_vmin):
        raise ValueError("color_vmin must be finite when provided")
    if color_vmax is not None and not np.isfinite(color_vmax):
        raise ValueError("color_vmax must be finite when provided")
    if color_vmin is not None and color_vmax is not None and color_vmin >= color_vmax:
        raise ValueError("color_vmin must be less than color_vmax")
    if point_size is not None and point_size <= 0:
        raise ValueError("point_size must be positive")
    if point_alpha is not None and not 0.0 < point_alpha <= 1.0:
        raise ValueError("point_alpha must be > 0 and <= 1")
    if figure_width is not None and figure_width <= 0:
        raise ValueError("figure_width must be positive")
    if figure_height is not None and figure_height <= 0:
        raise ValueError("figure_height must be positive")
    if colorbar_position is not None and colorbar_position not in {"bottom", "right"}:
        raise ValueError("colorbar_position must be bottom or right")
    if sort_points_by_color is not None and sort_points_by_color not in {
        "none",
        "ascending",
        "descending",
    }:
        raise ValueError("sort_points_by_color must be none, ascending, or descending")
    if not formats:
        raise ValueError("At least one figure format is required")
    if max_points_2d is not None and max_points_2d < 1:
        raise ValueError("max_points_2d must be positive when provided")
    if max_points_3d < 1:
        raise ValueError("max_points_3d must be positive")


def _resolve_coordinate_sources(data, coordinate_source: str) -> list[str]:
    if coordinate_source == "all_available":
        sources = ["analysis"]
        if "coordinates_canonical" in data.columns:
            sources.append("canonical")
        return sources
    if coordinate_source == "analysis":
        return ["analysis"]
    if coordinate_source == "canonical":
        if "coordinates_canonical" not in data.columns:
            raise ValueError(
                "coordinate_source='canonical' requires coordinates_canonical"
            )
        return ["canonical"]
    if coordinate_source == "canonical_if_available":
        if "coordinates_canonical" in data.columns:
            return ["canonical"]
        return ["analysis"]
    raise ValueError(f"Unsupported visualization coordinate_source: {coordinate_source}")


def _resolve_representations(representation: str) -> list[str]:
    if representation == "both":
        return ["euler", "rotvec"]
    return [representation]


def _normalize_formats(formats: Sequence[str]) -> list[str]:
    normalized: list[str] = []
    for fmt in formats:
        for part in str(fmt).split(","):
            clean = part.strip().lower().lstrip(".")
            if clean:
                normalized.append(clean)
    if not normalized:
        raise ValueError("At least one figure format is required")
    unsupported = sorted(set(normalized) - {"png", "svg", "pdf"})
    if unsupported:
        raise ValueError(f"Unsupported figure format(s): {unsupported}")
    return normalized


def _resolve_style_options(
    *,
    visual_style: str,
    color_map: str | None,
    point_size: float | None,
    point_alpha: float | None,
    figure_width: float | None,
    figure_height: float | None,
    colorbar_position: str | None,
    sort_points_by_color: str | None,
    axis_limits: Mapping[str, tuple[float | None, float | None] | None] | None,
    euler_axis_limits: tuple[float, float] | None,
) -> dict[str, Any]:
    preset = STYLE_PRESETS[visual_style]
    resolved_axis_limits: dict[str, tuple[float, float]] = {
        axis: tuple(bounds)
        for axis, bounds in preset["axis_limits"].items()
        if bounds is not None
    }
    if euler_axis_limits is not None and axis_limits is None:
        if len(euler_axis_limits) != 2 or euler_axis_limits[0] >= euler_axis_limits[1]:
            raise ValueError("euler_axis_limits must be a (lower, upper) pair")
        for axis in EULER_AXIS_NAMES:
            resolved_axis_limits[axis] = tuple(float(value) for value in euler_axis_limits)
    if axis_limits is not None:
        for axis, bounds in axis_limits.items():
            if bounds is None:
                resolved_axis_limits.pop(axis, None)
                continue
            lower, upper = bounds
            if lower is None or upper is None:
                raise ValueError("axis limits require both lower and upper bounds")
            if lower >= upper:
                raise ValueError("axis limit lower bound must be less than upper bound")
            resolved_axis_limits[axis] = (float(lower), float(upper))
    figure_size_2d = tuple(preset["figure_size_2d"])
    if figure_width is not None or figure_height is not None:
        figure_size_2d = (
            float(figure_width if figure_width is not None else figure_size_2d[0]),
            float(figure_height if figure_height is not None else figure_size_2d[1]),
        )
    return {
        "color_map": color_map or preset["color_map"],
        "point_size": float(point_size if point_size is not None else preset["point_size"]),
        "point_alpha": float(point_alpha if point_alpha is not None else preset["point_alpha"]),
        "figure_size_2d": figure_size_2d,
        "figure_size_3d": tuple(preset["figure_size_3d"]),
        "colorbar_position": colorbar_position or preset["colorbar_position"],
        "aspect": preset["aspect"],
        "sort_points_by_color": sort_points_by_color or preset["sort_points_by_color"],
        "axis_limits": resolved_axis_limits,
    }


def _validate_axis_limits(axis_limits: Mapping[str, tuple[float, float]], representations: list[str]) -> None:
    allowed = set(EULER_AXIS_NAMES) | set(ROT_VECTOR_AXIS_NAMES)
    unknown = sorted(axis for axis in axis_limits if axis not in allowed)
    if unknown:
        raise ValueError(f"Unsupported display axis limit for representation: {unknown}")


def _unused_axis_limits(
    axis_limits: Mapping[str, tuple[float, float]],
    representations: list[str],
) -> dict[str, tuple[float, float]]:
    active_axes = set()
    if "euler" in representations:
        active_axes.update(EULER_AXIS_NAMES)
    if "rotvec" in representations:
        active_axes.update(ROT_VECTOR_AXIS_NAMES)
    return {
        axis: bounds
        for axis, bounds in axis_limits.items()
        if axis not in active_axes
    }


def _json_axis_limits(axis_limits: Mapping[str, tuple[float, float]]) -> dict[str, list[float]]:
    return {axis: [float(bounds[0]), float(bounds[1])] for axis, bounds in axis_limits.items()}


def _density_report_field(landscape: Landscape, field_name: str):
    report = landscape.density_report
    if report is None:
        return None
    return getattr(report, field_name, None)


def _display_outlier_mask(data, color_field: str) -> np.ndarray | None:
    if color_field != "sld_display" or "sld_display_is_outlier" not in data.columns:
        return None
    return np.asarray(data["sld_display_is_outlier"], dtype=bool)


def _color_scale_mode(data, color_field: str) -> str:
    outlier_mask = _display_outlier_mask(data, color_field)
    if outlier_mask is not None and outlier_mask.any():
        return "tail_jump_exclude_outliers"
    return "auto"


def _resolve_color_scale(
    values,
    *,
    display_outlier_mask: np.ndarray | None,
    visual_style: str,
    color_vmin: float | None,
    color_vmax: float | None,
    display_threshold: float | None,
) -> tuple[float | None, float | None, str, str]:
    resolved_vmin = color_vmin
    resolved_vmax = color_vmax
    vmin_source = "explicit" if color_vmin is not None else "auto"
    vmax_source = "explicit" if color_vmax is not None else "auto"
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    scale_values = finite_values
    if display_outlier_mask is not None and display_outlier_mask.shape == np.asarray(values).shape:
        raw_values = np.asarray(values, dtype=float)
        scale_mask = np.isfinite(raw_values) & ~display_outlier_mask
        scale_values = raw_values[scale_mask]
        if scale_values.size and display_outlier_mask.any() and resolved_vmax is None:
            resolved_vmax = float(np.max(scale_values))
            vmax_source = "tail_jump_display_outliers"
    if visual_style == "legacy" and resolved_vmin is None and display_threshold is not None:
        resolved_vmin = display_threshold
        vmin_source = "legacy_display_threshold"
    if visual_style == "legacy" and resolved_vmax is None:
        if scale_values.size:
            resolved_vmax = math.ceil(float(scale_values.max()) * 10.0) / 10.0
            vmax_source = "legacy_color_max_ceiling"
    return resolved_vmin, resolved_vmax, vmin_source, vmax_source


def _resolve_display_filter(
    values,
    *,
    display_top_fraction: float | None,
    display_sld_threshold: float | None,
    display_filter_mode: str | None,
    display_max_divisor: float,
) -> tuple[str, float | None]:
    if display_filter_mode is None:
        if display_sld_threshold is not None:
            display_filter_mode = "threshold"
        elif display_top_fraction is not None:
            display_filter_mode = "top_fraction"
        else:
            display_filter_mode = "none"
    if display_filter_mode == "none":
        return display_filter_mode, None
    if display_filter_mode == "top_fraction":
        return display_filter_mode, None
    if display_filter_mode == "threshold":
        return display_filter_mode, display_sld_threshold
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if not finite_values.size:
        raise ValueError("legacy_max_divisor display filter requires finite density values")
    threshold = math.ceil(float(finite_values.max()) / display_max_divisor * 10.0) / 10.0
    return "legacy_max_divisor", threshold


def _validate_range_axes(range_bounds: Mapping[str, object], representations: list[str]) -> None:
    allowed = set()
    if "euler" in representations:
        allowed.update(EULER_AXIS_NAMES)
    if "rotvec" in representations:
        allowed.update(ROT_VECTOR_AXIS_NAMES)
    unknown = sorted(axis for axis in range_bounds if axis not in allowed)
    if unknown:
        raise ValueError(f"Unsupported display range axis for representation: {unknown}")


def _expected_visualization_paths(
    output_dir: Path,
    coordinate_sources: list[str],
    representations: list[str],
    formats: Sequence[str],
    *,
    write_projection_csvs: bool,
    output_prefix: str,
    generate_histograms: bool,
    generate_axis_direction_map: bool,
    artifact_layout: str,
    display_table_filename: str | None,
) -> list[Path]:
    paths = [
        output_dir / f"{output_prefix}visualization_report.json",
        output_dir / (display_table_filename or f"{output_prefix}display_coordinates.csv"),
    ]
    for source in coordinate_sources:
        for rep in representations:
            for fmt in formats:
                paths.append(
                    output_dir
                    / _triptych_figure_name(
                        coordinate_source=source,
                        representation=rep,
                        fmt=fmt,
                        output_prefix=output_prefix,
                        artifact_layout=artifact_layout,
                    )
                )
                cloud_name = (
                    f"landscape_3d_{rep}.{fmt}"
                    if artifact_layout == "run_bundle"
                    else f"{output_prefix}{source}_{rep}_cloud_3d.{fmt}"
                )
                paths.append(output_dir / cloud_name)
                if artifact_layout == "run_bundle":
                    projections = EULER_PROJECTIONS if rep == "euler" else ROT_VECTOR_PROJECTIONS
                    for projection_name, _, _ in projections:
                        paths.append(
                            output_dir
                            / _individual_projection_name(
                                representation=rep,
                                projection_name=projection_name,
                                fmt=fmt,
                            )
                        )
            if write_projection_csvs:
                projections = EULER_PROJECTIONS if rep == "euler" else ROT_VECTOR_PROJECTIONS
                for projection_name, _, _ in projections:
                    prefix = "euler_projection" if rep == "euler" else "projection"
                    paths.append(output_dir / f"{output_prefix}{source}_{prefix}_{projection_name}.csv")
        if generate_histograms:
            for fmt in formats:
                paths.append(output_dir / f"{output_prefix}{source}_sld_histogram.{fmt}")
                paths.append(output_dir / f"{output_prefix}{source}_rotation_angle_histogram.{fmt}")
                paths.append(output_dir / f"{output_prefix}{source}_euler_histograms.{fmt}")
        if generate_axis_direction_map:
            for fmt in formats:
                paths.append(output_dir / f"{output_prefix}{source}_axis_direction_azimuth_elevation.{fmt}")
    return paths


def _ensure_visualization_paths_available(paths: list[Path], *, overwrite: bool) -> None:
    if overwrite:
        return
    for path in paths:
        if path.exists():
            raise FileExistsError(f"Output path already exists: {path}")


def _display_indices(
    data,
    coordinates: np.ndarray,
    euler_coordinates: np.ndarray,
    *,
    representations: list[str],
    display_top_fraction: float | None,
    display_sld_threshold: float | None,
    display_density_field: str,
    range_bounds: Mapping[str, tuple[float | None, float | None] | None],
) -> np.ndarray:
    mask = np.ones(len(data), dtype=bool)
    if display_sld_threshold is not None:
        mask &= np.asarray(data[display_density_field], dtype=float) >= display_sld_threshold
    selected = np.where(mask)[0]
    if display_top_fraction is not None and selected.size:
        n_keep = max(1, int(math.ceil(selected.size * display_top_fraction)))
        values = np.asarray(data.iloc[selected][display_density_field], dtype=float)
        order = np.argsort(values, kind="mergesort")
        selected = selected[order[-n_keep:]]
        selected = np.sort(selected)
    if range_bounds and selected.size:
        range_mask = np.ones(selected.size, dtype=bool)
        if "euler" in representations:
            range_mask &= _range_mask(
                euler_coordinates[selected],
                axis_names=EULER_AXIS_NAMES,
                range_bounds=range_bounds,
                wraparound=True,
            )
        if "rotvec" in representations:
            range_mask &= _range_mask(
                coordinates[selected],
                axis_names=ROT_VECTOR_AXIS_NAMES,
                range_bounds=range_bounds,
                wraparound=False,
            )
        selected = selected[range_mask]
    return np.asarray(selected, dtype=int)


def _range_mask(
    coordinates: np.ndarray,
    *,
    axis_names: tuple[str, str, str],
    range_bounds: Mapping[str, tuple[float | None, float | None] | None],
    wraparound: bool,
) -> np.ndarray:
    mask = np.ones(coordinates.shape[0], dtype=bool)
    for axis_index, axis_name in enumerate(axis_names):
        bounds = range_bounds.get(axis_name)
        if bounds is None:
            continue
        lower, upper = bounds
        if lower is None and upper is None:
            continue
        values = coordinates[:, axis_index]
        axis_mask = np.ones(coordinates.shape[0], dtype=bool)
        if lower is not None and upper is not None and wraparound and lower > upper:
            axis_mask &= (values >= lower) | (values <= upper)
        else:
            if lower is not None:
                axis_mask &= values >= lower
            if upper is not None:
                axis_mask &= values <= upper
        mask &= axis_mask
    return mask


def _display_table(
    data,
    rotvec_coordinates: np.ndarray,
    euler_coordinates: np.ndarray,
    *,
    coordinate_source: str,
    representations: list[str],
) -> pd.DataFrame:
    output = pd.DataFrame({"particle_key": data["particle_key"].to_numpy()})
    output["coordinate_source"] = coordinate_source
    if "rotvec" in representations:
        output["rotvec_x"] = rotvec_coordinates[:, 0]
        output["rotvec_y"] = rotvec_coordinates[:, 1]
        output["rotvec_z"] = rotvec_coordinates[:, 2]
    if "euler" in representations:
        output["euler_alpha"] = euler_coordinates[:, 0]
        output["euler_beta"] = euler_coordinates[:, 1]
        output["euler_gamma"] = euler_coordinates[:, 2]
    for field in SLD_FIELDS:
        output[field] = data[field].to_numpy()
    return output


def _coordinates_for_representation(
    representation: str,
    rotvec_coordinates: np.ndarray,
    euler_coordinates: np.ndarray,
):
    if representation == "euler":
        return euler_coordinates, EULER_PROJECTIONS, EULER_AXIS_NAMES
    return rotvec_coordinates, ROT_VECTOR_PROJECTIONS, ROT_VECTOR_AXIS_NAMES


def _write_2d_figure(
    data,
    coordinates: np.ndarray,
    output_dir: Path,
    *,
    coordinate_source: str,
    representation: str,
    projections,
    axis_names: tuple[str, str, str],
    color_field: str,
    formats: Sequence[str],
    max_points_2d: int | None,
    random_seed: int,
    output_prefix: str,
    style_options: Mapping[str, Any],
    color_vmin: float | None,
    color_vmax: float | None,
    artifact_layout: str,
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    indices = _plot_indices(
        data,
        _downsample_indices(len(data), max_points_2d, random_seed),
        color_field=color_field,
        sort_points_by_color=style_options["sort_points_by_color"],
    )
    figure, axes = plt.subplots(
        1,
        3,
        figsize=style_options["figure_size_2d"],
        constrained_layout=True,
    )
    colors = np.asarray(data.iloc[indices][color_field], dtype=float)
    scatter = None
    for axis, (projection_name, axis_a, axis_b) in zip(axes, projections):
        scatter = axis.scatter(
            coordinates[indices, axis_a],
            coordinates[indices, axis_b],
            c=colors,
            s=style_options["point_size"],
            cmap=style_options["color_map"],
            alpha=style_options["point_alpha"],
            vmin=color_vmin,
            vmax=color_vmax,
            rasterized=True,
        )
        axis.set_title(projection_name.replace("_", "-"))
        axis.set_xlabel(axis_names[axis_a])
        axis.set_ylabel(axis_names[axis_b])
        if style_options["aspect"] == "equal":
            axis.set_aspect("equal", adjustable="box")
        axis_a_limits = style_options["axis_limits"].get(axis_names[axis_a])
        axis_b_limits = style_options["axis_limits"].get(axis_names[axis_b])
        if axis_a_limits is not None:
            axis.set_xlim(axis_a_limits)
        if axis_b_limits is not None:
            axis.set_ylim(axis_b_limits)
    if scatter is not None:
        if style_options["colorbar_position"] == "bottom":
            figure.colorbar(
                scatter,
                ax=axes.ravel().tolist(),
                label=color_field,
                orientation="horizontal",
                fraction=0.08,
                pad=0.12,
            )
        else:
            figure.colorbar(scatter, ax=axes.ravel().tolist(), label=color_field)
    generated = {}
    for fmt in formats:
        path = output_dir / _triptych_figure_name(
            coordinate_source=coordinate_source,
            representation=representation,
            fmt=fmt,
            output_prefix=output_prefix,
            artifact_layout=artifact_layout,
        )
        figure.savefig(path)
        generated[f"{coordinate_source}_{representation}_projection_2d_{fmt}"] = str(path)
    plt.close(figure)
    if artifact_layout == "run_bundle":
        generated.update(
            _write_individual_projection_figures(
                data,
                coordinates,
                output_dir,
                projections=projections,
                axis_names=axis_names,
                representation=representation,
                color_field=color_field,
                formats=formats,
                indices=indices,
                colors=colors,
                style_options=style_options,
                color_vmin=color_vmin,
                color_vmax=color_vmax,
            )
        )
    return generated


def _write_3d_figure(
    data,
    coordinates: np.ndarray,
    output_dir: Path,
    *,
    coordinate_source: str,
    representation: str,
    axis_names: tuple[str, str, str],
    color_field: str,
    formats: Sequence[str],
    max_points_3d: int,
    random_seed: int,
    output_prefix: str,
    style_options: Mapping[str, Any],
    color_vmin: float | None,
    color_vmax: float | None,
    artifact_layout: str,
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    indices = _plot_indices(
        data,
        _downsample_indices(len(data), max_points_3d, random_seed),
        color_field=color_field,
        sort_points_by_color=style_options["sort_points_by_color"],
    )
    figure = plt.figure(figsize=style_options["figure_size_3d"], constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    colors = np.asarray(data.iloc[indices][color_field], dtype=float)
    scatter = axis.scatter(
        coordinates[indices, 0],
        coordinates[indices, 1],
        coordinates[indices, 2],
        c=colors,
        s=style_options["point_size"],
        cmap=style_options["color_map"],
        alpha=style_options["point_alpha"],
        vmin=color_vmin,
        vmax=color_vmax,
        rasterized=True,
    )
    axis.set_xlabel(axis_names[0])
    axis.set_ylabel(axis_names[1])
    axis.set_zlabel(axis_names[2])
    axis.set_title(f"{coordinate_source} {representation} cloud")
    if style_options["colorbar_position"] == "bottom":
        figure.colorbar(
            scatter,
            ax=axis,
            label=color_field,
            orientation="horizontal",
            fraction=0.08,
            pad=0.12,
        )
    else:
        figure.colorbar(scatter, ax=axis, label=color_field)
    generated = {}
    for fmt in formats:
        filename = (
            f"landscape_3d_{representation}.{fmt}"
            if artifact_layout == "run_bundle"
            else f"{output_prefix}{coordinate_source}_{representation}_cloud_3d.{fmt}"
        )
        path = output_dir / filename
        figure.savefig(path)
        generated[f"{coordinate_source}_{representation}_cloud_3d_{fmt}"] = str(path)
    plt.close(figure)
    return generated


def _write_individual_projection_figures(
    data,
    coordinates: np.ndarray,
    output_dir: Path,
    *,
    projections,
    axis_names: tuple[str, str, str],
    representation: str,
    color_field: str,
    formats: Sequence[str],
    indices: np.ndarray,
    colors: np.ndarray,
    style_options: Mapping[str, Any],
    color_vmin: float | None,
    color_vmax: float | None,
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    generated = {}
    for projection_name, axis_a, axis_b in projections:
        figure, axis = plt.subplots(figsize=(6, 5), constrained_layout=True)
        scatter = axis.scatter(
            coordinates[indices, axis_a],
            coordinates[indices, axis_b],
            c=colors,
            s=style_options["point_size"],
            cmap=style_options["color_map"],
            alpha=style_options["point_alpha"],
            vmin=color_vmin,
            vmax=color_vmax,
            rasterized=True,
        )
        axis.set_title(projection_name.replace("_", "-"))
        axis.set_xlabel(axis_names[axis_a])
        axis.set_ylabel(axis_names[axis_b])
        if style_options["aspect"] == "equal":
            axis.set_aspect("equal", adjustable="box")
        axis_a_limits = style_options["axis_limits"].get(axis_names[axis_a])
        axis_b_limits = style_options["axis_limits"].get(axis_names[axis_b])
        if axis_a_limits is not None:
            axis.set_xlim(axis_a_limits)
        if axis_b_limits is not None:
            axis.set_ylim(axis_b_limits)
        if style_options["colorbar_position"] == "bottom":
            figure.colorbar(
                scatter,
                ax=axis,
                label=color_field,
                orientation="horizontal",
                fraction=0.08,
                pad=0.12,
            )
        else:
            figure.colorbar(scatter, ax=axis, label=color_field)
        for fmt in formats:
            path = output_dir / _individual_projection_name(
                representation=representation,
                projection_name=projection_name,
                fmt=fmt,
            )
            figure.savefig(path)
            generated[f"{representation}_{projection_name}_{fmt}"] = str(path)
        plt.close(figure)
    return generated


def _triptych_figure_name(
    *,
    coordinate_source: str,
    representation: str,
    fmt: str,
    output_prefix: str,
    artifact_layout: str,
) -> str:
    if artifact_layout == "run_bundle":
        return f"{representation}_3view_projection.{fmt}"
    return f"{output_prefix}{coordinate_source}_{representation}_projection_2d.{fmt}"


def _individual_projection_name(
    *,
    representation: str,
    projection_name: str,
    fmt: str,
) -> str:
    return f"{representation}_{projection_name}.{fmt}"


def _write_projection_csvs(
    data,
    coordinates: np.ndarray,
    output_dir: Path,
    *,
    coordinate_source: str,
    representation: str,
    projections,
    axis_names: tuple[str, str, str],
    color_field: str,
    overwrite: bool,
    output_prefix: str,
) -> dict[str, str]:
    generated = {}
    prefix = "euler_projection" if representation == "euler" else "projection"
    for projection_name, axis_a, axis_b in projections:
        filename = f"{output_prefix}{coordinate_source}_{prefix}_{projection_name}.csv"
        path = output_dir / filename
        _write_projection_csv(
            path,
            data,
            coordinates,
            axis_a=axis_a,
            axis_b=axis_b,
            axis_names=axis_names,
            color_field=color_field,
            overwrite=overwrite,
        )
        generated[f"{coordinate_source}_{prefix}_{projection_name}_csv"] = str(path)
    return generated


def _write_projection_csv(
    path: Path,
    data,
    coordinates: np.ndarray,
    *,
    axis_a: int,
    axis_b: int,
    axis_names: tuple[str, str, str],
    color_field: str,
    overwrite: bool,
) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["particle_key", axis_names[axis_a], axis_names[axis_b], color_field])
        for row_index, particle_key in enumerate(data["particle_key"]):
            writer.writerow(
                [
                    particle_key,
                    float(coordinates[row_index, axis_a]),
                    float(coordinates[row_index, axis_b]),
                    float(data.iloc[row_index][color_field]),
                ]
            )


def _write_histograms(
    data,
    rotvec_coordinates: np.ndarray,
    euler_coordinates: np.ndarray,
    output_dir: Path,
    *,
    coordinate_source: str,
    color_field: str,
    formats: Sequence[str],
    output_prefix: str,
    style_options: Mapping[str, Any],
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    generated: dict[str, str] = {}
    weights = _percentage_weights(len(data))
    histogram_specs = (
        (
            "sld_histogram",
            np.asarray(data[color_field], dtype=float),
            f"{color_field} histogram",
            color_field,
        ),
        (
            "rotation_angle_histogram",
            np.linalg.norm(rotvec_coordinates, axis=1),
            "rotation angle histogram",
            "rotation angle",
        ),
    )
    for name, values, title, xlabel in histogram_specs:
        figure, axis = plt.subplots(figsize=(6, 4), constrained_layout=True)
        axis.hist(values[np.isfinite(values)], bins=40, weights=weights[np.isfinite(values)])
        axis.set_title(title)
        axis.set_xlabel(xlabel)
        axis.set_ylabel("particles (%)")
        for fmt in formats:
            path = output_dir / f"{output_prefix}{coordinate_source}_{name}.{fmt}"
            figure.savefig(path)
            generated[f"{coordinate_source}_{name}_{fmt}"] = str(path)
        plt.close(figure)

    figure, axes = plt.subplots(1, 3, figsize=style_options["figure_size_2d"], constrained_layout=True)
    for axis, axis_name, values in zip(axes, EULER_AXIS_NAMES, euler_coordinates.T):
        finite = np.isfinite(values)
        axis.hist(values[finite], bins=40, weights=weights[finite])
        axis.set_title(axis_name)
        axis.set_xlabel(axis_name)
        axis.set_ylabel("particles (%)")
        axis_limits = style_options["axis_limits"].get(axis_name)
        if axis_limits is not None:
            axis.set_xlim(axis_limits)
    for fmt in formats:
        path = output_dir / f"{output_prefix}{coordinate_source}_euler_histograms.{fmt}"
        figure.savefig(path)
        generated[f"{coordinate_source}_euler_histograms_{fmt}"] = str(path)
    plt.close(figure)
    return generated


def _write_axis_direction_map(
    data,
    rotvec_coordinates: np.ndarray,
    output_dir: Path,
    *,
    coordinate_source: str,
    color_field: str,
    formats: Sequence[str],
    output_prefix: str,
    style_options: Mapping[str, Any],
    color_vmin: float | None,
    color_vmax: float | None,
) -> dict[str, str]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    norms = np.linalg.norm(rotvec_coordinates, axis=1)
    directions = np.zeros_like(rotvec_coordinates, dtype=float)
    nonzero = norms > 0
    directions[nonzero] = rotvec_coordinates[nonzero] / norms[nonzero, None]
    azimuth = np.degrees(np.arctan2(directions[:, 1], directions[:, 0]))
    elevation = np.degrees(np.arcsin(np.clip(directions[:, 2], -1.0, 1.0)))
    indices = _plot_indices(
        data,
        np.arange(len(data), dtype=int),
        color_field=color_field,
        sort_points_by_color=style_options["sort_points_by_color"],
    )
    figure, axis = plt.subplots(figsize=(8, 5), constrained_layout=True)
    scatter = axis.scatter(
        azimuth[indices],
        elevation[indices],
        c=np.asarray(data.iloc[indices][color_field], dtype=float),
        s=style_options["point_size"],
        cmap=style_options["color_map"],
        alpha=style_options["point_alpha"],
        vmin=color_vmin,
        vmax=color_vmax,
        rasterized=True,
    )
    axis.set_title("axis direction azimuth-elevation")
    axis.set_xlabel("azimuth")
    axis.set_ylabel("elevation")
    axis.set_xlim((-180.0, 180.0))
    axis.set_ylim((-90.0, 90.0))
    if style_options["colorbar_position"] == "bottom":
        figure.colorbar(
            scatter,
            ax=axis,
            label=color_field,
            orientation="horizontal",
            fraction=0.08,
            pad=0.12,
        )
    else:
        figure.colorbar(scatter, ax=axis, label=color_field)
    generated = {}
    for fmt in formats:
        path = output_dir / f"{output_prefix}{coordinate_source}_axis_direction_azimuth_elevation.{fmt}"
        figure.savefig(path)
        generated[f"{coordinate_source}_axis_direction_azimuth_elevation_{fmt}"] = str(path)
    plt.close(figure)
    return generated


def _percentage_weights(n_points: int) -> np.ndarray:
    if n_points == 0:
        return np.asarray([], dtype=float)
    return np.ones(n_points, dtype=float) * (100.0 / n_points)


def _plot_indices(
    data,
    indices: np.ndarray,
    *,
    color_field: str,
    sort_points_by_color: str,
) -> np.ndarray:
    if sort_points_by_color == "none" or indices.size == 0:
        return indices
    values = np.asarray(data.iloc[indices][color_field], dtype=float)
    order = np.argsort(values, kind="mergesort")
    if sort_points_by_color == "descending":
        order = order[::-1]
    return indices[order]


def _downsample_indices(n_points: int, max_points: int | None, random_seed: int) -> np.ndarray:
    if max_points is None or n_points <= max_points:
        return np.arange(n_points, dtype=int)
    rng = np.random.default_rng(random_seed)
    return np.sort(rng.choice(n_points, size=max_points, replace=False))


def _stack_coordinates(column, column_name: str) -> np.ndarray:
    coordinates = np.vstack([np.asarray(value, dtype=float) for value in column])
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(f"{column_name} must contain length-3 coordinates")
    return coordinates
