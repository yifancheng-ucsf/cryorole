"""Display-only visualization workflow service."""

from __future__ import annotations

import csv
from dataclasses import dataclass, fields, replace
import json
from pathlib import Path
import shutil
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from cryorole.core.display_policy import (
    ColorScale,
    EULER_AXIS_NAMES,
    coordinate_range_mask,
    resolve_display_indices,
)
from cryorole.core.euler_conventions import DEFAULT_EULER_CONVENTION, LEGACY_MISSING_EULER_CONVENTION_SOURCE, RAW_EULER_ANGLE_COLUMNS, resolve_euler_convention
from cryorole.export import read_landscape, read_selection_json, write_json_artifact, write_landscape_visualizations
from cryorole.io.landscape_resolver import resolve_landscape_path
from cryorole.io.writers.landscape_store import landscape_from_arrays, read_landscape_npz_arrays
from cryorole.models.landscape import Landscape
from cryorole.models.landscape_arrays import LandscapeArrays
from cryorole.run_bundle import validate_completed_run_bundle
from cryorole.logs import warn_user


@dataclass(frozen=True)
class VisualizationRequest:
    run_dir: str
    space: str = "raw"
    canonical_id: str = "default"
    selection_id: str | None = None
    use_selected_landscape: bool = False
    visual_id: str = "default"
    view: str | None = None
    representation: str = "both"
    colormap: str = "rainbow_r"
    range_bound: tuple[Any, ...] = ()
    axis_limit: tuple[Any, ...] = ()
    top_fraction: float | None = None
    threshold: float | None = None
    all_particles: bool = False
    formats: tuple[str, ...] | list[str] | str | None = None
    vmin: float | None = None
    vmax: float | None = None
    point_size: float | None = None
    alpha: float | None = None
    max_points: int | None = None
    bins: str | int = "auto"
    hist_mode: str = "percent"
    kde: bool = False
    kde_bandwidth: str | float | None = None
    three_d_mode: str = "interactive"
    overwrite: bool = False
    coordinate_source: str = "analysis"
    euler_convention: str | None = None
    # How omitted --run-dir / --canonical-id / --selection-id were filled in.
    resolved_by: dict[str, dict[str, str]] | None = None

    @classmethod
    def from_namespace(cls, namespace: Any) -> "VisualizationRequest":
        values = vars(namespace)
        return cls(**{field.name: values[field.name] for field in fields(cls) if field.name in values})


@dataclass(frozen=True)
class VisualizationResult:
    output_dir: Path
    report: dict[str, Any]


@dataclass(frozen=True)
class QuickLookRequest:
    """Typed compact rendering request used by upstream workflow services."""

    output_dir: str | Path
    euler_convention: str
    euler_convention_source: str
    color_field: str = "sld_raw"
    display_density_field: str = "sld_raw"
    color_map: str = "rainbow_r"
    color_vmin: float | None = None
    color_vmax: float | None = None
    max_points_2d: int | None = 500_000
    random_seed: int = 0
    overwrite: bool = False
    point_size: float | None = None
    point_alpha: float | None = None
    figure_width: float | None = None
    figure_height: float | None = None
    colorbar_position: str | None = None
    sort_points_by_color: str | None = None
    axis_limits: dict[str, tuple[float | None, float | None] | None] | None = None
    tail_jump_threshold: float | None = None


@dataclass(frozen=True)
class QuickLookResult:
    output_dir: Path
    report: dict[str, Any]


def write_quicklook(
    landscape: Landscape | LandscapeArrays,
    request: QuickLookRequest,
) -> QuickLookResult:
    """Write the compact, flat run quick-look from full landscape values."""

    output_dir = Path(request.output_dir)
    sld_raw = _quicklook_sld_raw(landscape)
    subset_indices, top_cutoff = _quicklook_subset_indices(sld_raw)
    generated_files: dict[str, str] = {}
    subset_counts: dict[str, int] = {}
    rendered_counts: dict[str, int] = {}
    base_report: dict[str, Any] | None = None

    for subset_name, indices in subset_indices.items():
        rendered = _deterministic_subset_sample(
            indices,
            max_points=request.max_points_2d,
            random_seed=request.random_seed,
        )
        subset_landscape = _quicklook_subset_landscape(landscape, rendered)
        subset_report = write_landscape_visualizations(
            subset_landscape,
            output_dir,
            overwrite=request.overwrite,
            coordinate_source="analysis",
            representation="both",
            color_field=request.color_field,
            euler_convention=request.euler_convention,
            euler_convention_source=request.euler_convention_source,
            euler_degrees=True,
            display_density_field=request.display_density_field,
            formats=("png",),
            max_points_2d=None,
            max_points_3d=1,
            random_seed=request.random_seed,
            write_projection_csvs=False,
            output_prefix=f"{subset_name}_",
            artifact_layout="run_bundle",
            visual_style="legacy",
            color_map=request.color_map,
            color_vmin=request.color_vmin,
            color_vmax=request.color_vmax,
            point_size=request.point_size,
            point_alpha=request.point_alpha,
            figure_width=request.figure_width,
            figure_height=request.figure_height,
            colorbar_position=request.colorbar_position,
            sort_points_by_color=request.sort_points_by_color,
            axis_limits=request.axis_limits,
            display_filter_mode="none",
            output_profile="quicklook",
        )
        if base_report is None:
            base_report = subset_report
        subset_counts[subset_name] = int(indices.size)
        rendered_counts[subset_name] = int(rendered.size)
        for key, path in subset_report["generated_files"].items():
            generated_files[f"{subset_name}_{key}"] = path

    distribution_path = output_dir / "sld_log_distribution.png"
    distribution = _write_sld_log_distribution(
        sld_raw,
        distribution_path,
        tail_jump_threshold=request.tail_jump_threshold,
        overwrite=request.overwrite,
    )
    generated_files["sld_log_distribution_png"] = str(distribution_path)
    report = dict(base_report or {})
    report.update(
        {
            "artifact_type": "run_quicklook",
            "output_profile": "quicklook",
            "report_path": None,
            "generated_files": generated_files,
            "generated_filenames": [Path(path).name for path in generated_files.values()],
            "subset_counts": subset_counts,
            "rendered_subset_counts": rendered_counts,
            "top_40pct_cutoff_sld": top_cutoff,
            "top_40pct_tie_policy": "include_all_rows_at_or_above_cutoff",
            "sld_distribution": distribution,
            "sampling_method": "deterministic_random_without_replacement",
        }
    )
    return QuickLookResult(output_dir=output_dir, report=report)


def _quicklook_sld_raw(landscape: Landscape | LandscapeArrays) -> np.ndarray:
    values = (
        landscape.sld_raw
        if isinstance(landscape, LandscapeArrays)
        else landscape.data["sld_raw"].to_numpy(dtype=float)
    )
    return np.asarray(values, dtype=float)


def _quicklook_subset_indices(
    sld_raw: np.ndarray,
) -> tuple[dict[str, np.ndarray], float]:
    values = np.asarray(sld_raw, dtype=float)
    finite = np.flatnonzero(np.isfinite(values))
    if finite.size == 0:
        raise ValueError("quick-look requires at least one finite sld_raw value")
    count = max(1, int(np.ceil(finite.size * 0.40)))
    ranked = np.sort(values[finite], kind="mergesort")
    cutoff = float(ranked[-count])
    return (
        {
            "all": np.arange(values.size, dtype=int),
            "sld_ge_1": np.flatnonzero(np.isfinite(values) & (values >= 1.0)),
            "top_40pct": np.flatnonzero(np.isfinite(values) & (values >= cutoff)),
        },
        cutoff,
    )


def _deterministic_subset_sample(
    indices: np.ndarray,
    *,
    max_points: int | None,
    random_seed: int,
) -> np.ndarray:
    selected = np.asarray(indices, dtype=int)
    if max_points is None or selected.size <= max_points:
        return selected
    sampled = np.random.default_rng(random_seed).choice(
        selected,
        size=max_points,
        replace=False,
    )
    return np.sort(sampled)


def _quicklook_subset_landscape(
    landscape: Landscape | LandscapeArrays,
    indices: np.ndarray,
) -> Landscape:
    if isinstance(landscape, LandscapeArrays):
        return landscape_from_arrays(landscape, indices)
    return Landscape(
        data=landscape.data.iloc[indices].copy(deep=True).reset_index(drop=True),
        canonical_transform=landscape.canonical_transform,
        active_policies=landscape.active_policies,
        density_report=(
            landscape.density_report if len(indices) == len(landscape.data) else None
        ),
        canonicalization_report=landscape.canonicalization_report,
    )


def _write_sld_log_distribution(
    sld_raw: np.ndarray,
    output_path: Path,
    *,
    tail_jump_threshold: float | None,
    overwrite: bool,
) -> dict[str, Any]:
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {output_path}")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    values = np.asarray(sld_raw, dtype=float)
    valid = np.isfinite(values) & (values > 0.0)
    positive = values[valid]
    if positive.size == 0:
        raise ValueError("SLD log distribution requires a positive finite sld_raw value")
    log_values = np.log10(positive)
    p99 = float(np.percentile(positive, 99.0))
    bin_count = min(72, max(10, int(np.ceil(np.sqrt(positive.size)))))
    figure, axis = plt.subplots(figsize=(9.0, 5.0), constrained_layout=True)
    weights = np.full(positive.size, 100.0 / positive.size)
    axis.hist(log_values, bins=bin_count, weights=weights, color="#4C78A8", alpha=0.85)
    axis.set_xlabel("log10(SLD raw)")
    axis.set_ylabel("Particles (%)")
    axis.set_title("SLD distribution (full landscape)")
    markers = (
        (1.0, "SLD = 1", "#444444"),
        (100.0, "SLD = 100", "#E45756"),
        (p99, "P99", "#54A24B"),
    )
    lower, upper = float(log_values.min()), float(log_values.max())
    outside_markers = []
    for value, label, color in markers:
        location = float(np.log10(value))
        if lower <= location <= upper:
            axis.axvline(location, color=color, linewidth=1.4, linestyle="--", label=label)
        else:
            outside_markers.append(f"{label} outside observed range")
    if tail_jump_threshold is not None and np.isfinite(tail_jump_threshold) and tail_jump_threshold > 0:
        location = float(np.log10(tail_jump_threshold))
        if lower <= location <= upper:
            axis.axvline(location, color="#B279A2", linewidth=1.6, linestyle=":", label="Tail jump")
        else:
            outside_markers.append("Tail jump outside observed range")
    handles, labels = axis.get_legend_handles_labels()
    if handles:
        axis.legend(handles, labels, frameon=False)
    if outside_markers:
        axis.text(
            0.99,
            0.98,
            "\n".join(outside_markers),
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=8,
        )
    axis.grid(axis="y", alpha=0.2)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=160)
    plt.close(figure)
    return {
        "input_count": int(values.size),
        "positive_finite_count": int(positive.size),
        "excluded_nonfinite_count": int((~np.isfinite(values)).sum()),
        "excluded_nonpositive_count": int((np.isfinite(values) & (values <= 0.0)).sum()),
        "p99_sld_raw": p99,
        "bin_count": bin_count,
        "tail_jump_threshold": tail_jump_threshold,
        "markers": ["sld_1", "sld_100", "p99"] + (["tail_jump"] if tail_jump_threshold is not None else []),
        "source": "full_landscape_sld_raw",
    }


def visualize(request: VisualizationRequest) -> VisualizationResult:
    args = request
    if not args.run_dir:
        raise ValueError("visualize requires --run-dir")
    validate_completed_run_bundle(args.run_dir)
    if args.use_selected_landscape and not args.selection_id:
        raise ValueError("visualize --use-selected-landscape requires --selection-id")
    if args.all_particles and (args.top_fraction is not None or args.threshold is not None):
        raise ValueError("Use only one display density filter: --all, --top-fraction, or --sld-threshold")
    if args.top_fraction is not None and args.threshold is not None:
        raise ValueError("Use only one display density filter: --top-fraction or --threshold")
    if args.top_fraction is not None and not 0.0 < args.top_fraction <= 1.0:
        raise ValueError("--top-fraction must be > 0 and <= 1")
    if args.threshold is not None and not np.isfinite(args.threshold):
        raise ValueError("--sld-threshold must be finite")
    if args.point_size is not None and (not np.isfinite(args.point_size) or args.point_size <= 0):
        raise ValueError("--point-size must be positive and finite")
    if args.alpha is not None and (not np.isfinite(args.alpha) or not 0.0 <= args.alpha <= 1.0):
        raise ValueError("--opacity must be finite and in [0, 1]")
    if args.max_points is not None and (
        isinstance(args.max_points, bool) or int(args.max_points) != args.max_points or args.max_points <= 0
    ):
        raise ValueError("--max-points must be a positive integer")
    if args.hist_mode not in {"count", "percent"}:
        raise ValueError("--hist-mode must be count or percent")
    if args.three_d_mode not in {"interactive", "static"}:
        raise ValueError("--3d-mode must be interactive or static")
    views = _parse_views(args.view)
    threshold = (
        None
        if args.all_particles or args.top_fraction is not None
        else (1.0 if args.threshold is None else float(args.threshold))
    )
    if args.kde_bandwidth is not None and not args.kde:
        raise ValueError("--kde-bandwidth requires --kde")
    kde_bandwidth = _parse_kde_bandwidth(args.kde_bandwidth) if args.kde else None
    bins = _parse_histogram_bins(args.bins)
    formats = _parse_formats(args.formats)
    _validate_visualize_color_field_args(args)
    axis_limits = dict(args.axis_limit or ())
    _validate_axis_limit_mapping(axis_limits)
    selected_landscape_report = None
    selected_landscape_path = None
    if args.use_selected_landscape:
        selected_landscape_path = _resolve_selected_landscape_path(
            Path(args.run_dir),
            str(args.selection_id),
        )
        selected_landscape_report = _read_selected_landscape_report(
            Path(args.run_dir),
            str(args.selection_id),
        )
    effective_space = _effective_visualization_space(
        args,
        selected_landscape_report=selected_landscape_report,
    )
    effective_args = replace(
        args,
        space=effective_space,
        threshold=threshold,
        view=",".join(views),
        formats=formats,
        bins=bins,
        kde_bandwidth=kde_bandwidth,
    )
    landscape_source = args.run_dir
    output_dir = _visualization_output_dir(effective_args)
    euler_metadata = _resolve_visualization_euler_metadata(
        effective_args,
        selected_landscape_report=selected_landscape_report,
    )
    selection_metadata = None
    array_native_stats = None
    resolved_color_scale = None
    if args.use_selected_landscape:
        source_landscape_path = selected_landscape_path
    else:
        source_landscape_path = resolve_landscape_path(
            landscape_source,
            space=effective_args.space,
            canonical_id=effective_args.canonical_id,
        )
    range_bounds = dict(args.range_bound or ())
    _validate_range_mapping(range_bounds)
    if Path(source_landscape_path).suffix.casefold() == ".npz":
        (
            landscape,
            selection_metadata,
            array_native_stats,
            resolved_color_scale,
        ) = _prepare_array_native_visualization(
            effective_args,
            Path(source_landscape_path),
            euler_metadata=euler_metadata,
            range_bounds=range_bounds,
            selected_landscape_report=selected_landscape_report,
        )
    else:
        landscape = read_landscape(
            source_landscape_path,
        )
        if args.use_selected_landscape:
            selection_metadata = _selected_landscape_visualization_metadata(
                Path(args.run_dir),
                args.selection_id,
                Path(source_landscape_path),
                selected_count=len(landscape.data),
                report=selected_landscape_report,
            )
        elif args.selection_id:
            landscape, selection_metadata = _filter_landscape_by_selection(effective_args, landscape)
        landscape, legacy_stats, resolved_color_scale = _prepare_legacy_visualization(
            effective_args,
            landscape,
            euler_metadata=euler_metadata,
            range_bounds=range_bounds,
            views=views,
        )
        array_native_stats = legacy_stats
    _validate_visualize_color_field(landscape, "sld_display")
    _prepare_visualization_output_dir(Path(output_dir), overwrite=args.overwrite)
    generated_files: dict[str, str] = {}
    view_reports: dict[str, Any] = {}
    report: dict[str, Any] = {}
    if "2d" in views:
        report = write_landscape_visualizations(
            landscape,
            output_dir,
            overwrite=False,
            coordinate_source=_coordinate_source_for_landscape(effective_args),
            representation=args.representation,
            color_field="sld_display",
            euler_convention=str(euler_metadata["euler_convention"]),
            euler_convention_source=str(euler_metadata["euler_convention_source"]),
            euler_degrees=True,
            formats=formats,
            max_points_2d=_max_points_for_view(effective_args, "2d"),
            max_points_3d=1,
            random_seed=0,
            write_projection_csvs=False,
            artifact_layout="run_bundle",
            selection_metadata=selection_metadata,
            visual_style="legacy",
            color_map=args.colormap,
            color_vmin=resolved_color_scale.vmin,
            color_vmax=resolved_color_scale.vmax,
            point_size=args.point_size,
            point_alpha=args.alpha,
            axis_limits=axis_limits or None,
            output_profile="quicklook",
        )
        generated_files.update(report.get("generated_files") or {})
        view_reports["2d"] = {
            "rendered_count": _view_count(array_native_stats, "2d"),
            "formats": list(formats),
        }
    if "1d" in views:
        one_d_report = _write_combined_1d_distributions(
            landscape,
            output_dir=Path(output_dir),
            coordinate_source=_coordinate_source_for_landscape(effective_args),
            representation=args.representation,
            euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
            formats=formats,
            bins=bins,
            hist_mode=args.hist_mode,
            kde=args.kde,
            kde_bandwidth=kde_bandwidth,
            axis_limits=axis_limits,
        )
        generated_files.update(one_d_report["generated_files"])
        view_reports["1d"] = one_d_report
    if "3d" in views:
        if args.three_d_mode == "interactive":
            three_d_report = _write_interactive_3d(
                landscape,
                output_dir=Path(output_dir),
                coordinate_source=_coordinate_source_for_landscape(effective_args),
                representation=args.representation,
                euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
                colormap=args.colormap,
                color_vmin=resolved_color_scale.vmin,
                color_vmax=resolved_color_scale.vmax,
                max_points=_max_points_for_view(effective_args, "3d"),
                random_seed=0,
            )
        else:
            three_d_report = _write_static_3d(
                landscape,
                output_dir=Path(output_dir),
                coordinate_source=_coordinate_source_for_landscape(effective_args),
                representation=args.representation,
                euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
                formats=formats,
                colormap=args.colormap,
                color_vmin=resolved_color_scale.vmin,
                color_vmax=resolved_color_scale.vmax,
                point_size=args.point_size,
                alpha=args.alpha,
                axis_limits=axis_limits,
                max_points=_max_points_for_view(effective_args, "3d"),
                random_seed=0,
            )
        generated_files.update(three_d_report["generated_files"])
        view_reports["3d"] = three_d_report
    report = _visualization_public_report(
        report,
        args=effective_args,
        source_landscape_path=source_landscape_path,
        euler_metadata=euler_metadata,
        formats=formats,
        requested_views=views,
        view_reports=view_reports,
        generated_files=generated_files,
        preparation_stats=array_native_stats,
        resolved_color_scale=resolved_color_scale,
        axis_limits=axis_limits,
        selection_metadata=selection_metadata,
    )
    write_json_artifact(
        report,
        Path(output_dir) / "visualization_report.json",
        overwrite=True,
    )
    return VisualizationResult(output_dir=Path(output_dir), report=report)


def _prepare_array_native_visualization(
    args: VisualizationRequest,
    source_path: Path,
    *,
    euler_metadata: dict[str, object],
    range_bounds: dict[str, tuple[float | None, float | None] | None],
    selected_landscape_report: dict[str, object] | None,
) -> tuple[Landscape, dict[str, object] | None, dict[str, object], Any]:
    """Filter compact arrays before constructing the plotting DataFrame."""

    arrays = read_landscape_npz_arrays(source_path)
    full_input_count = arrays.n_points
    candidate_indices = np.arange(arrays.n_points, dtype=int)
    selection_metadata = None
    if args.use_selected_landscape:
        selection_metadata = _selected_landscape_visualization_metadata(
            Path(args.run_dir),
            str(args.selection_id),
            source_path,
            selected_count=arrays.n_points,
            report=selected_landscape_report,
        )
    elif args.selection_id:
        selection_keys, selection_metadata = _read_visualization_selection_keys(
            Path(args.run_dir),
            args.selection_id,
            landscape_count=arrays.n_points,
        )
        selected_set = {str(value) for value in selection_keys}
        available_set = set(arrays.particle_key.tolist())
        missing = sorted(selected_set.difference(available_set))
        if missing:
            preview = ", ".join(missing[:5])
            suffix = "" if len(missing) <= 5 else f", ... ({len(missing)} total)"
            raise ValueError(
                "Selection contains particle_key values absent from the requested landscape: "
                f"{preview}{suffix}"
            )
        candidate_indices = candidate_indices[
            np.fromiter(
                (key in selected_set for key in arrays.particle_key),
                dtype=bool,
                count=arrays.n_points,
            )
        ]
        if not candidate_indices.size:
            raise ValueError("Selection filtering produced an empty landscape")
        selection_metadata.update(
            {
                "selection_filter_applied": True,
                "selected_landscape_used": False,
                "landscape_count_before_selection_filter": arrays.n_points,
                "landscape_count_after_selection_filter": int(candidate_indices.size),
                "missing_selected_key_count": 0,
            }
        )
    scientific_input_count = int(candidate_indices.size)
    density = arrays.sld_display[candidate_indices]
    relative_display_indices = resolve_display_indices(
        density,
        threshold=args.threshold,
        top_fraction=args.top_fraction,
    )
    display_indices = candidate_indices[relative_display_indices]
    coordinates = _array_visualization_coordinates(arrays, args.space)
    if range_bounds and display_indices.size:
        range_mask = np.ones(display_indices.size, dtype=bool)
        if args.representation in {"both", "euler"}:
            euler = Rotation.from_rotvec(coordinates[display_indices]).as_euler(
                str(euler_metadata["scipy_euler_sequence"]),
                degrees=True,
            )
            range_mask &= coordinate_range_mask(
                euler,
                axis_names=EULER_AXIS_NAMES,
                range_bounds=range_bounds,
                wraparound=True,
            )
        if args.representation in {"both", "rotvec"}:
            range_mask &= coordinate_range_mask(
                coordinates[display_indices],
                axis_names=("x", "y", "z"),
                range_bounds=range_bounds,
                wraparound=False,
            )
        display_indices = display_indices[range_mask]
    if not display_indices.size:
        raise ValueError("Display filtering produced an empty landscape")
    color_scale = _resolve_visualize_color_scale(
        arrays.sld_display[candidate_indices],
        display_outlier_mask=arrays.sld_display_is_outlier[candidate_indices],
        color_vmin=args.vmin,
        color_vmax=args.vmax,
        display_threshold=args.threshold,
    )
    views = _parse_views(args.view)
    materialized_indices = display_indices
    if "1d" not in views:
        requested_caps = [_max_points_for_view(args, view) for view in views]
        materialized_indices = _deterministic_subset_sample(
            display_indices,
            max_points=max(requested_caps) if requested_caps else None,
            random_seed=0,
        )
    landscape = landscape_from_arrays(arrays, materialized_indices)
    source = "canonical" if args.space == "canonical" else "analysis"
    filtered_count = int(display_indices.size)
    count_2d = min(filtered_count, _max_points_for_view(args, "2d")) if "2d" in views else 0
    count_3d = min(filtered_count, _max_points_for_view(args, "3d")) if "3d" in views else 0
    stats: dict[str, object] = {
        "array_native_preparation": True,
        "n_points_full_parent": int(full_input_count),
        "n_points_input": scientific_input_count,
        "n_points_after_display_filter": {source: filtered_count},
        "n_points_1d": {source: filtered_count if "1d" in views else 0},
        "n_points_2d": {source: count_2d},
        "n_points_3d": {source: count_3d},
        "materialized_plotting_rows": int(materialized_indices.size),
        "downsampled": {
            source: {
                "2d": bool("2d" in views and count_2d < filtered_count),
                "3d": bool("3d" in views and count_3d < filtered_count),
            }
        },
        "sampling_method": "deterministic_random_without_replacement",
        "sampling_policy": {
            "filter_scope": "full_parent_arrays",
            "max_points_2d": _max_points_for_view(args, "2d"),
            "max_points_3d": _max_points_for_view(args, "3d"),
            "random_seed": 0,
            "one_d_uses_full_filtered_rows": True,
            "independent_2d_3d_limits": True,
        },
        "n_sld_display_outliers": int(
            arrays.sld_display_is_outlier[candidate_indices].sum()
        ),
        "fraction_sld_display_outliers": float(
            arrays.sld_display_is_outlier[candidate_indices].sum()
            / scientific_input_count
        ),
    }
    return landscape, selection_metadata, stats, color_scale


def _array_visualization_coordinates(arrays: LandscapeArrays, space: str) -> np.ndarray:
    if space == "canonical":
        if arrays.coordinates_canonical is None:
            raise ValueError("Canonical landscape is missing coordinates_canonical")
        return arrays.coordinates_canonical
    return arrays.coordinates_analysis


def _prepare_legacy_visualization(
    args: VisualizationRequest,
    landscape: Landscape,
    *,
    euler_metadata: dict[str, object],
    range_bounds: dict[str, tuple[float | None, float | None] | None],
    views: tuple[str, ...],
) -> tuple[Landscape, dict[str, object], ColorScale]:
    """Apply display filters before materializing legacy plotting rows."""

    data = landscape.data
    coordinates = _landscape_coordinates(landscape, _coordinate_source_for_landscape(args))
    density = data["sld_display"].to_numpy(dtype=float)
    display_indices = resolve_display_indices(
        density,
        threshold=args.threshold,
        top_fraction=args.top_fraction,
    )
    if range_bounds and display_indices.size:
        mask = np.ones(display_indices.size, dtype=bool)
        if args.representation in {"both", "euler"}:
            euler = Rotation.from_rotvec(coordinates[display_indices]).as_euler(
                str(euler_metadata["scipy_euler_sequence"]), degrees=True
            )
            mask &= coordinate_range_mask(
                euler,
                axis_names=EULER_AXIS_NAMES,
                range_bounds=range_bounds,
                wraparound=True,
            )
        if args.representation in {"both", "rotvec"}:
            mask &= coordinate_range_mask(
                coordinates[display_indices],
                axis_names=("x", "y", "z"),
                range_bounds=range_bounds,
                wraparound=False,
            )
        display_indices = display_indices[mask]
    if not display_indices.size:
        raise ValueError("Display filtering produced an empty landscape")
    materialized = display_indices
    if "1d" not in views:
        cap = max(_max_points_for_view(args, view) for view in views)
        materialized = _deterministic_subset_sample(display_indices, max_points=cap, random_seed=0)
    subset = Landscape(
        data=data.iloc[materialized].copy(deep=True).reset_index(drop=True),
        canonical_transform=landscape.canonical_transform,
        active_policies=landscape.active_policies,
        density_report=landscape.density_report,
        canonicalization_report=landscape.canonicalization_report,
    )
    source = "canonical" if args.space == "canonical" else "analysis"
    outlier = (
        data["sld_display_is_outlier"].to_numpy(dtype=bool)
        if "sld_display_is_outlier" in data.columns
        else np.zeros(len(data), dtype=bool)
    )
    filtered_count = int(display_indices.size)
    stats = {
        "array_native_preparation": False,
        "n_points_full_parent": int(len(data)),
        "n_points_input": int(len(data)),
        "n_points_after_display_filter": {source: filtered_count},
        "n_points_1d": {source: filtered_count if "1d" in views else 0},
        "n_points_2d": {
            source: min(filtered_count, _max_points_for_view(args, "2d")) if "2d" in views else 0
        },
        "n_points_3d": {
            source: min(filtered_count, _max_points_for_view(args, "3d")) if "3d" in views else 0
        },
        "materialized_plotting_rows": int(len(materialized)),
        "n_sld_display_outliers": int(outlier.sum()),
        "fraction_sld_display_outliers": float(outlier.sum() / len(data)) if len(data) else 0.0,
        "sampling_method": "deterministic_random_without_replacement",
        "sampling_policy": {
            "filter_scope": "full_parent_rows",
            "max_points_2d": _max_points_for_view(args, "2d"),
            "max_points_3d": _max_points_for_view(args, "3d"),
            "random_seed": 0,
            "one_d_uses_full_filtered_rows": True,
        },
    }
    scale = _resolve_visualize_color_scale(
        density,
        display_outlier_mask=outlier,
        color_vmin=args.vmin,
        color_vmax=args.vmax,
        display_threshold=args.threshold,
    )
    return subset, stats, scale


def _write_combined_1d_distributions(
    landscape: Landscape,
    *,
    output_dir: Path,
    coordinate_source: str,
    representation: str,
    euler_sequence: str,
    formats: tuple[str, ...],
    bins: str | int,
    hist_mode: str,
    kde: bool,
    kde_bandwidth: str | float | None,
    axis_limits: dict[str, tuple[float | None, float | None] | None],
) -> dict[str, Any]:
    """Write one flat three-panel marginal-distribution figure per representation."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    coordinates = _landscape_coordinates(landscape, coordinate_source)
    generated: dict[str, str] = {}
    statistics: dict[str, dict[str, Any]] = {}
    resolved_bins: dict[str, int] = {}
    warnings: list[str] = []
    for rep in _representations(representation):
        if rep == "euler":
            values = Rotation.from_rotvec(coordinates).as_euler(euler_sequence, degrees=True)
            names = ("alpha", "beta", "gamma")
            labels = ("alpha (deg)", "beta (deg)", "gamma (deg)")
            if kde:
                warnings.append(
                    "Euler KDE is a coordinate KDE; it does not correct periodic boundaries or gimbal singularity."
                )
        else:
            values = coordinates
            names = ("x", "y", "z")
            labels = ("rotvec x (rad)", "rotvec y (rad)", "rotvec z (rad)")
        figure, axes = plt.subplots(1, 3, figsize=(18.0, 5.0), constrained_layout=True)
        for index, (axis, name, label) in enumerate(zip(axes, names, labels)):
            finite = values[:, index]
            finite = finite[np.isfinite(finite)]
            bin_count = _resolved_bin_count(finite, bins)
            resolved_bins[f"{rep}_{name}"] = bin_count
            weights = np.full(finite.size, 100.0 / finite.size) if hist_mode == "percent" else None
            counts, edges, _ = axis.hist(
                finite,
                bins=bin_count,
                weights=weights,
                color="silver",
                edgecolor="black",
                alpha=0.65,
            )
            kde_status = "disabled"
            if kde:
                kde_status = _draw_coordinate_kde(
                    axis,
                    finite,
                    edges=edges,
                    hist_mode=hist_mode,
                    bandwidth=kde_bandwidth or "scott",
                )
                if kde_status != "ok":
                    warnings.append(f"{rep}_{name}: {kde_status}")
            limit = axis_limits.get(name)
            if limit is not None and limit[0] is not None and limit[1] is not None:
                axis.set_xlim(float(limit[0]), float(limit[1]))
            axis.set_xlabel(label)
            axis.set_ylabel("Particles (%)" if hist_mode == "percent" else "Particles")
            axis.grid(True, alpha=0.25)
            statistics[f"{rep}_{name}"] = _distribution_stats(finite, kde_status=kde_status)
        for fmt in formats:
            path = output_dir / f"{rep}_1d_distribution.{fmt}"
            figure.savefig(path)
            generated[f"{rep}_1d_distribution_{fmt}"] = str(path)
        plt.close(figure)
    return {
        "input_row_count": int(len(landscape.data)),
        "bins": bins,
        "resolved_bins": resolved_bins,
        "hist_mode": hist_mode,
        "kde": kde,
        "kde_bandwidth": kde_bandwidth,
        "coordinate_density_only": True,
        "statistics": statistics,
        "warnings": warnings,
        "generated_files": generated,
    }


def _draw_coordinate_kde(axis, values, *, edges, hist_mode: str, bandwidth) -> str:
    from scipy.stats import gaussian_kde

    try:
        if values.size < 2 or np.isclose(float(np.std(values)), 0.0):
            raise ValueError("not enough spread for KDE")
        lower, upper = float(np.min(values)), float(np.max(values))
        grid = np.linspace(lower, upper, 512)
        curve = gaussian_kde(values, bw_method=bandwidth)(grid)
        width = float(edges[1] - edges[0]) if len(edges) > 1 else 1.0
        curve *= (100.0 if hist_mode == "percent" else values.size) * width
        axis.plot(grid, curve, color="black", linewidth=1.8)
        return "ok"
    except Exception as exc:
        return f"skipped: {exc}"


def _distribution_stats(values: np.ndarray, *, kde_status: str) -> dict[str, Any]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    result: dict[str, Any] = {"n": int(finite.size), "kde_status": kde_status}
    if not finite.size:
        result.update({key: None for key in ("min", "max", "median", "q05", "q25", "q75", "q95")})
        return result
    result.update(
        {
            "min": float(np.min(finite)),
            "max": float(np.max(finite)),
            "median": float(np.median(finite)),
            "q05": float(np.quantile(finite, 0.05)),
            "q25": float(np.quantile(finite, 0.25)),
            "q75": float(np.quantile(finite, 0.75)),
            "q95": float(np.quantile(finite, 0.95)),
        }
    )
    return result


def _write_static_3d(
    landscape: Landscape,
    *,
    output_dir: Path,
    coordinate_source: str,
    representation: str,
    euler_sequence: str,
    formats: tuple[str, ...],
    colormap: str,
    color_vmin: float | None,
    color_vmax: float | None,
    point_size: float | None,
    alpha: float | None,
    axis_limits: dict[str, tuple[float | None, float | None] | None],
    max_points: int,
    random_seed: int,
) -> dict[str, Any]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    coordinates = _landscape_coordinates(landscape, coordinate_source)
    indices = _deterministic_subset_sample(
        np.arange(len(coordinates)), max_points=max_points, random_seed=random_seed
    )
    density = landscape.data["sld_display"].to_numpy(dtype=float)[indices]
    generated: dict[str, str] = {}
    for rep in _representations(representation):
        rep_values = (
            Rotation.from_rotvec(coordinates[indices]).as_euler(euler_sequence, degrees=True)
            if rep == "euler"
            else coordinates[indices]
        )
        names = ("alpha", "beta", "gamma") if rep == "euler" else ("x", "y", "z")
        figure = plt.figure(figsize=(8.0, 7.0), constrained_layout=True)
        axis = figure.add_subplot(111, projection="3d")
        scatter = axis.scatter(
            rep_values[:, 0], rep_values[:, 1], rep_values[:, 2],
            c=density, cmap=colormap, vmin=color_vmin, vmax=color_vmax,
            s=point_size or 1.0, alpha=1.0 if alpha is None else alpha,
            rasterized=True,
        )
        for setter, name in zip((axis.set_xlim, axis.set_ylim, axis.set_zlim), names):
            limit = axis_limits.get(name)
            if limit is not None and limit[0] is not None and limit[1] is not None:
                setter(float(limit[0]), float(limit[1]))
        axis.set_xlabel(names[0])
        axis.set_ylabel(names[1])
        axis.set_zlabel(names[2])
        spans = tuple(
            max(float(upper - lower), np.finfo(float).eps)
            for lower, upper in (axis.get_xlim3d(), axis.get_ylim3d(), axis.get_zlim3d())
        )
        axis.set_box_aspect(spans)
        figure.colorbar(scatter, ax=axis, label="SLD", orientation="horizontal", pad=0.1)
        for fmt in formats:
            path = output_dir / f"landscape_3d_{rep}.{fmt}"
            figure.savefig(path)
            generated[f"landscape_3d_{rep}_{fmt}"] = str(path)
        plt.close(figure)
    return {
        "mode": "static",
        "rendered_count": int(len(indices)),
        "max_points": max_points,
        "generated_files": generated,
        "warnings": _three_d_representation_warnings(representation),
    }


def _write_interactive_3d(
    landscape: Landscape,
    *,
    output_dir: Path,
    coordinate_source: str,
    representation: str,
    euler_sequence: str,
    colormap: str,
    color_vmin: float | None,
    color_vmax: float | None,
    max_points: int,
    random_seed: int,
) -> dict[str, Any]:
    """Write a self-contained, selection-free offline canvas 3D viewer."""

    import matplotlib
    from matplotlib.colors import Normalize, to_hex

    coordinates = _landscape_coordinates(landscape, coordinate_source)
    indices = _deterministic_subset_sample(
        np.arange(len(coordinates)), max_points=max_points, random_seed=random_seed
    )
    density = landscape.data["sld_display"].to_numpy(dtype=float)[indices]
    finite = density[np.isfinite(density)]
    lower = float(np.min(finite)) if color_vmin is None else float(color_vmin)
    upper = float(np.max(finite)) if color_vmax is None else float(color_vmax)
    if np.isclose(lower, upper):
        upper = lower + 1.0
    cmap = matplotlib.colormaps.get_cmap(colormap)
    colors = [to_hex(cmap(Normalize(lower, upper, clip=True)(value))) for value in density]
    reps: dict[str, Any] = {}
    if representation in {"both", "rotvec"}:
        reps["rotvec"] = {
            "labels": ["RV x (rad)", "RV y (rad)", "RV z (rad)"],
            "values": coordinates[indices].tolist(),
        }
    if representation in {"both", "euler"}:
        reps["euler"] = {
            "labels": ["alpha (deg)", "beta (deg)", "gamma (deg)"],
            "values": Rotation.from_rotvec(coordinates[indices]).as_euler(
                euler_sequence, degrees=True
            ).tolist(),
        }
    payload = {
        "representations": reps,
        "particle_keys": landscape.data["particle_key"].astype(str).to_numpy()[indices].tolist(),
        "sld": density.tolist(),
        "colors": colors,
        "colormap": colormap,
        "vmin": lower,
        "vmax": upper,
        "display_only": True,
        "selection_enabled": False,
    }
    payload_json = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    html = _INTERACTIVE_3D_HTML.replace("__CRYOROLE_DATA__", payload_json)
    path = output_dir / "landscape_3d.html"
    path.write_text(html, encoding="utf-8")
    return {
        "mode": "interactive_offline_html",
        "rendered_count": int(len(indices)),
        "max_points": max_points,
        "offline": True,
        "public_cdn": False,
        "selection_enabled": False,
        "generated_files": {"landscape_3d_html": str(path)},
        "warnings": _three_d_representation_warnings(representation),
    }


def _three_d_representation_warnings(representation: str) -> list[str]:
    warnings = ["Interactive 3D is a display embedding, not an undistorted representation of SO(3)."]
    if representation in {"both", "euler"}:
        warnings.append("Euler coordinates have periodic seams and gimbal singularities.")
    if representation in {"both", "rotvec"}:
        warnings.append("Rotation-vector coordinates have a representation boundary at rotation angle pi.")
    return warnings


_INTERACTIVE_3D_HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>cryoROLE interactive 3D</title><style>
body{margin:0;font:14px system-ui;background:#f7f7f8;color:#202124}header{padding:10px 16px;background:#fff;border-bottom:1px solid #ddd;display:flex;gap:14px;align-items:center}canvas{display:block;width:100vw;height:calc(100vh - 92px);background:#fff}.note{padding:6px 16px;color:#555}#tip{position:fixed;display:none;padding:6px 8px;background:#fff;border:1px solid #aaa;pointer-events:none;white-space:pre}</style></head>
<body><header><strong>cryoROLE 3D · display only</strong><label>Representation <select id="rep"></select></label><button id="reset">Reset view</button><span id="summary"></span></header>
<div class="note">Drag to rotate · wheel to zoom · hover for particle key, SLD, and coordinates. This viewer cannot create a Selection.</div><canvas id="plot"></canvas><div id="tip"></div>
<script id="data" type="application/json">__CRYOROLE_DATA__</script><script>
"use strict";const d=JSON.parse(document.getElementById("data").textContent),c=document.getElementById("plot"),x=c.getContext("2d"),sel=document.getElementById("rep"),tip=document.getElementById("tip");Object.keys(d.representations).forEach(k=>{const o=document.createElement("option");o.value=k;o.textContent=k;sel.appendChild(o)});let ax=-.35,ay=.55,z=1,drag=false,lx=0,ly=0,screen=[];function size(){c.width=innerWidth*devicePixelRatio;c.height=(innerHeight-92)*devicePixelRatio;draw()}function project(p){let X=p[0],Y=p[1],Z=p[2],cy=Math.cos(ay),sy=Math.sin(ay),cx=Math.cos(ax),sx=Math.sin(ax),u=cy*X+sy*Z,v=sx*sy*X+cx*Y-sx*cy*Z,w=-cx*sy*X+sx*Y+cx*cy*Z;return[u,v,w]}function draw(){const r=d.representations[sel.value],v=r.values;if(!v.length)return;x.clearRect(0,0,c.width,c.height);const q=v.map(project);let m=1e-9;for(const p of q)m=Math.max(m,Math.abs(p[0]),Math.abs(p[1]));const s=.42*Math.min(c.width,c.height)*z/m;screen=q.map((p,i)=>[c.width/2+p[0]*s,c.height/2-p[1]*s,p[2],i]).sort((a,b)=>a[2]-b[2]);for(const p of screen){x.globalAlpha=.78;x.fillStyle=d.colors[p[3]];x.beginPath();x.arc(p[0],p[1],2.2*devicePixelRatio,0,Math.PI*2);x.fill()}x.globalAlpha=1;document.getElementById("summary").textContent=`${v.length.toLocaleString()} displayed · ${d.colormap} · SLD ${d.vmin.toPrecision(4)}–${d.vmax.toPrecision(4)}`}
c.onpointerdown=e=>{drag=true;lx=e.clientX;ly=e.clientY;c.setPointerCapture(e.pointerId)};c.onpointermove=e=>{if(drag){ay+=(e.clientX-lx)*.008;ax+=(e.clientY-ly)*.008;lx=e.clientX;ly=e.clientY;draw();tip.style.display="none";return}const rect=c.getBoundingClientRect(),px=(e.clientX-rect.left)*c.width/rect.width,py=(e.clientY-rect.top)*c.height/rect.height;let b=null,bd=100*devicePixelRatio;for(const p of screen){const dd=Math.hypot(p[0]-px,p[1]-py);if(dd<bd){bd=dd;b=p}}if(!b){tip.style.display="none";return}const i=b[3],r=d.representations[sel.value];tip.textContent=`${d.particle_keys[i]}\nSLD ${d.sld[i]}\n${r.labels.map((n,j)=>n+" "+r.values[i][j].toPrecision(6)).join("\n")}`;tip.style.left=e.clientX+12+"px";tip.style.top=e.clientY+12+"px";tip.style.display="block"};c.onpointerup=c.onpointercancel=c.onlostpointercapture=()=>drag=false;c.onwheel=e=>{e.preventDefault();z*=Math.exp(-e.deltaY*.001);draw()};sel.onchange=draw;document.getElementById("reset").onclick=()=>{ax=-.35;ay=.55;z=1;draw()};addEventListener("resize",size);size();
</script></body></html>"""


def _parse_views(value: str | None) -> tuple[str, ...]:
    requested = ["2d"] if value is None else [part.strip().lower() for part in value.split(",")]
    views = tuple(dict.fromkeys(part for part in requested if part))
    invalid = sorted(set(views).difference({"1d", "2d", "3d"}))
    if invalid or not views:
        detail = ", ".join(invalid) if invalid else "empty view list"
        raise ValueError(f"--view supports only 2d,1d,3d: {detail}")
    return views


def _representations(value: str) -> tuple[str, ...]:
    if value == "both":
        return ("euler", "rotvec")
    if value in {"euler", "rotvec"}:
        return (value,)
    raise ValueError("--representation must be both, euler, or rotvec")


def _parse_histogram_bins(value: str | int) -> str | int:
    if isinstance(value, int):
        if value <= 0:
            raise ValueError("--bins must be auto or a positive integer")
        return value
    normalized = str(value).strip().lower()
    if normalized == "auto":
        return "auto"
    try:
        parsed = int(normalized)
    except ValueError as exc:
        raise ValueError("--bins must be auto or a positive integer") from exc
    if parsed <= 0:
        raise ValueError("--bins must be auto or a positive integer")
    return parsed


def _resolved_bin_count(values: np.ndarray, bins: str | int) -> int:
    if isinstance(bins, int):
        return bins
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size < 2 or np.isclose(float(np.ptp(finite)), 0.0):
        return 20
    edges = np.histogram_bin_edges(finite, bins="fd")
    return min(100, max(20, int(len(edges) - 1)))


def _max_points_for_view(args: VisualizationRequest, view: str) -> int:
    if args.max_points is not None:
        return int(args.max_points)
    return 50_000 if view == "3d" else 100_000


def _landscape_coordinates(landscape: Landscape, coordinate_source: str) -> np.ndarray:
    column = f"coordinates_{coordinate_source}"
    if column not in landscape.data.columns:
        raise ValueError(f"visualize requires landscape coordinate field: {column}")
    values = np.stack(landscape.data[column].to_numpy())
    if values.ndim != 2 or values.shape[1] != 3 or not np.isfinite(values).all():
        raise ValueError(f"{column} must contain finite three-component coordinates")
    return np.asarray(values, dtype=float)


def _resolve_visualize_color_scale(
    values,
    *,
    display_outlier_mask: np.ndarray | None,
    color_vmin: float | None,
    color_vmax: float | None,
    display_threshold: float | None,
) -> ColorScale:
    """Resolve the shared run-like legacy scale without changing stored density."""

    raw = np.asarray(values, dtype=float)
    finite = raw[np.isfinite(raw)]
    if not finite.size:
        raise ValueError("visualize requires at least one finite sld_display value")
    vmin = color_vmin if color_vmin is not None else display_threshold
    vmax = color_vmax if color_vmax is not None else min(float(np.max(finite)), 100.0)
    vmax_source = "explicit" if color_vmax is not None else "full_candidate_cap_100"
    if vmin is not None and vmax <= vmin:
        if color_vmax is not None:
            raise ValueError("--vmin must be less than --vmax")
        vmax = float(vmin) + max(abs(float(vmin)) * 0.01, 0.1)
        vmax_source = "constant_range_expansion"
    return ColorScale(
        None if vmin is None else float(vmin),
        float(vmax),
        "explicit" if color_vmin is not None else ("display_threshold" if vmin is not None else "auto"),
        vmax_source,
    )


def _validate_axis_limit_mapping(
    limits: dict[str, tuple[float | None, float | None] | None],
) -> None:
    valid = {"alpha", "beta", "gamma", "x", "y", "z"}
    unknown = sorted(set(limits).difference(valid))
    if unknown:
        raise ValueError(f"Unknown --axis-limit axis: {', '.join(unknown)}")
    for axis, bound in limits.items():
        if bound is None or bound[0] is None or bound[1] is None:
            raise ValueError(f"--axis-limit {axis} requires both LOWER and UPPER")
        if not np.isfinite(bound).all() or bound[0] >= bound[1]:
            raise ValueError(f"--axis-limit {axis} requires finite LOWER < UPPER")


def _validate_range_mapping(
    ranges: dict[str, tuple[float | None, float | None] | None],
) -> None:
    valid = {"alpha", "beta", "gamma", "x", "y", "z"}
    unknown = sorted(set(ranges).difference(valid))
    if unknown:
        raise ValueError(f"Unknown --range axis: {', '.join(unknown)}")
    for axis, bound in ranges.items():
        if bound is None or (bound[0] is None and bound[1] is None):
            raise ValueError(f"--range {axis} requires at least one finite bound")
        finite_bounds = [value for value in bound if value is not None]
        if not np.isfinite(finite_bounds).all():
            raise ValueError(f"--range {axis} bounds must be finite")


def _prepare_visualization_output_dir(output_dir: Path, *, overwrite: bool) -> None:
    if output_dir.exists():
        if not overwrite:
            raise FileExistsError(f"Visualization output directory already exists: {output_dir}")
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)


def _view_count(stats: dict[str, Any], view: str) -> int:
    counts = stats.get(f"n_points_{view}") or {}
    return int(next(iter(counts.values()), 0)) if isinstance(counts, dict) else int(counts)


def _visualization_output_dir(args) -> Path:
    run_dir = Path(args.run_dir)
    visual_id = _safe_artifact_id(getattr(args, "visual_id", None) or "default", "visual-id")
    if getattr(args, "selection_id", None):
        if getattr(args, "use_selected_landscape", False):
            return run_dir / "visualizations" / "selections" / args.selection_id / "selected_landscape" / visual_id
        if args.space == "canonical":
            return (
                run_dir
                / "visualizations"
                / "selections"
                / args.selection_id
                / "parent_canonical"
                / args.canonical_id
                / visual_id
            )
        return run_dir / "visualizations" / "selections" / args.selection_id / "parent_raw" / visual_id
    if args.space == "canonical":
        return run_dir / "visualizations" / "canonical" / args.canonical_id / visual_id
    return run_dir / "visualizations" / "raw" / visual_id


def _safe_artifact_id(value: str, option: str) -> str:
    if not value or value in {".", ".."} or Path(value).name != value:
        raise ValueError(f"--{option} must be a single safe path component")
    return value


def _validate_visualize_color_field_args(args) -> None:
    if args.vmin is not None and not np.isfinite(args.vmin):
        raise ValueError("--vmin must be finite")
    if args.vmax is not None and not np.isfinite(args.vmax):
        raise ValueError("--vmax must be finite")
    if args.vmin is not None and args.vmax is not None and args.vmin >= args.vmax:
        raise ValueError("--vmin must be less than --vmax")


def _validate_visualize_color_field(landscape: Landscape, color_field: str) -> None:
    if color_field not in landscape.data.columns:
        raise ValueError(f"visualize requires landscape color field: {color_field}")
    values = pd.to_numeric(landscape.data[color_field], errors="coerce")
    if values.notna().sum() == 0:
        raise ValueError(f"visualize color field must be numeric or boolean-like: {color_field}")


def _parse_kde_bandwidth(value: str | float | None) -> str | float:
    if value is None:
        return "scott"
    if isinstance(value, (int, float)):
        parsed = float(value)
        if parsed <= 0 or not np.isfinite(parsed):
            raise ValueError("--kde-bandwidth float must be positive")
        return parsed
    normalized = str(value).strip().lower()
    if normalized in {"scott", "silverman"}:
        return normalized
    try:
        parsed = float(normalized)
    except ValueError as exc:
        raise ValueError("--kde-bandwidth must be scott, silverman, or a positive float") from exc
    if parsed <= 0 or not np.isfinite(parsed):
        raise ValueError("--kde-bandwidth float must be positive")
    return parsed


def _visualization_public_report(
    report: dict[str, Any],
    *,
    args: VisualizationRequest,
    source_landscape_path: Path,
    euler_metadata: dict[str, object],
    formats: tuple[str, ...],
    requested_views: tuple[str, ...],
    view_reports: dict[str, Any],
    generated_files: dict[str, str],
    preparation_stats: dict[str, Any],
    resolved_color_scale: ColorScale,
    axis_limits: dict[str, tuple[float | None, float | None] | None],
    selection_metadata: dict[str, object] | None,
) -> dict[str, Any]:
    report = dict(report)
    filter_mode = (
        "top_fraction"
        if args.top_fraction is not None
        else ("threshold" if args.threshold is not None else "all")
    )
    report.update(
        {
            "artifact_type": "visualization_report",
            "schema_version": "2",
            "status": "ok",
            "source_landscape_path": str(source_landscape_path),
            "space": args.space,
            "canonical_id": args.canonical_id if args.space == "canonical" else None,
            "resolved_by": dict(getattr(args, "resolved_by", None) or {}),
            "visual_id": args.visual_id,
            "coordinate_source_resolved": [_coordinate_source_for_landscape(args)],
            "coordinate_sources": [_coordinate_source_for_landscape(args)],
            "requested_views": list(requested_views),
            "view_reports": view_reports,
            "representation": args.representation,
            "representations": list(_representations(args.representation)),
            "inherited_euler_convention": euler_metadata["euler_convention"],
            "inherited_euler_convention_source": euler_metadata["euler_convention_source"],
            "display_density_field": "sld_display",
            "color_field": "sld_display",
            "color_field_usage": "display_only",
            "colormap": args.colormap,
            "color_map": args.colormap,
            "display_filter_mode": filter_mode,
            "display_top_fraction": args.top_fraction,
            "display_sld_threshold": args.threshold,
            "display_filter_threshold": args.threshold,
            "display_range_filter_applied": bool(args.range_bound),
            "range_bounds": dict(args.range_bound or ()),
            "top_fraction": args.top_fraction,
            "threshold": args.threshold,
            "formats": list(formats),
            "vmin": args.vmin,
            "vmax": args.vmax,
            "color_vmin": resolved_color_scale.vmin,
            "color_vmax": resolved_color_scale.vmax,
            "display_color_vmax": resolved_color_scale.vmax,
            "color_vmin_source": resolved_color_scale.vmin_source,
            "color_vmax_source": resolved_color_scale.vmax_source,
            "point_size": args.point_size,
            "point_alpha": args.alpha,
            "axis_limits": {
                key: [value[0], value[1]] if value is not None else None
                for key, value in axis_limits.items()
            },
            "distribution_1d_policy": {
                "enabled": "1d" in requested_views,
                "bins": args.bins,
                "hist_mode": args.hist_mode,
                "kde": args.kde,
                "kde_bandwidth": args.kde_bandwidth,
                "full_filtered_rows": True,
                "coordinate_density_only": True,
            },
            "three_d_policy": {
                "enabled": "3d" in requested_views,
                "mode": args.three_d_mode if "3d" in requested_views else None,
                "display_only": True,
                "selection_enabled": False,
            },
            "generated_files": generated_files,
        }
    )
    report.update(preparation_stats)
    if selection_metadata is not None:
        report.update(selection_metadata)
    return report


def _resolve_selected_landscape_path(run_dir: Path, selection_id: str) -> Path:
    selection_dir = run_dir / "selections" / selection_id / "selected_landscape"
    npz_path = selection_dir / "landscape.npz"
    csv_path = selection_dir / "landscape.csv"
    if npz_path.exists():
        return npz_path
    if csv_path.exists():
        return csv_path
    raise ValueError(
        "Selected-derived landscape not found; expected landscape.npz or "
        f"landscape.csv under {selection_dir}"
    )


def _selected_landscape_visualization_metadata(
    run_dir: Path,
    selection_id: str,
    selected_landscape_path: Path,
    *,
    selected_count: int,
    report: dict[str, object] | None = None,
) -> dict[str, object]:
    report_path = run_dir / "selections" / selection_id / "selected_landscape" / "landscape_report.json"
    report = report if report is not None else (_read_json_if_exists(report_path) or {})
    output_paths = report.get("output_paths") if isinstance(report, dict) else None
    return {
        "selection_id": selection_id,
        "selected_count": selected_count,
        "total_count": report.get("parent_landscape_row_count"),
        "selection_filter_applied": False,
        "selected_landscape_used": True,
        "selected_landscape_path": str(selected_landscape_path),
        "selected_landscape_report_path": str(report_path) if report_path.exists() else None,
        "selected_landscape_parent_path": report.get("parent_landscape_path"),
        "parent_landscape_path": report.get("parent_landscape_path"),
        "selected_landscape_coordinate_space": report.get("coordinate_space"),
        "coordinate_space": report.get("coordinate_space"),
        "selected_landscape_density_source": report.get("density_source"),
        "density_source": report.get("density_source"),
        "selected_landscape_sld_recomputed": report.get("sld_recomputed"),
        "sld_recomputed": report.get("sld_recomputed"),
        "selected_landscape_output_paths": output_paths if isinstance(output_paths, dict) else None,
    }


def _read_selected_landscape_report(run_dir: Path, selection_id: str) -> dict[str, object] | None:
    return _read_json_if_exists(
        run_dir / "selections" / selection_id / "selected_landscape" / "landscape_report.json"
    )


def _filter_landscape_by_selection(args, landscape: Landscape) -> tuple[Landscape, dict[str, object]]:
    if not getattr(args, "run_dir", None):
        raise ValueError("visualize --selection-id requires --run-dir")
    selection_keys, selection_metadata = _read_visualization_selection_keys(
        Path(args.run_dir),
        args.selection_id,
        landscape_count=len(landscape.data),
    )
    if "particle_key" not in landscape.data.columns:
        raise ValueError("Landscape must contain particle_key for --selection-id filtering")

    selected_key_strings = [str(key) for key in selection_keys]
    selected_key_set = set(selected_key_strings)
    landscape_key_strings = landscape.data["particle_key"].map(str)
    missing_keys = sorted(selected_key_set.difference(set(landscape_key_strings)))
    if missing_keys:
        preview = ", ".join(missing_keys[:5])
        suffix = "" if len(missing_keys) <= 5 else f", ... ({len(missing_keys)} total)"
        raise ValueError(
            "Selection contains particle_key values absent from the requested landscape: "
            f"{preview}{suffix}"
        )

    filtered_data = landscape.data.loc[landscape_key_strings.isin(selected_key_set)].copy(deep=True)
    if filtered_data.empty:
        raise ValueError("Selection filtering produced an empty landscape")
    filtered_landscape = Landscape(
        data=filtered_data.reset_index(drop=True),
        canonical_transform=landscape.canonical_transform,
        active_policies=landscape.active_policies,
        density_report=None,
        canonicalization_report=landscape.canonicalization_report,
    )
    selection_metadata.update(
        {
            "selection_filter_applied": True,
            "selected_landscape_used": False,
            "landscape_count_before_selection_filter": int(len(landscape.data)),
            "landscape_count_after_selection_filter": int(len(filtered_data)),
            "missing_selected_key_count": 0,
        }
    )
    return filtered_landscape, selection_metadata


def _read_visualization_selection_keys(
    run_dir: Path,
    selection_id: str,
    *,
    landscape_count: int,
) -> tuple[list[object], dict[str, object]]:
    selection_dir = run_dir / "selections" / selection_id
    if not selection_dir.exists():
        raise ValueError(f"Selection directory does not exist: {selection_dir}")
    if not selection_dir.is_dir():
        raise ValueError(f"Selection path is not a directory: {selection_dir}")

    selected_keys_csv = selection_dir / "selected_particle_keys.csv"
    selection_json = selection_dir / "selection.json"
    if selected_keys_csv.exists():
        keys = _read_selected_particle_keys_csv(selected_keys_csv)
        metadata: dict[str, object] = {
            "selection_id": selection_id,
            "selected_count": int(len(keys)),
            "total_count": int(landscape_count),
            "selection_source_path": str(selected_keys_csv),
            "selection_count_source": "selected_particle_keys_csv",
            "selection_total_count_source": "landscape_row_count",
        }
        if selection_json.exists():
            counts = _read_selection_json_counts(selection_json)
            if counts.get("selected_count") is not None:
                metadata["selected_count"] = int(counts["selected_count"])
                if len(keys) != counts["selected_count"]:
                    raise ValueError(
                        "selected_particle_keys.csv row count does not match "
                        "selection.json selected_count"
                    )
            if counts.get("total_count") is not None:
                metadata["total_count"] = int(counts["total_count"])
                metadata["selection_total_count_source"] = "selection_json"
            if counts.get("selection_id") not in (None, selection_id):
                raise ValueError(
                    "selection.json selection_id does not match requested --selection-id"
                )
        return keys, metadata

    if selection_json.exists():
        selection = read_selection_json(selection_json)
        return list(selection.selected_particle_keys), {
            "selection_id": selection_id,
            "selected_count": int(selection.selected_count),
            "total_count": int(selection.total_count),
            "selection_source_path": str(selection_json),
            "selection_filter_applied": True,
            "selection_count_source": "selection_json",
            "selection_total_count_source": "selection_json",
        }

    raise ValueError(
        "Selection artifacts not found; expected selected_particle_keys.csv "
        f"or selection.json under {selection_dir}"
    )


def _read_selected_particle_keys_csv(path: Path) -> list[str]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"Selected particle keys CSV is empty: {path}")
        if "particle_key" not in reader.fieldnames:
            raise ValueError("selected_particle_keys.csv must contain a particle_key column")
        keys = [row["particle_key"] for row in reader if row.get("particle_key") not in (None, "")]
    if not keys:
        raise ValueError(f"Selected particle keys CSV contains no particle keys: {path}")
    return keys


def _read_selection_json_counts(path: Path) -> dict[str, object]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Malformed selection JSON: {path}") from exc
    if not isinstance(payload, dict):
        raise ValueError("Selection JSON must be an object")

    counts: dict[str, object] = {"selection_id": payload.get("selection_id")}
    for field_name in ("selected_count", "total_count"):
        value = payload.get(field_name)
        if value is None:
            counts[field_name] = None
            continue
        if not isinstance(value, int) or value < 0:
            raise ValueError(f"Selection JSON {field_name} must be a non-negative integer")
        counts[field_name] = value
    selected_count = counts.get("selected_count")
    total_count = counts.get("total_count")
    if (
        isinstance(selected_count, int)
        and isinstance(total_count, int)
        and selected_count > total_count
    ):
        raise ValueError("Selection JSON selected_count must be <= total_count")
    return counts


def _coordinate_source_for_landscape(args) -> str:
    if args.run_dir and args.space == "raw":
        return "analysis"
    if args.run_dir and args.space == "canonical":
        return "canonical"
    return args.coordinate_source


def _effective_visualization_space(
    args,
    *,
    selected_landscape_report: dict[str, object] | None,
) -> str:
    if getattr(args, "use_selected_landscape", False) and selected_landscape_report:
        coordinate_space = selected_landscape_report.get("coordinate_space")
        if coordinate_space in {"raw", "canonical"}:
            return str(coordinate_space)
    return getattr(args, "space", "raw")


def _resolve_visualization_euler_metadata(
    args,
    *,
    selected_landscape_report: dict[str, object] | None,
) -> dict[str, object]:
    report_metadata = (
        selected_landscape_report.get("euler_metadata")
        if isinstance(selected_landscape_report, dict)
        else None
    )
    if isinstance(report_metadata, dict):
        convention = report_metadata.get("euler_convention")
        sequence = report_metadata.get("scipy_euler_sequence")
        if convention or sequence:
            resolved = resolve_euler_convention(
                convention if isinstance(convention, str) else None,
                scipy_euler_sequence=sequence if isinstance(sequence, str) else None,
                source="inherited",
            )
            return resolved.metadata(euler_angle_columns=RAW_EULER_ANGLE_COLUMNS)
    return _resolve_landscape_euler_metadata(
        args,
        columns=RAW_EULER_ANGLE_COLUMNS,
    )


def _resolve_landscape_euler_metadata(
    args,
    *,
    columns=RAW_EULER_ANGLE_COLUMNS,
    explicit_sequence: str | None = None,
    warn_on_legacy_missing: bool = True,
) -> dict[str, object]:
    if getattr(args, "euler_convention", None):
        resolved = resolve_euler_convention(args.euler_convention, source="cli_override")
        return resolved.metadata(euler_angle_columns=columns)
    if explicit_sequence:
        resolved = resolve_euler_convention(
            scipy_euler_sequence=explicit_sequence,
            source="cli_override",
        )
        return resolved.metadata(euler_angle_columns=columns)
    inherited = _read_parent_euler_metadata(args, columns=columns)
    if inherited is not None:
        return inherited
    resolved = resolve_euler_convention(
        DEFAULT_EULER_CONVENTION,
        source=LEGACY_MISSING_EULER_CONVENTION_SOURCE,
    )
    metadata = resolved.metadata(euler_angle_columns=columns)
    if warn_on_legacy_missing:
        warn_user(
            "parent landscape lacks Euler convention metadata; "
            f"defaulting to {DEFAULT_EULER_CONVENTION} for derived Euler coordinates."
        )
    return metadata


def _read_parent_euler_metadata(args, *, columns) -> dict[str, object] | None:
    run_dir_raw = getattr(args, "run_dir", None)
    if not run_dir_raw:
        return None
    run_dir = Path(run_dir_raw)
    payloads: list[dict[str, object]] = []
    if getattr(args, "space", "raw") == "canonical":
        canonical_id = getattr(args, "canonical_id", "default")
        payload = _read_json_if_exists(
            run_dir / "canonical" / canonical_id / "canonicalize_summary.json"
        )
        if payload is not None:
            payloads.append(payload)
    run_summary = _read_json_if_exists(run_dir / "run_summary.json")
    if run_summary is not None:
        payloads.append(run_summary)
    for payload in payloads:
        convention = payload.get("euler_convention")
        sequence = payload.get("scipy_euler_sequence")
        if convention or sequence:
            resolved = resolve_euler_convention(
                convention if isinstance(convention, str) else None,
                scipy_euler_sequence=sequence if isinstance(sequence, str) else None,
                source="inherited",
            )
            return resolved.metadata(euler_angle_columns=columns)
    return None


def _read_json_if_exists(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        return None
    return payload


def _parse_formats(value: tuple[str, ...] | list[str] | str | None) -> tuple[str, ...]:
    raw_values = ["png"] if value is None else ([value] if isinstance(value, str) else list(value))
    formats = tuple(
        dict.fromkeys(
            part.strip().lower().lstrip(".")
            for raw in raw_values
            for part in str(raw).split(",")
            if part.strip()
        )
    )
    if not formats:
        raise ValueError("--format must include at least one format")
    unsupported = sorted(set(formats).difference({"png", "svg", "pdf"}))
    if unsupported:
        raise ValueError(f"Unsupported figure format: {', '.join(unsupported)}")
    return formats
