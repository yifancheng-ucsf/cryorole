"""Canonicalization workflow service independent of CLI parsing."""

from __future__ import annotations

from dataclasses import dataclass, fields
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from cryorole.canonicalize import canonicalize_landscape, canonicalize_landscape_arrays
from cryorole.canonicalize.transforms import apply_canonical_transform, validate_canonical_transform
from cryorole.core.display_policy import max_displayed_density
from cryorole.core.euler_conventions import CANONICAL_EULER_ANGLE_COLUMNS, DEFAULT_EULER_CONVENTION, resolve_euler_convention
from cryorole.export import (
    landscape_from_arrays,
    read_landscape,
    read_landscape_metadata,
    read_landscape_npz_arrays,
    write_canonical_landscape_csv,
    write_canonical_landscape_csv_from_npz,
    write_json_artifact,
    write_landscape_npz,
    write_landscape_npz_arrays,
    write_landscape_visualizations,
    write_report_json,
)
from cryorole.export.serialization import to_json_safe
from cryorole.io.landscape_resolver import resolve_landscape_path
from cryorole.models.landscape import Landscape
from cryorole.models.landscape_arrays import LandscapeArrays
from cryorole.models.policies import CanonicalizationPolicy
from cryorole.run_bundle import validate_completed_run_bundle
from cryorole.workflows.progress import current_rss_bytes

CANONICAL_FRAME_TRANSFORM_DIRECTION = "canonical_rv = raw_rv @ canonical_transform"
CANONICAL_FRAME_COORDINATE_SPACE = "rotvec_ro_radians"


@dataclass(frozen=True)
class CanonicalizeRequest:
    """Typed canonicalization request independent of argparse."""

    run_dir: str | None = None
    landscape: str | None = None
    output_dir: str | None = None
    canonical_id: str = "default"
    fit_top_fraction: float | None = 0.40
    positive_side: str | None = "low"
    legacy_positive_side: str | None = None
    use_frame: str | None = None
    no_visualize: bool = False
    no_csv: bool = False
    csv_chunk_size: int = 100_000
    profile_memory: bool = False
    overwrite: bool = False
    # How an omitted --run-dir was filled in (cryorole.workflow.resolve).
    resolved_by: dict[str, dict[str, str]] | None = None

    @classmethod
    def from_namespace(cls, namespace: Any) -> "CanonicalizeRequest":
        values = vars(namespace)
        return cls(
            **{
                field.name: values[field.name]
                for field in fields(cls)
                if field.name in values
            }
        )


@dataclass(frozen=True)
class RunCanonicalizationResult:
    """Canonical artifacts produced inside a transactional run staging tree."""

    landscape: Landscape
    landscape_npz_path: Path
    landscape_csv_path: Path
    report_path: Path


def canonicalize_run_artifacts(
    *,
    output_dir: Path,
    run_backend: str,
    policy: CanonicalizationPolicy,
    overwrite: bool,
    csv_chunk_size: int,
    euler_sequence: str,
    raw_arrays: LandscapeArrays | None,
    compatibility_result: object | None,
    runner: object,
    density_report: object,
    density_policy: object,
    preview_max_points: int,
) -> RunCanonicalizationResult:
    """Write hidden ``run --canonicalize`` artifacts through this service."""

    output_dir.mkdir(parents=True, exist_ok=True)
    if run_backend == "array_native":
        if raw_arrays is None:
            raise ValueError("array-native run canonicalization requires raw arrays")
        canonical_arrays, report = canonicalize_landscape_arrays(
            raw_arrays,
            policy=policy,
        )
        npz_path = write_landscape_npz_arrays(
            canonical_arrays,
            output_dir / "canonical_landscape.npz",
            overwrite=overwrite,
            artifact_type="canonical_landscape",
        )
        csv_path = write_canonical_landscape_csv_from_npz(
            npz_path,
            output_dir / "canonical_landscape.csv",
            overwrite=overwrite,
            chunk_size=csv_chunk_size,
            euler_sequence=euler_sequence,
        )
        if canonical_arrays.n_points <= preview_max_points:
            indices = np.arange(canonical_arrays.n_points, dtype=np.int64)
        else:
            indices = np.unique(
                np.linspace(
                    0,
                    canonical_arrays.n_points - 1,
                    num=preview_max_points,
                    dtype=np.int64,
                )
            )
        landscape = landscape_from_arrays(canonical_arrays, indices)
        landscape.active_policies = {"density_policy": density_policy}
        landscape.density_report = (
            density_report if len(indices) == canonical_arrays.n_points else None
        )
        landscape.canonicalization_report = report
    else:
        canonical_result = runner.canonicalize_density_result(
            compatibility_result,
            canonicalization_policy=policy,
        )
        landscape = canonical_result.landscape
        npz_path = write_landscape_npz(
            landscape,
            output_dir / "canonical_landscape.npz",
            overwrite=overwrite,
            artifact_type="canonical_landscape",
        )
        csv_path = write_canonical_landscape_csv(
            landscape,
            output_dir / "canonical_landscape.csv",
            overwrite=overwrite,
            euler_sequence=euler_sequence,
        )
    report_path = write_report_json(
        landscape.canonicalization_report,
        output_dir / "canonicalization_report.json",
        overwrite=overwrite,
        artifact_type="canonicalization_report",
    )
    return RunCanonicalizationResult(
        landscape=landscape,
        landscape_npz_path=npz_path,
        landscape_csv_path=csv_path,
        report_path=report_path,
    )


@dataclass(frozen=True)
class CanonicalizeResult:
    """Where ``canonicalize`` wrote the canonical frame and what it contains."""

    output_dir: Path
    summary: dict[str, Any]


def canonicalize(request: CanonicalizeRequest) -> int:
    """Compatibility wrapper returning an exit status; see ``canonicalize_bundle``."""

    canonicalize_bundle(request)
    return 0


def canonicalize_bundle(request: CanonicalizeRequest) -> CanonicalizeResult:
    """Derive (or apply) a canonical frame and write it; never prints."""

    args = request
    memory_profiler = _MemoryProfiler(enabled=args.profile_memory)
    memory_profiler.sample("start")
    if not args.run_dir and not (args.landscape and args.output_dir):
        raise ValueError("canonicalize requires --run-dir")
    if args.run_dir:
        validate_completed_run_bundle(args.run_dir)
    landscape_source = _landscape_source_from_args(args)
    _validate_canonicalize_frame_conflicts(args)
    if args.csv_chunk_size <= 0:
        raise ValueError("--csv-chunk-size must be positive")
    fit_top_fraction = _canonicalize_fit_top_fraction(args)
    if not 0.0 < fit_top_fraction <= 1.0:
        raise ValueError("--fit-top-fraction must be > 0 and <= 1")
    positive_side = _canonicalize_positive_side(args)
    output_dir = _prepare_output_dir(
        _canonicalize_output_dir(args),
        overwrite=args.overwrite,
        label="Canonicalize",
    )
    memory_profiler.sample("output_dir_prepared")
    source_path = resolve_landscape_path(landscape_source)
    source_metadata = read_landscape_metadata(source_path)
    euler_metadata = _canonicalize_euler_metadata()
    memory_profiler.sample("source_metadata_read")
    frame_source_metadata = _read_canonical_frame(args.use_frame) if args.use_frame else None
    policy = _canonicalization_policy_for_command(
        fit_top_fraction=fit_top_fraction,
        positive_side=positive_side,
        frame_metadata=frame_source_metadata,
    )
    frame_applied = frame_source_metadata is not None

    if source_path.suffix.lower() == ".npz":
        canonicalization_backend = "array_native"
        arrays = read_landscape_npz_arrays(source_path)
        memory_profiler.sample("array_landscape_loaded")
        if frame_source_metadata is None:
            canonical_arrays, canonicalization_report = canonicalize_landscape_arrays(
                arrays,
                policy=policy,
            )
        else:
            canonical_arrays = _apply_canonical_frame_arrays(
                arrays,
                frame_source_metadata["canonical_transform"],
            )
            canonicalization_report = _frame_applied_canonicalization_report(
                policy=policy,
                row_count=canonical_arrays.n_points,
                frame_metadata=frame_source_metadata,
            )
        memory_profiler.sample("array_canonicalization_complete")
        output_row_count = canonical_arrays.n_points
        canonical_landscape_path = write_landscape_npz_arrays(
            canonical_arrays,
            output_dir / "canonical_landscape.npz",
            overwrite=args.overwrite,
            artifact_type="canonical_landscape",
        )
        memory_profiler.sample("canonical_npz_written")
        canonical_landscape_csv_path = None
        csv_backend = "skipped"
        if not args.no_csv:
            canonical_landscape_csv_path = write_canonical_landscape_csv_from_npz(
                canonical_landscape_path,
                output_dir / "canonical_landscape.csv",
                overwrite=args.overwrite,
                chunk_size=args.csv_chunk_size,
                euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
            )
            csv_backend = "array_native_chunked"
            memory_profiler.sample("canonical_csv_written")
        visualization_landscape = _landscape_from_arrays_for_visualization(
            canonical_arrays,
            canonicalization_report=canonicalization_report,
        )
        canonical_transform = canonical_arrays.canonical_transform
    else:
        canonicalization_backend = "dataframe_compat"
        landscape = read_landscape(source_path)
        memory_profiler.sample("dataframe_landscape_loaded")
        if frame_source_metadata is None:
            canonical_landscape = canonicalize_landscape(landscape, policy=policy)
        else:
            canonical_landscape = _apply_canonical_frame_landscape(
                landscape,
                frame_source_metadata["canonical_transform"],
                policy=policy,
                frame_metadata=frame_source_metadata,
            )
        memory_profiler.sample("dataframe_canonicalization_complete")
        canonicalization_report = canonical_landscape.canonicalization_report
        output_row_count = len(canonical_landscape.data)
        canonical_landscape_path = write_landscape_npz(
            canonical_landscape,
            output_dir / "canonical_landscape.npz",
            overwrite=args.overwrite,
            artifact_type="canonical_landscape",
        )
        memory_profiler.sample("canonical_npz_written")
        canonical_landscape_csv_path = None
        csv_backend = "skipped"
        if not args.no_csv:
            canonical_landscape_csv_path = write_canonical_landscape_csv(
                canonical_landscape,
                output_dir / "canonical_landscape.csv",
                overwrite=args.overwrite,
                euler_sequence=str(euler_metadata["scipy_euler_sequence"]),
            )
            csv_backend = "dataframe_compat"
            memory_profiler.sample("canonical_csv_written")
        visualization_landscape = canonical_landscape
        canonical_transform = canonical_landscape.canonical_transform

    if canonical_transform is None:
        raise ValueError("canonicalization did not produce canonical_transform")
    frame_paths = _write_canonical_frame_artifacts(
        output_dir=output_dir,
        canonical_transform=canonical_transform,
        policy=policy,
        args=args,
        euler_metadata=euler_metadata,
        source_frame_metadata=frame_source_metadata,
        overwrite=args.overwrite,
    )
    canonicalization_report_payload = _canonicalization_report_payload(
        canonicalization_report,
        frame_paths=frame_paths,
        frame_applied=frame_applied,
        frame_source_metadata=frame_source_metadata,
    )
    canonicalization_report_payload["parent_run_id"] = (
        _run_id_from_bundle(args.run_dir) if args.run_dir else None
    )
    canonicalization_report_path = write_report_json(
        canonicalization_report_payload,
        output_dir / "canonicalization_report.json",
        overwrite=args.overwrite,
        artifact_type="canonicalization_report",
    )
    memory_profiler.sample("canonicalization_report_written")
    output_artifacts = {
        "canonical_landscape_npz": str(canonical_landscape_path),
        "canonical_landscape_csv": (
            str(canonical_landscape_csv_path)
            if canonical_landscape_csv_path is not None
            else None
        ),
        "canonicalization_report_json": str(canonicalization_report_path),
        "canonical_frame_json": frame_paths["canonical_frame_json"],
        "canonical_frame_npz": frame_paths["canonical_frame_npz"],
    }
    memory_profile_path = output_dir / "canonical_memory_profile.json"
    if args.profile_memory:
        output_artifacts["memory_profile_json"] = str(memory_profile_path)
    visualization_performed = False
    visualization_report_path = None
    if args.no_visualize:
        pass
    else:
        visualization_report = _write_canonical_default_visualization(
            args=args,
            canonical_landscape=visualization_landscape,
            canonical_landscape_csv_path=canonical_landscape_csv_path,
            euler_metadata=euler_metadata,
            fit_top_fraction=policy.fit_top_fraction,
            frame_source_metadata=frame_source_metadata,
            overwrite=args.overwrite,
        )
        visualization_performed = True
        visualization_report_path = visualization_report["report_path"]
        output_artifacts["canonical_visualization_report_json"] = str(visualization_report_path)

    canonicalize_summary = {
            "artifact_type": "canonicalize_summary",
            "schema_version": "1",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "source_landscape_path": source_metadata["path"],
            "source_run_dir": args.run_dir,
            "parent_run_id": _run_id_from_bundle(args.run_dir) if args.run_dir else None,
            "resolved_by": dict(args.resolved_by or {}),
            "source_landscape_artifact_type": source_metadata["artifact_type"],
            "source_landscape_schema_version": source_metadata["schema_version"],
            "source_landscape_row_count": source_metadata["row_count"],
            "canonicalization_policy": to_json_safe(policy),
            "fit_top_fraction": policy.fit_top_fraction,
            "canonicalization_performed": not frame_applied,
            "frame_applied": frame_applied,
            "frame_source_path": (
                frame_source_metadata["frame_source_path"] if frame_source_metadata else None
            ),
            "frame_metadata": (
                frame_source_metadata["metadata"] if frame_source_metadata else None
            ),
            "canonicalization_backend": canonicalization_backend,
            "canonical_landscape_npz": str(canonical_landscape_path),
            "csv_performed": canonical_landscape_csv_path is not None,
            "csv_backend": csv_backend,
            "csv_chunk_size": args.csv_chunk_size,
            "canonical_landscape_csv": (
                str(canonical_landscape_csv_path)
                if canonical_landscape_csv_path is not None
                else None
            ),
            "canonical_frame_json": frame_paths["canonical_frame_json"],
            "canonical_frame_npz": frame_paths["canonical_frame_npz"],
            "visualization_performed": visualization_performed,
            "visualization_report": visualization_report_path,
            "memory_profile_performed": args.profile_memory,
            "memory_profile": str(memory_profile_path) if args.profile_memory else None,
            "selection_performed": False,
            "updated_star_cs_export_performed": False,
            "composite_map_export_performed": False,
            "future_optional_exports": (
                "updated STAR/CS metadata and composite map export are planned "
                "future derived outputs and were not performed."
            ),
            "output_row_count": output_row_count,
            "output_artifacts": output_artifacts,
        }
    canonicalize_summary.update(euler_metadata)
    write_json_artifact(
        canonicalize_summary,
        output_dir / "canonicalize_summary.json",
        overwrite=args.overwrite,
    )
    memory_profiler.sample("canonicalize_summary_written")
    if args.profile_memory:
        memory_profiler.write(memory_profile_path, overwrite=args.overwrite)
    return CanonicalizeResult(output_dir=output_dir, summary=canonicalize_summary)


def _canonicalize_fit_top_fraction(args) -> float:
    return 0.40 if args.fit_top_fraction is None else float(args.fit_top_fraction)


def _canonicalize_positive_side(args) -> str:
    if getattr(args, "legacy_positive_side", None):
        return str(args.legacy_positive_side)
    if args.positive_side == "high":
        return "high_density_skew"
    return "low_density_skew"


def _validate_canonicalize_frame_conflicts(args) -> None:
    if not getattr(args, "use_frame", None):
        return
    conflicts = []
    if args.fit_top_fraction is not None:
        conflicts.append("--fit-top/--fit-top-fraction")
    if args.positive_side is not None:
        conflicts.append("--positive-side")
    if getattr(args, "legacy_positive_side", None) is not None:
        conflicts.append("--skewness-positive-side")
    if conflicts:
        raise ValueError("--use-frame conflicts with explicit fitting controls: " + ", ".join(conflicts))


def _canonicalize_euler_metadata() -> dict[str, object]:
    return resolve_euler_convention(
        DEFAULT_EULER_CONVENTION,
        source="public_default",
    ).metadata(euler_angle_columns=CANONICAL_EULER_ANGLE_COLUMNS)


def _canonicalization_policy_for_command(
    *,
    fit_top_fraction: float,
    positive_side: str,
    frame_metadata: dict[str, Any] | None,
) -> CanonicalizationPolicy:
    if frame_metadata is not None:
        metadata = frame_metadata["metadata"]
        fit_top_fraction = float(metadata.get("fit_top_fraction", fit_top_fraction))
        positive_side = str(metadata.get("positive_side", positive_side))
    return CanonicalizationPolicy(
        fit_top_fraction=fit_top_fraction,
        sign_rule="density_weighted_skewness",
        positive_side=positive_side,
    )


def _read_canonical_frame(path: str | None) -> dict[str, Any] | None:
    if not path:
        return None
    frame_path = Path(path)
    if not frame_path.exists():
        raise ValueError(f"Canonical frame does not exist: {frame_path}")
    if frame_path.suffix.lower() != ".json":
        raise ValueError("canonicalize --use-frame expects canonical_frame.json")
    with frame_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    if not isinstance(metadata, dict):
        raise ValueError("canonical_frame.json must contain a JSON object")
    transform = _canonical_transform_from_frame_metadata(metadata, frame_path)
    transform_direction = metadata.get("transform_direction")
    if transform_direction != CANONICAL_FRAME_TRANSFORM_DIRECTION:
        raise ValueError(
            "canonical_frame transform_direction must be "
            f"{CANONICAL_FRAME_TRANSFORM_DIRECTION!r}"
        )
    coordinate_space = metadata.get("coordinate_space")
    if coordinate_space != CANONICAL_FRAME_COORDINATE_SPACE:
        raise ValueError(
            "canonical_frame coordinate_space must be "
            f"{CANONICAL_FRAME_COORDINATE_SPACE!r}"
        )
    if metadata.get("euler_convention") not in {None, DEFAULT_EULER_CONVENTION}:
        raise ValueError(f"canonical_frame euler_convention must be {DEFAULT_EULER_CONVENTION!r}")
    return {
        "frame_source_path": str(frame_path),
        "canonical_transform": transform,
        "metadata": metadata,
    }


def _canonical_transform_from_frame_metadata(metadata: dict[str, Any], frame_path: Path) -> np.ndarray:
    if "canonical_transform" in metadata:
        transform = np.asarray(metadata["canonical_transform"], dtype=float)
    else:
        npz_value = metadata.get("canonical_frame_npz")
        if not npz_value:
            raise ValueError("canonical_frame.json requires canonical_transform or canonical_frame_npz")
        npz_path = Path(str(npz_value))
        if not npz_path.is_absolute():
            npz_path = frame_path.parent / npz_path
        with np.load(npz_path, allow_pickle=False) as payload:
            if "canonical_transform" not in payload:
                raise ValueError("canonical_frame.npz missing canonical_transform")
            transform = np.asarray(payload["canonical_transform"], dtype=float)
    validate_canonical_transform(transform)
    return transform


def _apply_canonical_frame_arrays(
    arrays: LandscapeArrays,
    canonical_transform: np.ndarray,
) -> LandscapeArrays:
    canonical_coordinates = apply_canonical_transform(
        arrays.coordinates_analysis,
        canonical_transform,
    )
    return LandscapeArrays(
        particle_key=arrays.particle_key,
        coordinates_analysis=arrays.coordinates_analysis,
        coordinates_display=arrays.coordinates_display,
        sld_unfloored=arrays.sld_unfloored,
        sld_raw=arrays.sld_raw,
        sld_display=arrays.sld_display,
        sld_display_is_outlier=arrays.sld_display_is_outlier,
        sld_was_floored=arrays.sld_was_floored,
        sld_local_k_mean=arrays.sld_local_k_mean,
        sld_effective_local_k_mean=arrays.sld_effective_local_k_mean,
        sld_distance_floor=arrays.sld_distance_floor,
        ref_source_row_id=arrays.ref_source_row_id,
        mov_source_row_id=arrays.mov_source_row_id,
        coordinates_canonical=canonical_coordinates,
        canonical_transform=canonical_transform,
    )


def _apply_canonical_frame_landscape(
    landscape: Landscape,
    canonical_transform: np.ndarray,
    *,
    policy: CanonicalizationPolicy,
    frame_metadata: dict[str, Any],
) -> Landscape:
    data = landscape.data.copy(deep=True)
    coordinates = np.vstack([np.asarray(value, dtype=float) for value in data["coordinates_analysis"]])
    canonical_coordinates = apply_canonical_transform(coordinates, canonical_transform)
    data["coordinates_canonical"] = [row.copy() for row in canonical_coordinates]
    active_policies = dict(landscape.active_policies or {})
    active_policies["canonicalization_policy"] = policy
    return Landscape(
        data=data,
        canonical_transform=canonical_transform,
        active_policies=active_policies,
        density_report=landscape.density_report,
        canonicalization_report=_frame_applied_canonicalization_report(
            policy=policy,
            row_count=len(data),
            frame_metadata=frame_metadata,
        ),
    )


def _frame_applied_canonicalization_report(
    *,
    policy: CanonicalizationPolicy,
    row_count: int,
    frame_metadata: dict[str, Any],
) -> dict[str, Any]:
    source_metadata = frame_metadata["metadata"]
    return {
        "fit_subset": "imported_frame",
        "n_fit_points": 0,
        "density_support_field": source_metadata.get("density_support_field", "sld_raw"),
        "singular_values": (),
        "explained_variance_ratios": (),
        "axis_degeneracy_detected": False,
        "fit_top_fraction": policy.fit_top_fraction,
        "axis_assignment": source_metadata.get("axis_assignment", policy.axis_assignment),
        "pca_axis_order": (),
        "assigned_coordinate_names": ("alpha", "beta", "gamma"),
        "assigned_rotvec_columns": ("z", "y", "x"),
        "sign_rule": source_metadata.get("sign_rule", policy.sign_rule),
        "positive_side": source_metadata.get("positive_side", policy.positive_side),
        "sign_weight_field": policy.sign_weight_field,
        "axis_weighted_skewness": (),
        "flipped_axes": (),
        "ambiguous_sign_axes": (),
        "handedness_rule": policy.handedness_rule,
        "handedness_adjusted_axes": (),
        "transform_direction": CANONICAL_FRAME_TRANSFORM_DIRECTION,
        "warnings": (),
        "frame_applied": True,
        "frame_source_path": frame_metadata["frame_source_path"],
        "input_row_count": int(row_count),
    }


def _write_canonical_frame_artifacts(
    *,
    output_dir: Path,
    canonical_transform: np.ndarray,
    policy: CanonicalizationPolicy,
    args,
    euler_metadata: dict[str, object],
    source_frame_metadata: dict[str, Any] | None,
    overwrite: bool,
) -> dict[str, str]:
    frame_npz_path = output_dir / "canonical_frame.npz"
    if frame_npz_path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {frame_npz_path}")
    np.savez(
        frame_npz_path,
        schema_version=np.asarray("1"),
        canonical_transform=np.asarray(canonical_transform, dtype=float),
    )
    frame_json_path = output_dir / "canonical_frame.json"
    metadata = {
        "artifact_type": "canonical_frame",
        "schema_version": "1",
        "canonical_transform": np.asarray(canonical_transform, dtype=float),
        "canonical_frame_npz": str(frame_npz_path),
        "transform_direction": CANONICAL_FRAME_TRANSFORM_DIRECTION,
        "coordinate_space": CANONICAL_FRAME_COORDINATE_SPACE,
        "axis_assignment": policy.axis_assignment,
        "sign_rule": policy.sign_rule,
        "positive_side": policy.positive_side,
        "fit_top_fraction": policy.fit_top_fraction,
        "density_support_field": policy.density_support_field,
        "euler_convention": euler_metadata["euler_convention"],
        "source_run_dir": args.run_dir,
        "parent_run_id": _run_id_from_bundle(args.run_dir) if args.run_dir else None,
        "source_canonical_id": args.canonical_id,
        "frame_applied": source_frame_metadata is not None,
        "frame_source_path": (
            source_frame_metadata["frame_source_path"] if source_frame_metadata else None
        ),
    }
    if source_frame_metadata is not None:
        metadata["source_frame_metadata"] = source_frame_metadata["metadata"]
    write_json_artifact(metadata, frame_json_path, overwrite=overwrite)
    return {
        "canonical_frame_json": str(frame_json_path),
        "canonical_frame_npz": str(frame_npz_path),
    }


def _canonicalization_report_payload(
    report: Any,
    *,
    frame_paths: dict[str, str],
    frame_applied: bool,
    frame_source_metadata: dict[str, Any] | None,
) -> dict[str, Any]:
    payload = dict(to_json_safe(report))
    payload.update(
        {
            "frame_applied": frame_applied,
            "frame_source_path": (
                frame_source_metadata["frame_source_path"] if frame_source_metadata else None
            ),
            "frame_metadata": (
                frame_source_metadata["metadata"] if frame_source_metadata else None
            ),
            "canonical_frame_json": frame_paths["canonical_frame_json"],
            "canonical_frame_npz": frame_paths["canonical_frame_npz"],
        }
    )
    return payload


def _landscape_from_arrays_for_visualization(
    arrays: LandscapeArrays,
    *,
    canonicalization_report: Any,
) -> Landscape:
    data = pd.DataFrame(
        {
            "particle_key": arrays.particle_key.astype(str),
            "coordinates_analysis": list(arrays.coordinates_analysis),
            "coordinates_display": list(
                arrays.coordinates_display
                if arrays.coordinates_display is not None
                else arrays.coordinates_analysis
            ),
            "coordinates_canonical": list(arrays.coordinates_canonical),
            "sld_unfloored": arrays.sld_unfloored,
            "sld_raw": arrays.sld_raw,
            "sld_display": arrays.sld_display,
            "sld_display_is_outlier": arrays.sld_display_is_outlier,
            "sld_was_floored": arrays.sld_was_floored,
            "sld_local_k_mean": arrays.sld_local_k_mean,
            "sld_effective_local_k_mean": arrays.sld_effective_local_k_mean,
            "sld_distance_floor": arrays.sld_distance_floor,
        }
    )
    if arrays.ref_source_row_id is not None:
        data["ref_source_row_id"] = arrays.ref_source_row_id
    if arrays.mov_source_row_id is not None:
        data["mov_source_row_id"] = arrays.mov_source_row_id
    return Landscape(
        data=data,
        canonical_transform=arrays.canonical_transform,
        canonicalization_report=canonicalization_report,
    )


def _write_canonical_default_visualization(
    *,
    args,
    canonical_landscape: Landscape,
    canonical_landscape_csv_path: Path | None,
    euler_metadata: dict[str, object],
    fit_top_fraction: float,
    frame_source_metadata: dict[str, Any] | None,
    overwrite: bool,
) -> dict[str, Any]:
    output_dir = _canonical_visualization_output_dir(args)
    top_group = f"filter_particles_by_top_sld_{_fit_fraction_pct_label(fit_top_fraction)}pct"
    warnings: list[str] = []
    if frame_source_metadata is not None and "fit_top_fraction" not in frame_source_metadata["metadata"]:
        warnings.append(
            "canonical frame metadata missing fit_top_fraction; "
            "using current canonicalize default/policy for top-SLD preview"
        )

    preview_specs = (
        {
            "label": "all_particles",
            "description": "all canonical landscape rows",
            "display_filter_mode": "none",
            "display_sld_threshold": None,
            "display_top_fraction": None,
        },
        {
            "label": "filter_particles_by_sld_gt_1p5",
            "description": "display-only preview of particles with sld_raw > 1.5",
            "display_filter_mode": "threshold",
            "display_sld_threshold": 1.5,
            "display_top_fraction": None,
        },
        {
            "label": top_group,
            "description": "display-only preview of the top-sld_raw canonical fit support",
            "display_filter_mode": "top_fraction",
            "display_sld_threshold": None,
            "display_top_fraction": fit_top_fraction,
        },
    )

    preview_reports: dict[str, dict[str, Any]] = {}
    generated_files: dict[str, str] = {}
    full_landscape_table = (
        str(canonical_landscape_csv_path) if canonical_landscape_csv_path is not None else None
    )
    for spec in preview_specs:
        label = str(spec["label"])
        max_displayed_sld = max_displayed_density(
            canonical_landscape.data["sld_raw"],
            threshold=spec["display_sld_threshold"],
            strict_threshold=True,
        )
        resolved_vmax = (
            min(max_displayed_sld, 100.0) if max_displayed_sld is not None else None
        )
        color_vmax = resolved_vmax
        if color_vmax is None and spec["display_sld_threshold"] is not None:
            color_vmax = float(spec["display_sld_threshold"])
        selection_metadata: dict[str, Any] = {
            "preview_group": label,
            "preview_description": spec["description"],
            "preview_only": True,
            "not_a_selection": True,
            "display_filter_mode": spec["display_filter_mode"],
            "display_vmax_cap": 100.0,
            "display_vmax_policy": "min(max_displayed_sld, 100)",
            "max_displayed_sld": max_displayed_sld,
            "resolved_vmax": resolved_vmax,
        }
        if spec["display_sld_threshold"] is not None:
            selection_metadata["display_sld_threshold"] = spec["display_sld_threshold"]
            selection_metadata["display_filter_label"] = "sld_raw > 1.5"
        if spec["display_top_fraction"] is not None:
            selection_metadata["top_sld_fraction"] = spec["display_top_fraction"]
            selection_metadata["fit_support_preview"] = True
        report = write_landscape_visualizations(
            canonical_landscape,
            output_dir / label,
            overwrite=overwrite,
            coordinate_source="canonical",
            representation="both",
            color_field="sld_raw",
            euler_convention=DEFAULT_EULER_CONVENTION,
            euler_convention_source=str(euler_metadata["euler_convention_source"]),
            display_top_fraction=spec["display_top_fraction"],
            display_sld_threshold=spec["display_sld_threshold"],
            display_density_field="sld_raw",
            formats=("png",),
            max_points_2d=500000,
            max_points_3d=50000,
            random_seed=0,
            write_projection_csvs=False,
            full_landscape_table=full_landscape_table,
            display_table_filename="display_table.csv",
            artifact_layout="run_bundle",
            visual_style="legacy",
            color_map="rainbow_r",
            color_vmax=color_vmax,
            display_filter_mode=str(spec["display_filter_mode"]),
            selection_metadata=selection_metadata,
        )
        preview_reports[label] = report
        generated_files[f"{label}_visualization_report_json"] = report["report_path"]
        for artifact_key, artifact_path in report["generated_files"].items():
            generated_files[f"{label}_{artifact_key}"] = artifact_path

    aggregate_report = {
        "artifact_type": "canonical_default_visualization_report",
        "schema_version": "1",
        "status": "ok",
        "canonical_id": args.canonical_id,
        "preview_only": True,
        "not_a_selection": True,
        "visual_style": "legacy_rainbow",
        "color_map": "rainbow_r",
        "color_field": "sld_raw",
        "display_vmax_cap": 100.0,
        "display_vmax_policy": "min(max_displayed_sld, 100)",
        "fit_top_fraction": fit_top_fraction,
        "top_sld_preview_group": top_group,
        "preview_groups": [str(spec["label"]) for spec in preview_specs],
        "preview_reports": {
            label: report["report_path"] for label, report in preview_reports.items()
        },
        "generated_files": generated_files,
        "warnings": warnings,
    }
    report_path = write_json_artifact(
        aggregate_report,
        output_dir / "visualization_report.json",
        overwrite=overwrite,
    )
    aggregate_report["report_path"] = str(report_path)
    return aggregate_report


def _canonical_visualization_output_dir(args) -> Path:
    if args.run_dir:
        return Path(args.run_dir) / "visualizations" / "canonical" / args.canonical_id
    return Path(args.output_dir) / "visualizations"


def _fit_fraction_pct_label(fraction: float) -> str:
    percent = float(fraction) * 100.0
    rounded = round(percent)
    if np.isclose(percent, rounded):
        return str(int(rounded))
    return f"{percent:.3g}".replace(".", "p")


def _landscape_source_from_args(args):
    if getattr(args, "run_dir", None):
        return args.run_dir
    if getattr(args, "landscape", None):
        return args.landscape
    raise ValueError("Provide --run-dir or --landscape")


def _canonicalize_output_dir(args) -> Path:
    if args.output_dir:
        return Path(args.output_dir)
    return Path(args.run_dir) / "canonical" / args.canonical_id


def _read_json_if_exists(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        return None
    return payload


def _run_id_from_bundle(run_dir: str | Path) -> str | None:
    for filename in ("run_summary.json", "run_manifest.json"):
        payload = _read_json_if_exists(Path(run_dir) / filename)
        if payload is not None and payload.get("run_id") is not None:
            return str(payload["run_id"])
    return None


def _prepare_output_dir(
    output_dir: str | Path,
    *,
    overwrite: bool,
    label: str,
) -> Path:
    path = Path(output_dir)
    if path.exists() and not path.is_dir():
        raise ValueError(f"{label} output path exists and is not a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(f"{label} output directory already exists: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


class _MemoryProfiler:
    """Small RSS sampler for optional CLI diagnostics."""

    def __init__(self, *, enabled: bool) -> None:
        self.enabled = enabled
        self.samples: list[dict[str, Any]] = []
        self.backend: str | None = None

    def sample(self, stage: str) -> None:
        if not self.enabled:
            return
        rss_bytes, backend = current_rss_bytes()
        if self.backend is None:
            self.backend = backend
        start_rss = self.samples[0]["rss_bytes"] if self.samples else rss_bytes
        self.samples.append(
            {
                "stage": stage,
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "rss_bytes": rss_bytes,
                "rss_mib": rss_bytes / (1024 * 1024),
                "delta_from_start_bytes": rss_bytes - start_rss,
                "delta_from_start_mib": (rss_bytes - start_rss) / (1024 * 1024),
            }
        )

    def write(self, path: Path, *, overwrite: bool) -> Path:
        peak_rss = max((sample["rss_bytes"] for sample in self.samples), default=None)
        payload = {
            "artifact_type": "canonical_memory_profile",
            "schema_version": "1",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "metric": "rss",
            "units": "bytes",
            "backend": self.backend,
            "sample_count": len(self.samples),
            "peak_rss_bytes": peak_rss,
            "peak_rss_mib": (
                peak_rss / (1024 * 1024) if peak_rss is not None else None
            ),
            "samples": self.samples,
        }
        return write_json_artifact(payload, path, overwrite=overwrite)
