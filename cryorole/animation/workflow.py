"""Phase 1-4 orchestration for ``cryorole animate``."""

from __future__ import annotations

import json
import shutil
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping

import numpy as np

from cryorole.animation.artifacts import (
    AnimationLandscape,
    load_animation_landscape,
    load_canonical_frame,
    resolve_recorded_euler_convention,
)
from cryorole.animation.chimerax import (
    ChimeraXExecutionError,
    execute_chimerax,
    generate_chimerax_scripts,
    normalize_model_id_groups,
    validate_chimerax_completion,
    validate_matching_session_baselines,
    validate_structure_frames,
    write_scene_transforms_csv,
)
from cryorole.animation.compositor import (
    CompositeResult,
    compose_frame_sequences,
    resolve_composite_layout,
    validate_dual_structure_horizontal_crop,
    validate_structure_crop_fraction,
)
from cryorole.animation.manifest import (
    ANIMATION_MANIFEST_SCHEMA_VERSION,
    write_animation_manifest,
)
from cryorole.animation.projection import ProjectionRenderer, resolve_axis_limits
from cryorole.animation.trajectory import (
    build_trajectory,
    parse_waypoint_csv,
    reverse_waypoints,
    write_trajectory_csv,
)
from cryorole.animation.transforms import (
    canonical_rotation_to_raw,
    conjugate_active_delta,
    load_explicit_scene_basis,
    parse_ro_assertion,
    pivoted_scene_delta,
    ro_target_to_active_delta,
)
from cryorole.animation.video import (
    FFmpegExecutionError,
    FFmpegOutputError,
    FFmpegResult,
    FFprobeExecutionError,
    MovieValidation,
    MovieValidationError,
    encode_mp4,
    validate_mp4,
)
from cryorole.models.policies import AnimationPolicy


@dataclass(frozen=True)
class AnimationResult:
    output_dir: Path
    manifest_path: Path
    status: str
    frame_count: int


def run_animation(args: SimpleNamespace | Any) -> AnimationResult:
    """Build a Phase 1-4 animation bundle without touching source artifacts."""

    run_dir = Path(args.run_dir).resolve()
    path_csv = Path(args.path_csv).resolve()
    session = Path(args.chimerax_session).resolve()
    secondary_session = (
        Path(args.secondary_chimerax_session).resolve()
        if getattr(args, "secondary_chimerax_session", None)
        else None
    )
    tertiary_session = (
        Path(args.tertiary_chimerax_session).resolve()
        if getattr(args, "tertiary_chimerax_session", None)
        else None
    )
    reference_model_ids, moving_model_ids = normalize_model_id_groups(
        args.reference_model_id,
        args.moving_model_id,
    )
    if not run_dir.is_dir():
        raise ValueError(f"Run directory does not exist: {run_dir}")
    if not path_csv.is_file():
        raise ValueError(f"Waypoint CSV does not exist: {path_csv}")
    if not session.is_file():
        raise ValueError(f"ChimeraX session does not exist: {session}")
    if secondary_session is not None and not secondary_session.is_file():
        raise ValueError(f"Secondary ChimeraX session does not exist: {secondary_session}")
    if tertiary_session is not None and secondary_session is None:
        raise ValueError(
            "--tertiary-chimerax-session requires --secondary-chimerax-session"
        )
    if tertiary_session is not None and not tertiary_session.is_file():
        raise ValueError(f"Tertiary ChimeraX session does not exist: {tertiary_session}")
    if secondary_session is not None and args.layout != "stacked":
        raise ValueError("multiple structure views require --layout stacked")
    view_sessions = {"primary": session}
    if secondary_session is not None:
        view_sessions["secondary"] = secondary_session
    if tertiary_session is not None:
        view_sessions["tertiary"] = tertiary_session
    multiple_views = len(view_sessions) > 1
    structure_horizontal_crop = validate_dual_structure_horizontal_crop(
        getattr(
            args,
            "structure_horizontal_crop",
            getattr(args, "dual_structure_horizontal_crop", 0.0),
        ),
        has_secondary=secondary_session is not None,
        layout_name=args.layout,
    )
    structure_vertical_crop = validate_structure_crop_fraction(
        getattr(args, "structure_vertical_crop", 0.0),
        axis_name="vertical",
        has_multiple_views=multiple_views,
        layout_name=args.layout,
    )
    structure_view_count = len(view_sessions)
    movie_name = _validated_movie_name(args.movie_name)
    output_dir = _prepare_output_directory(
        Path(args.output_dir),
        run_dir=run_dir,
        protected_paths=tuple(
            path
            for path in (path_csv, *view_sessions.values())
            if path is not None
        ),
        overwrite=bool(args.overwrite),
    )
    log_path = output_dir / "logs" / "cryorole_animate.log"
    manifest_path = output_dir / "manifest.json"
    warnings: list[str] = []
    return_code = None
    ffmpeg_return_code = None
    ffprobe_return_code = None
    structure_validation = None
    completion_status = None
    structure_validations: dict[str, object] = {}
    completion_statuses: dict[str, dict[str, object]] = {}
    return_codes: dict[str, int] = {}
    session_transform_parity = {
        "status": (
            (
                "not_checked_script_only"
                if args.render_mode == "script-only"
                else "pending"
            )
            if multiple_views
            else "not_applicable_single_view"
        ),
        "view_count": structure_view_count,
    }
    composite_result = None
    ffmpeg_result = None
    movie_validation = None
    active_stage = "trajectory"
    stage_statuses = {
        "trajectory": "not_started",
        "landscape": "not_started",
        "chimerax_scripts": "not_started",
        "structure_render": "not_started",
        "compositor": "not_started",
        "encoder": "not_started",
        "movie_validation": "not_started",
    }
    if multiple_views:
        stage_statuses.update(
            {
                f"structure_render_{view_name}": "not_started"
                for view_name in view_sessions
            }
        )
    manifest = _base_manifest(
        args,
        run_dir,
        path_csv,
        session,
        secondary_session=secondary_session,
        tertiary_session=tertiary_session,
        reference_model_ids=reference_model_ids,
        moving_model_ids=moving_model_ids,
    )
    try:
        if args.render_mode == "execute" and not args.no_encode:
            active_stage = "encoder"
            _validate_encoder_preflight(args)
            active_stage = "trajectory"
        _log(log_path, "resolve", "Resolving landscape and Euler provenance")
        euler = resolve_recorded_euler_convention(
            run_dir,
            coordinate_set=args.coordinate_set,
            canonical_id=args.canonical_id,
            requested=args.euler_convention,
        )
        canonical_frame = (
            load_canonical_frame(run_dir, args.canonical_id or "default")
            if args.coordinate_set == "canonical" or args.map_frame == "canonical"
            else None
        )
        waypoints, parse_warnings = parse_waypoint_csv(
            path_csv,
            path_space=args.path_space,
            scipy_euler_sequence=euler.scipy_euler_sequence,
            default_segment_frames=args.frames_per_segment,
            default_hold_frames=args.hold_frames,
        )
        warnings.extend(parse_warnings)
        if args.reverse:
            waypoints = reverse_waypoints(waypoints)
        frames = build_trajectory(waypoints, fps=args.fps, ping_pong=args.ping_pong)
        write_trajectory_csv(
            frames,
            output_dir / "trajectory.csv",
            scipy_euler_sequence=euler.scipy_euler_sequence,
        )
        raw_rotations = tuple(
            canonical_rotation_to_raw(frame.rotation, canonical_frame.transform)
            if args.coordinate_set == "canonical"
            else frame.rotation
            for frame in frames
        )
        baseline, baseline_metadata = parse_ro_assertion(
            args.baseline_ro,
            scipy_euler_sequence=euler.scipy_euler_sequence,
            first_waypoint_raw=raw_rotations[0],
        )
        raw_to_scene, map_provenance = _resolve_map_frame(args, canonical_frame)
        pivot = np.asarray(args.pivot, dtype=float)
        if pivot.shape != (3,) or not np.isfinite(pivot).all():
            raise ValueError("pivot must contain exactly three finite scene coordinates")
        scene_transforms = np.stack(
            [
                pivoted_scene_delta(
                    conjugate_active_delta(
                        ro_target_to_active_delta(rotation, baseline).as_matrix(),
                        raw_to_scene,
                    ),
                    pivot,
                )
                for rotation in raw_rotations
            ]
        )
        transform_csv = write_scene_transforms_csv(
            scene_transforms,
            output_dir / "chimerax" / "frame_transforms.csv",
        )
        stage_statuses["trajectory"] = "complete"
        _log(log_path, "trajectory", f"Generated {len(frames)} SO(3) frames")

        active_stage = "landscape"
        landscape = load_animation_landscape(
            run_dir,
            coordinate_set=args.coordinate_set,
            canonical_id=args.canonical_id,
            scipy_euler_sequence=euler.scipy_euler_sequence,
            sld_threshold=args.threshold,
            top_fraction=args.top_fraction,
            range_bounds=dict(args.range_bound or ()),
            color_vmin=args.vmin,
            color_vmax=args.vmax,
        )
        trajectory_euler = np.vstack(
            [
                frame.rotation.as_euler(euler.scipy_euler_sequence, degrees=True)
                for frame in frames
            ]
        )
        explicit_limits = _read_axis_limits(args.axis_limits)
        range_bounds = dict(args.range_bound or ())
        axis_limits, expansions = resolve_axis_limits(
            landscape.euler_degrees,
            trajectory_euler,
            explicit=explicit_limits,
            range_bounds=range_bounds,
        )
        if expansions:
            warning = "Axis limits expanded for trajectory panels: " + ", ".join(expansions)
            warnings.append(warning)
            _log(log_path, "landscape", warning)
        planned_layout = resolve_composite_layout(
            canvas_size=(args.composite_width, args.composite_height),
            landscape_dimensions=(args.canvas_width, args.canvas_height),
            structure_dimensions=(args.structure_width, args.structure_height),
            secondary_structure_dimensions=(
                (args.structure_width, args.structure_height)
                if secondary_session is not None
                else None
            ),
            tertiary_structure_dimensions=(
                (args.structure_width, args.structure_height)
                if tertiary_session is not None
                else None
            ),
            layout_name=args.layout,
            background_color=args.background_color,
            dual_structure_horizontal_crop=structure_horizontal_crop,
            structure_vertical_crop=structure_vertical_crop,
        )
        landscape_destination = planned_layout.landscape_destination
        composite_scale = min(
            landscape_destination[2] / args.canvas_width,
            landscape_destination[3] / args.canvas_height,
        )
        renderer = ProjectionRenderer(
            landscape_euler=landscape.euler_degrees,
            sld_display=landscape.sld_display,
            sld_display_is_outlier=landscape.sld_display_is_outlier,
            trajectory_frames=frames,
            scipy_euler_sequence=euler.scipy_euler_sequence,
            axis_limits=axis_limits,
            canvas_width=args.canvas_width,
            canvas_height=args.canvas_height,
            point_size=args.point_size,
            landscape_alpha=args.landscape_alpha,
            colormap=args.colormap,
            color_vmin=args.vmin,
            color_vmax=args.vmax,
            display_threshold=args.threshold,
            resolved_color_scale=landscape.color_scale,
            trail_frames=args.trail_frames,
            show_current_values=args.show_current_values,
            composite_scale=composite_scale,
        )
        landscape_paths = renderer.render(output_dir / "frames" / "landscape")
        if len(landscape_paths) != len(frames):
            raise RuntimeError("Landscape frame count does not match trajectory")
        _log(log_path, "landscape", f"Rendered {len(landscape_paths)} landscape frames")
        stage_statuses["landscape"] = "complete"

        active_stage = "chimerax_scripts"
        if not multiple_views:
            scripts_by_view = {
                "primary": generate_chimerax_scripts(
                    output_dir=output_dir,
                    session_path=session,
                    transform_csv=transform_csv,
                    reference_model_id=reference_model_ids,
                    moving_model_id=moving_model_ids,
                    structure_width=args.structure_width,
                    structure_height=args.structure_height,
                )
            }
        else:
            scripts_by_view = {
                view_name: generate_chimerax_scripts(
                    output_dir=output_dir,
                    session_path=view_session,
                    transform_csv=transform_csv,
                    reference_model_id=reference_model_ids,
                    moving_model_id=moving_model_ids,
                    structure_width=args.structure_width,
                    structure_height=args.structure_height,
                    view_name=view_name,
                )
                for view_name, view_session in view_sessions.items()
            }
        stage_statuses["chimerax_scripts"] = "complete"
        status = "scripts_ready"
        if args.render_mode == "execute":
            if not args.chimerax_bin:
                raise ValueError("--chimerax-bin is required for --render-mode execute")
            for view_name, scripts in scripts_by_view.items():
                active_stage = (
                    f"structure_render_{view_name}"
                    if multiple_views
                    else "structure_render"
                )
                render_log = output_dir / "logs" / (
                    f"chimerax_render_{view_name}.log"
                    if multiple_views
                    else "chimerax_render.log"
                )
                _log(log_path, "chimerax", f"Starting {view_name} ChimeraX subprocess")
                result = execute_chimerax(
                    chimerax_bin=args.chimerax_bin,
                    cxc_script=scripts.cxc_script,
                    log_path=render_log,
                )
                return_codes[view_name] = result.returncode
                if view_name == "primary":
                    return_code = result.returncode
                completion_statuses[view_name] = validate_chimerax_completion(
                    scripts.status_path,
                    expected_count=len(frames),
                    expected_reference_model_ids=reference_model_ids,
                    expected_moving_model_ids=moving_model_ids,
                )
                structure_dir = (
                    output_dir / "frames" / "structure" / view_name
                    if multiple_views
                    else output_dir / "frames" / "structure"
                )
                structure_validations[view_name] = validate_structure_frames(
                    structure_dir,
                    expected_count=len(frames),
                )
                if active_stage in stage_statuses:
                    stage_statuses[active_stage] = "complete"
            return_code = return_codes["primary"]
            completion_status = completion_statuses["primary"]
            structure_validation = structure_validations["primary"]
            if multiple_views:
                try:
                    session_transform_parity = validate_matching_session_baselines(
                        completion_statuses["primary"],
                        completion_statuses["secondary"],
                        completion_statuses.get("tertiary"),
                    )
                except RuntimeError as exc:
                    session_transform_parity = {
                        "status": "failed",
                        "matched": False,
                        "error": str(exc),
                    }
                    raise
            status = "structure_rendered"
            stage_statuses["structure_render"] = "complete"
            _log(log_path, "chimerax", "All structure views validated")
            active_stage = "compositor"
            composite_result = compose_frame_sequences(
                landscape_dir=output_dir / "frames" / "landscape",
                structure_dir=(
                    output_dir / "frames" / "structure" / "primary"
                    if multiple_views
                    else output_dir / "frames" / "structure"
                ),
                secondary_structure_dir=(
                    output_dir / "frames" / "structure" / "secondary"
                    if secondary_session is not None
                    else None
                ),
                tertiary_structure_dir=(
                    output_dir / "frames" / "structure" / "tertiary"
                    if tertiary_session is not None
                    else None
                ),
                output_dir=output_dir / "frames" / "composite",
                expected_count=len(frames),
                canvas_size=(args.composite_width, args.composite_height),
                layout_name=args.layout,
                background_color=args.background_color,
                dual_structure_horizontal_crop=structure_horizontal_crop,
                structure_vertical_crop=structure_vertical_crop,
            )
            stage_statuses["compositor"] = "complete"
            status = "composite_rendered"
            _log(
                log_path,
                "compositor",
                f"Rendered {composite_result.frame_count} composite frames",
            )
            if args.no_encode:
                stage_statuses["encoder"] = "not_requested"
                stage_statuses["movie_validation"] = "not_requested"
            else:
                active_stage = "encoder"
                ffmpeg_result = encode_mp4(
                    ffmpeg_bin=args.ffmpeg_bin,
                    composite_dir=output_dir / "frames" / "composite",
                    output_path=output_dir / movie_name,
                    fps=args.fps,
                    crf=args.crf,
                    log_path=output_dir / "logs" / "ffmpeg_encode.log",
                )
                ffmpeg_return_code = ffmpeg_result.return_code
                stage_statuses["encoder"] = "complete"
                active_stage = "movie_validation"
                movie_validation = validate_mp4(
                    ffprobe_bin=args.ffprobe_bin,
                    movie_path=ffmpeg_result.output_path,
                    expected_dimensions=composite_result.dimensions,
                    expected_fps=args.fps,
                    expected_frame_count=len(frames),
                    log_path=output_dir / "logs" / "ffprobe_validate.log",
                )
                ffprobe_return_code = movie_validation.return_code
                stage_statuses["movie_validation"] = "complete"
                status = "movie_encoded"
                _log(log_path, "encoder", f"Validated MP4: {ffmpeg_result.output_path}")
        else:
            stage_statuses["structure_render"] = "not_requested"
            if multiple_views:
                for view_name in view_sessions:
                    stage_statuses[f"structure_render_{view_name}"] = "not_requested"
            stage_statuses["compositor"] = "not_requested"
            stage_statuses["encoder"] = "not_requested"
            stage_statuses["movie_validation"] = "not_requested"
        manifest.update(
            _success_manifest(
                args=args,
                status=status,
                euler=euler,
                canonical_frame=canonical_frame,
                baseline_metadata=baseline_metadata,
                raw_to_scene=raw_to_scene,
                map_provenance=map_provenance,
                pivot=pivot,
                frames=frames,
                landscape=landscape,
                axis_limits=axis_limits,
                expansions=expansions,
                stage_statuses=stage_statuses,
                return_code=return_code,
                return_codes=return_codes,
                output_dir=output_dir,
                structure_validation=structure_validation,
                structure_validations=structure_validations,
                completion_status=completion_status,
                completion_statuses=completion_statuses,
                session_transform_parity=session_transform_parity,
                primary_session=session,
                secondary_session=secondary_session,
                tertiary_session=tertiary_session,
                reference_model_ids=reference_model_ids,
                moving_model_ids=moving_model_ids,
                composite_result=composite_result,
                ffmpeg_result=ffmpeg_result,
                movie_validation=movie_validation,
                ffmpeg_return_code=ffmpeg_return_code,
                ffprobe_return_code=ffprobe_return_code,
                movie_name=movie_name,
                renderer=renderer,
                warnings=warnings,
            )
        )
        write_animation_manifest(manifest, manifest_path)
        _log(log_path, "complete", f"Animation bundle status: {status}")
        return AnimationResult(output_dir, manifest_path, status, len(frames))
    except Exception as exc:
        if isinstance(exc, ChimeraXExecutionError):
            return_code = exc.return_code
        if isinstance(exc, (FFmpegExecutionError, FFmpegOutputError)):
            ffmpeg_return_code = exc.return_code
        if isinstance(exc, (FFprobeExecutionError, MovieValidationError)):
            ffprobe_return_code = exc.return_code
        if active_stage in stage_statuses:
            stage_statuses[active_stage] = "failed"
        if active_stage.startswith("structure_render"):
            stage_statuses["structure_render"] = "failed"
        manifest.update(
            {
                "status": "failed",
                "structure_view_count": structure_view_count,
                "stage_statuses": stage_statuses,
                "ChimeraX_return_code": return_code,
                "ChimeraX_return_codes": return_codes,
                "ChimeraX_completions": completion_statuses,
                "structure_validations": {
                    name: {
                        "frame_count": validation.frame_count,
                        "dimensions": validation.dimensions,
                    }
                    for name, validation in structure_validations.items()
                },
                "session_transform_parity": session_transform_parity,
                "FFmpeg_return_code": ffmpeg_return_code,
                "FFprobe_return_code": ffprobe_return_code,
                "FFmpeg_arguments": getattr(exc, "arguments", None)
                if active_stage == "encoder"
                else None,
                "FFprobe_arguments": getattr(exc, "arguments", None)
                if active_stage == "movie_validation"
                else None,
                "compositor": _compositor_manifest(composite_result),
                "output_paths": {
                    "manifest": str(manifest_path),
                    "trajectory": str(output_dir / "trajectory.csv"),
                    "landscape_frames": str(output_dir / "frames" / "landscape"),
                    "structure_frames": (
                        str(output_dir / "frames" / "structure")
                        if (output_dir / "frames" / "structure").is_dir()
                        else None
                    ),
                    "composite_frames": (
                        str(output_dir / "frames" / "composite")
                        if (output_dir / "frames" / "composite").is_dir()
                        else None
                    ),
                    "movie": (
                        str(output_dir / movie_name)
                        if (output_dir / movie_name).is_file()
                        else None
                    ),
                    "ffmpeg_log": (
                        str(output_dir / "logs" / "ffmpeg_encode.log")
                        if (output_dir / "logs" / "ffmpeg_encode.log").is_file()
                        else None
                    ),
                    "ffprobe_log": (
                        str(output_dir / "logs" / "ffprobe_validate.log")
                        if (output_dir / "logs" / "ffprobe_validate.log").is_file()
                        else None
                    ),
                },
                "warnings": warnings,
                "errors": [f"{type(exc).__name__}: {exc}"],
            }
        )
        write_animation_manifest(manifest, manifest_path)
        _log(log_path, "failed", f"{type(exc).__name__}: {exc}")
        raise


def _resolve_map_frame(args, canonical_frame):
    if args.map_frame == "raw":
        if args.map_frame_transform:
            raise ValueError("--map-frame-transform is only valid with --map-frame explicit")
        return np.eye(3), {"kind": "raw", "source": "identity"}
    if args.map_frame == "canonical":
        if args.map_frame_transform:
            raise ValueError("--map-frame-transform is only valid with --map-frame explicit")
        if canonical_frame is None:
            raise ValueError("canonical map frame requires a validated canonical frame artifact")
        if not canonical_frame.physical_change_of_basis:
            raise ValueError(
                "canonical_frame is a validated RV coordinate frame but is not "
                "explicitly audited as a physical map change-of-basis; use "
                "--map-frame explicit with an audited raw_to_scene transform"
            )
        return canonical_frame.transform.T, {
            "kind": "canonical",
            "source": str(canonical_frame.source_path),
            "column_vector_convention": canonical_frame.raw_to_canonical_column_convention,
        }
    if args.map_frame == "explicit":
        if not args.map_frame_transform:
            raise ValueError("--map-frame explicit requires --map-frame-transform")
        return load_explicit_scene_basis(args.map_frame_transform), {
            "kind": "explicit",
            "source": str(Path(args.map_frame_transform).resolve()),
            "column_vector_convention": "x_scene = raw_to_scene @ x_raw",
        }
    raise ValueError("map_frame must be raw, canonical, or explicit")


def _prepare_output_directory(
    path: Path,
    *,
    run_dir: Path,
    protected_paths: tuple[Path, ...],
    overwrite: bool,
) -> Path:
    output = path.resolve()
    if output == run_dir or output in run_dir.parents:
        raise ValueError("Animation output directory must not be the run directory or its parent")
    for source in protected_paths:
        resolved_source = source.resolve()
        if output == resolved_source or output in resolved_source.parents:
            raise ValueError(
                "Animation output directory must not replace or contain a source artifact"
            )
    if output.exists():
        if not overwrite:
            raise FileExistsError(f"Animation output directory already exists: {output}")
        if output == Path(output.anchor) or len(output.parts) < 3:
            raise ValueError(f"Refusing unsafe animation output replacement: {output}")
        shutil.rmtree(output)
    output.mkdir(parents=True)
    return output


def _read_axis_limits(value: str | None) -> Mapping[str, list[float]] | None:
    if value is None:
        return None
    path = Path(value)
    if not path.is_file():
        raise ValueError(f"Axis-limit JSON does not exist: {path}")
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("Axis-limit JSON must contain an object")
    return payload


def _validated_movie_name(value: str) -> str:
    name = str(value).strip()
    path = Path(name)
    if not name or path.name != name or path.suffix.lower() != ".mp4":
        raise ValueError("--movie-name must be a basename ending in .mp4")
    return name


def _validate_encoder_preflight(args) -> None:
    if not args.ffmpeg_bin:
        raise ValueError("--ffmpeg-bin is required unless --no-encode is used")
    if not args.ffprobe_bin:
        raise ValueError("--ffprobe-bin is required unless --no-encode is used")
    if not Path(args.ffmpeg_bin).is_file():
        raise ValueError(f"FFmpeg executable does not exist: {args.ffmpeg_bin}")
    if not Path(args.ffprobe_bin).is_file():
        raise ValueError(f"FFprobe executable does not exist: {args.ffprobe_bin}")
    if not isinstance(args.crf, int) or not 0 <= args.crf <= 51:
        raise ValueError("--crf must be an integer in [0, 51]")


def _base_manifest(
    args,
    run_dir: Path,
    path_csv: Path,
    session: Path,
    *,
    secondary_session: Path | None,
    tertiary_session: Path | None,
    reference_model_ids: tuple[str, ...],
    moving_model_ids: tuple[str, ...],
) -> dict[str, object]:
    try:
        package_version = version("cryorole")
    except PackageNotFoundError:
        package_version = "unknown"
    resolved_args = {
        key: str(value) if isinstance(value, Path) else value
        for key, value in vars(args).items()
        if key != "handler"
    }
    return {
        "schema_version": ANIMATION_MANIFEST_SCHEMA_VERSION,
        "status": "building",
        "command": "cryorole animate",
        "cryorole_version": package_version,
        "git_commit": None,
        "creation_timestamp": datetime.now(timezone.utc).isoformat(),
        "run_dir": str(run_dir),
        "path_csv": str(path_csv),
        "ChimeraX_session": str(session),
        "reference_model_ids": list(reference_model_ids),
        "moving_model_ids": list(moving_model_ids),
        "reference_model_count": len(reference_model_ids),
        "moving_model_count": len(moving_model_ids),
        "structure_view_count": (
            1 + int(secondary_session is not None) + int(tertiary_session is not None)
        ),
        "dual_structure_horizontal_crop_fraction": (
            getattr(
                args,
                "structure_horizontal_crop",
                getattr(args, "dual_structure_horizontal_crop", 0.0),
            )
        ),
        "structure_horizontal_crop_fraction": getattr(
            args,
            "structure_horizontal_crop",
            getattr(args, "dual_structure_horizontal_crop", 0.0),
        ),
        "structure_vertical_crop_fraction": getattr(
            args, "structure_vertical_crop", 0.0
        ),
        "ChimeraX_sessions": {
            name: str(path)
            for name, path in (
                ("primary", session),
                ("secondary", secondary_session),
                ("tertiary", tertiary_session),
            )
            if path is not None
        },
        "resolved_CLI_arguments": resolved_args,
        "warnings": [],
        "errors": [],
    }


def _success_manifest(
    *,
    args,
    status: str,
    euler,
    canonical_frame,
    baseline_metadata,
    raw_to_scene,
    map_provenance,
    pivot,
    frames,
    landscape: AnimationLandscape,
    axis_limits,
    expansions,
    stage_statuses,
    return_code,
    return_codes,
    output_dir,
    structure_validation,
    structure_validations,
    completion_status,
    completion_statuses,
    session_transform_parity,
    primary_session,
    secondary_session,
    tertiary_session,
    reference_model_ids,
    moving_model_ids,
    composite_result: CompositeResult | None,
    ffmpeg_result: FFmpegResult | None,
    movie_validation: MovieValidation | None,
    ffmpeg_return_code,
    ffprobe_return_code,
    movie_name: str,
    renderer,
    warnings,
) -> dict[str, object]:
    animation_policy = AnimationPolicy(
        coordinate_set=args.coordinate_set,
        path_space=args.path_space,
        euler_convention=euler.euler_convention,
        frames_per_segment=args.frames_per_segment,
        fps=args.fps,
        hold_frames=args.hold_frames,
        reverse=args.reverse,
        ping_pong=args.ping_pong,
        map_frame=args.map_frame,
        render_mode=args.render_mode,
        composite_width=args.composite_width,
        composite_height=args.composite_height,
        composite_layout=args.layout,
        composite_background_color=args.background_color,
        encode_movie=args.render_mode == "execute" and not args.no_encode,
        crf=args.crf,
        movie_name=movie_name,
    )
    view_sessions = {"primary": primary_session}
    if secondary_session is not None:
        view_sessions["secondary"] = secondary_session
    if tertiary_session is not None:
        view_sessions["tertiary"] = tertiary_session
    multiple_views = len(view_sessions) > 1
    return {
        "status": status,
        "structure_view_count": len(view_sessions),
        "dual_structure_horizontal_crop_fraction": (
            getattr(
                args,
                "structure_horizontal_crop",
                getattr(args, "dual_structure_horizontal_crop", 0.0),
            )
        ),
        "structure_horizontal_crop_fraction": getattr(
            args,
            "structure_horizontal_crop",
            getattr(args, "dual_structure_horizontal_crop", 0.0),
        ),
        "structure_vertical_crop_fraction": getattr(
            args, "structure_vertical_crop", 0.0
        ),
        "ChimeraX_sessions": {
            view_name: str(view_session)
            for view_name, view_session in view_sessions.items()
        },
        "source_landscape": str(landscape.source_path),
        "canonical_frame": (
            {
                "path": str(canonical_frame.source_path),
                "canonical_id": canonical_frame.canonical_id,
                "transform": canonical_frame.transform,
                "transform_direction": canonical_frame.transform_direction,
                "physical_change_of_basis": canonical_frame.physical_change_of_basis,
            }
            if canonical_frame is not None
            else None
        ),
        "waypoint_coordinate_set": args.coordinate_set,
        "resolved_ro_space": "raw",
        "euler_convention": euler.euler_convention,
        "scipy_euler_sequence": euler.scipy_euler_sequence,
        "path_space": args.path_space,
        "quaternion_storage_order": "wxyz",
        "quaternion_internal_order": "xyzw",
        "baseline_ro_assertion": baseline_metadata,
        "map_session_frame": args.map_frame,
        "raw_to_scene_transform": raw_to_scene,
        "raw_to_scene_provenance": map_provenance,
        "matrix_vector_convention": "column vectors; x_scene = S @ x_raw",
        "pivot": pivot,
        "pivot_coordinate_frame": "ChimeraX scene",
        "reference_model_ID": (
            reference_model_ids[0] if len(reference_model_ids) == 1 else None
        ),
        "moving_model_ID": moving_model_ids[0] if len(moving_model_ids) == 1 else None,
        "reference_model_ids": list(reference_model_ids),
        "moving_model_ids": list(moving_model_ids),
        "reference_model_count": len(reference_model_ids),
        "moving_model_count": len(moving_model_ids),
        "model_group_policy": {
            "transform": "shared_absolute_scene_delta_from_each_model_baseline",
            "shared_trajectory": True,
            "shared_baseline_ro": True,
            "shared_pivot": True,
            "independent_model_motion": False,
        },
        "trajectory_policies": {
            "interpolation": "SO(3) quaternion SLERP",
            "frames_per_segment": args.frames_per_segment,
            "hold_frames": args.hold_frames,
            "reverse": args.reverse,
            "ping_pong": args.ping_pong,
            "segment_endpoints": "inclusive_shared_waypoint_deduplicated",
        },
        "animation_policy": asdict(animation_policy),
        "fps": args.fps,
        "frame_count": len(frames),
        "display_filter": {
            "threshold": args.threshold,
            "top_fraction": args.top_fraction,
            "range": dict(args.range_bound or ()),
            "display_only": True,
        },
        "resolved_visualize_display_policy": {
            "color_map": args.colormap,
            "point_size": args.point_size,
            "point_alpha": args.landscape_alpha,
            "color_vmin": args.vmin,
            "color_vmax": args.vmax,
            "resolved_color_vmin": renderer.resolved_vmin,
            "resolved_color_vmax": renderer.resolved_vmax,
            "color_vmin_source": renderer.color_vmin_source,
            "color_vmax_source": renderer.color_vmax_source,
            "aspect": "equal",
            "colorbar_position": "bottom",
            "sort_points_by_color": "ascending",
            "projection_order": ["alpha_beta", "beta_gamma", "alpha_gamma"],
            "projection_layout_version": renderer.projection_layout_version,
            "projection_panel_rectangles": renderer.panel_rectangles,
            "projection_rectangle_convention": "top_left_origin_xywh_pixels",
            "projection_gap_pixels": renderer.projection_gap_pixels,
            "projection_gap_fraction": renderer.projection_gap_fraction,
            "annotation_policy": renderer.annotation_policy,
            "shared_colorbar_policy": renderer.shared_colorbar_policy,
            "visible_colorbar_label": renderer.colorbar_label,
            "underlying_color_field": renderer.sld_field,
            "resolved_marker_diameter_pixels": (
                renderer.resolved_marker_diameter_pixels
            ),
            "resolved_marker_size_points2": renderer.marker_size_points2,
            "resolved_current_value_font_size_pixels": (
                renderer.resolved_current_value_font_size_pixels
            ),
            "resolved_current_value_font_size_points": (
                renderer.current_value_font_size_points
            ),
        },
        "resolved_SLD_field": "sld_display",
        "total_landscape_rows": landscape.total_rows,
        "displayed_landscape_rows": landscape.displayed_rows,
        "axis_limit_policy": (
            "explicit_json" if args.axis_limits else "range" if args.range_bound else "legacy_euler"
        ),
        "axis_limits": axis_limits,
        "axis_expansions": expansions,
        "canvas_size": [args.canvas_width, args.canvas_height],
        "structure_frame_size": [args.structure_width, args.structure_height],
        "composite_canvas_size": [args.composite_width, args.composite_height],
        "interpretation": (
            "Rigid-body rendering from the composite density; not an independently "
            "reconstructed density at each trajectory position."
        ),
        "stage_statuses": stage_statuses,
        "ChimeraX_return_code": return_code,
        "ChimeraX_completion": completion_status,
        "ChimeraX_return_codes": return_codes,
        "ChimeraX_completions": completion_statuses,
        "model_group_validation": _model_group_validation_manifest(
            completion_status,
            reference_model_ids=reference_model_ids,
            moving_model_ids=moving_model_ids,
        ),
        "session_transform_parity": session_transform_parity,
        "compositor": _compositor_manifest(composite_result),
        "encoding": _encoding_manifest(
            args=args,
            result=ffmpeg_result,
            return_code=ffmpeg_return_code,
        ),
        "movie_validation": _movie_validation_manifest(
            movie_validation,
            return_code=ffprobe_return_code,
        ),
        "structure_validation": (
            {
                "frame_count": structure_validation.frame_count,
                "dimensions": structure_validation.dimensions,
            }
            if structure_validation is not None
            else None
        ),
        "structure_validations": {
            name: {
                "frame_count": validation.frame_count,
                "dimensions": validation.dimensions,
            }
            for name, validation in structure_validations.items()
        },
        "output_paths": {
            "manifest": str(output_dir / "manifest.json"),
            "trajectory": str(output_dir / "trajectory.csv"),
            "landscape_frames": str(output_dir / "frames" / "landscape"),
            "structure_frames": (
                (
                    {
                        view_name: str(
                            output_dir / "frames" / "structure" / view_name
                        )
                        for view_name in view_sessions
                    }
                    if multiple_views
                    else str(output_dir / "frames" / "structure")
                )
                if args.render_mode == "execute"
                else None
            ),
            "composite_frames": (
                str(output_dir / "frames" / "composite")
                if composite_result is not None
                else None
            ),
            "movie": (
                str(output_dir / movie_name)
                if movie_validation is not None
                else None
            ),
            "chimerax_python": (
                {
                    view_name: str(
                        output_dir
                        / "chimerax"
                        / view_name
                        / "render_structure.py"
                    )
                    for view_name in view_sessions
                }
                if multiple_views
                else str(output_dir / "chimerax" / "render_structure.py")
            ),
            "chimerax_cxc": (
                {
                    view_name: str(
                        output_dir / "chimerax" / view_name / "run_render.cxc"
                    )
                    for view_name in view_sessions
                }
                if multiple_views
                else str(output_dir / "chimerax" / "run_render.cxc")
            ),
            "scene_transforms": str(output_dir / "chimerax" / "frame_transforms.csv"),
            "chimerax_render_status": (
                {
                    view_name: str(
                        output_dir
                        / "logs"
                        / f"chimerax_render_{view_name}_status.json"
                    )
                    for view_name in view_sessions
                }
                if multiple_views
                else str(output_dir / "logs" / "chimerax_render_status.json")
            ) if args.render_mode == "execute" else None,
            "chimerax_render_log": (
                {
                    view_name: str(
                        output_dir / "logs" / f"chimerax_render_{view_name}.log"
                    )
                    for view_name in view_sessions
                }
                if multiple_views
                else str(output_dir / "logs" / "chimerax_render.log")
            ) if args.render_mode == "execute" else None,
            "ffmpeg_log": (
                str(output_dir / "logs" / "ffmpeg_encode.log")
                if ffmpeg_result is not None
                else None
            ),
            "ffprobe_log": (
                str(output_dir / "logs" / "ffprobe_validate.log")
                if movie_validation is not None
                else None
            ),
        },
        "warnings": warnings,
        "errors": [],
    }


def _model_group_validation_manifest(
    completion: Mapping[str, object] | None,
    *,
    reference_model_ids,
    moving_model_ids,
) -> dict[str, object]:
    base = {
        "reference_model_ids": list(reference_model_ids),
        "moving_model_ids": list(moving_model_ids),
    }
    if completion is None:
        return {"status": "not_checked_script_only", **base}
    fields = (
        "reference_baseline_matrices",
        "moving_baseline_matrices",
        "reference_unchanged_by_model",
        "stationary_unchanged_by_model",
        "moving_restored_by_model",
        "moving_group_relative_transform_checks",
        "moving_group_relative_transforms_preserved",
    )
    return {
        "status": "validated",
        **base,
        **{field: completion.get(field) for field in fields},
    }


def _compositor_manifest(result: CompositeResult | None) -> dict[str, object] | None:
    if result is None:
        return None
    layout = result.layout
    return {
        "policy": (
            "fixed_horizontal_source_crop_then_contain"
            if layout.dual_structure_horizontal_crop_fraction > 0
            else "full_artifact_contain_letterbox"
        ),
        "layout": layout.name,
        "layout_version": layout.layout_version,
        "structure_view_count": layout.structure_view_count,
        "canvas_size": layout.canvas_size,
        "landscape_source_dimensions": result.landscape_validation.dimensions,
        "structure_source_dimensions": result.structure_validation.dimensions,
        "secondary_structure_source_dimensions": (
            result.secondary_structure_validation.dimensions
            if result.secondary_structure_validation is not None
            else None
        ),
        "tertiary_structure_source_dimensions": (
            result.tertiary_structure_validation.dimensions
            if result.tertiary_structure_validation is not None
            else None
        ),
        "dual_structure_horizontal_crop_fraction": (
            layout.dual_structure_horizontal_crop_fraction
        ),
        "structure_horizontal_crop_fraction": (
            layout.structure_horizontal_crop_fraction
        ),
        "structure_vertical_crop_fraction": layout.structure_vertical_crop_fraction,
        "structure_source_crop_rectangles": tuple(
            crop
            for crop in (
                layout.primary_structure_source_crop,
                layout.secondary_structure_source_crop,
                layout.tertiary_structure_source_crop,
            )
            if crop is not None
        ),
        "structure_cropped_dimensions": tuple(
            dimensions
            for dimensions in (
                layout.primary_structure_cropped_dimensions,
                layout.secondary_structure_cropped_dimensions,
                layout.tertiary_structure_cropped_dimensions,
            )
            if dimensions is not None
        ),
        "structure_alignments": tuple(
            alignment
            for alignment in (
                layout.primary_structure_alignment,
                layout.secondary_structure_alignment,
                layout.tertiary_structure_alignment,
            )
            if alignment is not None
        ),
        "landscape_region": layout.landscape_region,
        "structure_region": layout.structure_region,
        "landscape_destination": layout.landscape_destination,
        "structure_destination": layout.structure_destination,
        "structure_regions": tuple(
            region
            for region in (
                layout.primary_structure_region,
                layout.secondary_structure_region,
                layout.tertiary_structure_region,
            )
            if region is not None
        ),
        "structure_destinations": tuple(
            destination
            for destination in (
                layout.primary_structure_destination,
                layout.secondary_structure_destination,
                layout.tertiary_structure_destination,
            )
            if destination is not None
        ),
        "dual_view_gap_pixels": layout.dual_view_gap_pixels,
        "structure_view_gap_pixels": layout.structure_view_gap_pixels,
        "dual_view_gap_fraction": (
            layout.dual_view_gap_pixels / layout.canvas_size[0]
        ),
        "aspect_policy": layout.aspect_policy,
        "visible_colorbar_label": "SLD",
        "underlying_SLD_field": "sld_display",
        "background_color": layout.background_color,
        "resize_filter": layout.resize_filter,
        "color_mode": "RGB",
        "padding_pixels": layout.padding_pixels,
        "gap_pixels": layout.gap_pixels,
        "structure_fraction": layout.structure_fraction,
        "landscape_fraction": layout.landscape_fraction,
        "frame_count": result.frame_count,
        "dimensions": result.dimensions,
    }


def _encoding_manifest(
    *,
    args,
    result: FFmpegResult | None,
    return_code,
) -> dict[str, object]:
    return {
        "requested": args.render_mode == "execute" and not args.no_encode,
        "codec": "libx264",
        "pixel_format": "yuv420p",
        "fps": args.fps,
        "crf": args.crf,
        "faststart": True,
        "ffmpeg_bin": args.ffmpeg_bin,
        "arguments": result.arguments if result is not None else None,
        "return_code": return_code,
        "output_path": str(result.output_path) if result is not None else None,
        "log_path": str(result.log_path) if result is not None else None,
    }


def _movie_validation_manifest(
    result: MovieValidation | None,
    *,
    return_code,
) -> dict[str, object] | None:
    if result is None:
        return None
    return {
        "ffprobe_bin": result.arguments[0],
        "arguments": result.arguments,
        "return_code": return_code,
        "codec": result.codec,
        "pixel_format": result.pixel_format,
        "dimensions": result.dimensions,
        "fps": result.fps,
        "frame_count": result.frame_count,
        "duration_sec": result.duration_sec,
        "validation_method": result.validation_method,
        "frame_count_tolerance": result.frame_count_tolerance,
        "probe_payload": result.probe_payload,
        "log_path": str(result.log_path),
    }


def _log(path: Path, stage: str, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(
            json.dumps(
                {
                    "timestamp": datetime.now(timezone.utc).isoformat(),
                    "stage": stage,
                    "message": message,
                }
            )
            + "\n"
        )
