"""Camera-only ChimeraX views along validated canonical axes."""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from PIL import Image, UnidentifiedImageError

from cryorole.animation.artifacts import load_canonical_frame
from cryorole.animation.chimerax import ChimeraXExecutionError, execute_chimerax
from cryorole.animation.manifest import write_animation_manifest
from cryorole.animation.transforms import load_explicit_scene_basis, validate_scene_basis
from cryorole.canonicalize.transforms import validate_canonical_transform


CANONICAL_VIEWS_SCHEMA_VERSION = "1"


@dataclass(frozen=True)
class CanonicalViewBasis:
    """One right-handed screen basis expressed in ChimeraX scene coordinates."""

    name: str
    file_name: str
    outward: np.ndarray
    right: np.ndarray
    up: np.ndarray

    @property
    def camera_rotation(self) -> np.ndarray:
        return np.column_stack((self.right, self.up, self.outward))

    def manifest_record(self) -> dict[str, object]:
        return {
            "name": self.name,
            "file_name": self.file_name,
            "toward_viewer": self.outward,
            "screen_right": self.right,
            "screen_up": self.up,
            "camera_rotation": self.camera_rotation,
        }


@dataclass(frozen=True)
class CanonicalViewScripts:
    python_script: Path
    cxc_script: Path
    status_path: Path


@dataclass(frozen=True)
class CanonicalViewsResult:
    output_dir: Path
    manifest_path: Path
    status: str


def resolve_canonical_view_bases(
    canonical_transform: np.ndarray,
    raw_to_scene: np.ndarray,
) -> tuple[CanonicalViewBasis, ...]:
    """Resolve deterministic +X/+Y/+Z camera bases from ``axes_scene = S @ C``."""

    canonical = np.asarray(canonical_transform, dtype=float)
    validate_canonical_transform(canonical)
    scene = validate_scene_basis(raw_to_scene)
    axes = scene @ canonical
    x_axis, y_axis, z_axis = (axes[:, index] for index in range(3))
    views = (
        CanonicalViewBasis("canonical_x_plus", "canonical_x_plus.png", x_axis, y_axis, z_axis),
        CanonicalViewBasis("canonical_y_plus", "canonical_y_plus.png", y_axis, x_axis, -z_axis),
        CanonicalViewBasis("canonical_z_plus", "canonical_z_plus.png", z_axis, x_axis, y_axis),
    )
    for view in views:
        rotation = view.camera_rotation
        if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-8):
            raise ValueError(f"{view.name} camera basis is not orthonormal")
        if not np.isclose(np.linalg.det(rotation), 1.0, atol=1e-8):
            raise ValueError(f"{view.name} camera basis is not right-handed")
    return views


def generate_canonical_view_scripts(
    *,
    output_dir: str | Path,
    session_path: str | Path,
    views: Sequence[CanonicalViewBasis],
    width: int,
    height: int,
) -> CanonicalViewScripts:
    """Generate camera-only scripts without requiring ChimeraX."""

    root = Path(output_dir)
    session = Path(session_path).resolve()
    resolved_views = tuple(views)
    if not session.is_file():
        raise ValueError(f"ChimeraX session does not exist: {session}")
    if len(resolved_views) != 3:
        raise ValueError("canonical view export requires exactly three views")
    if width <= 0 or height <= 0:
        raise ValueError("canonical view dimensions must be positive")
    script_dir = root / "chimerax"
    script_dir.mkdir(parents=True, exist_ok=True)
    python_path = script_dir / "render_canonical_views.py"
    cxc_path = script_dir / "run_render.cxc"
    status_path = root / "logs" / "chimerax_render_status.json"
    config = {
        "session_path": str(session),
        "view_dir": str((root / "views").resolve()),
        "status_path": str(status_path.resolve()),
        "width": int(width),
        "height": int(height),
        "views": [view.manifest_record() for view in resolved_views],
    }
    python_path.write_text(_python_script_text(config), encoding="utf-8")
    cxc_path.write_text(
        f"runscript {_quote_cxc_path(python_path.resolve())}\nexit\n",
        encoding="utf-8",
    )
    return CanonicalViewScripts(python_path, cxc_path, status_path)


def run_canonical_views(args: Any) -> CanonicalViewsResult:
    """Prepare or execute a non-destructive canonical-axis view bundle."""

    run_dir = Path(args.run_dir).resolve()
    session = Path(args.chimerax_session).resolve()
    if not run_dir.is_dir():
        raise ValueError(f"Run directory does not exist: {run_dir}")
    if not session.is_file():
        raise ValueError(f"ChimeraX session does not exist: {session}")
    output_dir = _prepare_output_directory(
        Path(args.output_dir),
        protected=(run_dir, session),
        overwrite=bool(args.overwrite),
    )
    manifest_path = output_dir / "canonical_views.json"
    stage_statuses = {"scripts": "not_started", "render": "not_started"}
    manifest: dict[str, object] = {
        "schema_version": CANONICAL_VIEWS_SCHEMA_VERSION,
        "status": "building",
        "command": "cryorole canonical-views",
        "creation_timestamp": datetime.now(timezone.utc).isoformat(),
        "run_dir": str(run_dir),
        "canonical_id": args.canonical_id,
        "ChimeraX_session": str(session),
        "stage_statuses": stage_statuses,
        "warnings": [],
        "errors": [],
        "resolved_CLI_arguments": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(args).items()
            if key != "handler"
        },
    }
    return_code = None
    try:
        frame = load_canonical_frame(run_dir, args.canonical_id)
        raw_to_scene, map_provenance = _resolve_raw_to_scene(args)
        views = resolve_canonical_view_bases(frame.transform, raw_to_scene)
        scripts = generate_canonical_view_scripts(
            output_dir=output_dir,
            session_path=session,
            views=views,
            width=args.width,
            height=args.height,
        )
        stage_statuses["scripts"] = "complete"
        status = "scripts_ready"
        completion = None
        image_validation = None
        if args.render_mode == "execute":
            if not args.chimerax_bin:
                raise ValueError("--chimerax-bin is required for --render-mode execute")
            result = execute_chimerax(
                chimerax_bin=args.chimerax_bin,
                cxc_script=scripts.cxc_script,
                log_path=output_dir / "logs" / "chimerax_render.log",
            )
            return_code = result.returncode
            completion = validate_canonical_view_completion(
                scripts.status_path,
                expected_count=len(views),
            )
            image_validation = validate_canonical_view_images(
                output_dir / "views",
                expected_names=tuple(view.file_name for view in views),
                expected_dimensions=(args.width, args.height),
            )
            stage_statuses["render"] = "complete"
            status = "rendered"
        else:
            stage_statuses["render"] = "not_requested"
        manifest.update(
            {
                "status": status,
                "canonical_frame": str(frame.source_path),
                "canonical_transform": frame.transform,
                "transform_direction": frame.transform_direction,
                "map_frame": args.map_frame,
                "raw_to_scene": raw_to_scene,
                "raw_to_scene_provenance": map_provenance,
                "canonical_axes_raw": frame.transform,
                "canonical_axes_scene": raw_to_scene @ frame.transform,
                "views": [view.manifest_record() for view in views],
                "camera_policy": "absolute_from_saved_baseline_no_model_motion",
                "image_size": [args.width, args.height],
                "ChimeraX_return_code": return_code,
                "completion": completion,
                "image_validation": image_validation,
                "output_paths": _output_paths(output_dir, scripts, status),
                "stage_statuses": stage_statuses,
            }
        )
        write_animation_manifest(manifest, manifest_path)
        return CanonicalViewsResult(output_dir, manifest_path, status)
    except Exception as exc:
        if isinstance(exc, ChimeraXExecutionError):
            return_code = exc.return_code
        if stage_statuses["scripts"] == "not_started":
            stage_statuses["scripts"] = "failed"
        elif stage_statuses["render"] == "not_started":
            stage_statuses["render"] = "failed"
        manifest.update(
            {
                "status": "failed",
                "stage_statuses": stage_statuses,
                "ChimeraX_return_code": return_code,
                "errors": [f"{type(exc).__name__}: {exc}"],
            }
        )
        write_animation_manifest(manifest, manifest_path)
        raise


def validate_canonical_view_completion(
    path: str | Path,
    *,
    expected_count: int,
) -> dict[str, object]:
    """Require renderer success, unchanged models, and restored camera."""

    status_path = Path(path)
    if not status_path.is_file():
        raise RuntimeError(f"ChimeraX renderer completion artifact is missing: {status_path}")
    try:
        payload = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Invalid ChimeraX completion artifact: {status_path}") from exc
    if not isinstance(payload, dict) or payload.get("status") != "success":
        detail = payload.get("error", "unknown renderer failure") if isinstance(payload, dict) else ""
        raise RuntimeError(f"ChimeraX canonical-view renderer failed: {detail}")
    if payload.get("frame_count") != expected_count:
        raise RuntimeError("ChimeraX canonical-view count does not match expected count")
    if payload.get("models_unchanged") is not True:
        raise RuntimeError("ChimeraX canonical-view renderer changed a model transform")
    if payload.get("camera_restored") is not True:
        raise RuntimeError("ChimeraX canonical-view renderer did not restore the camera")
    return payload


def validate_canonical_view_images(
    directory: str | Path,
    *,
    expected_names: Sequence[str],
    expected_dimensions: tuple[int, int],
) -> dict[str, object]:
    """Validate exact, non-empty, fixed-size canonical-view PNGs."""

    root = Path(directory)
    paths = tuple(sorted(root.glob("*.png"))) if root.is_dir() else ()
    expected = tuple(sorted(expected_names))
    if tuple(path.name for path in paths) != expected:
        raise RuntimeError("Canonical-view PNG names do not match the required three views")
    for path in paths:
        if path.stat().st_size <= 0:
            raise RuntimeError(f"Canonical-view PNG is empty: {path}")
        try:
            with Image.open(path) as image:
                if image.format != "PNG" or image.size != expected_dimensions:
                    raise RuntimeError(f"Canonical-view PNG has invalid format or size: {path}")
                image.verify()
        except (OSError, UnidentifiedImageError) as exc:
            raise RuntimeError(f"Canonical-view output is not a valid PNG: {path}") from exc
    return {
        "count": len(paths),
        "dimensions": expected_dimensions,
        "files": [str(path) for path in paths],
    }


def _resolve_raw_to_scene(args: Any) -> tuple[np.ndarray, dict[str, object]]:
    if args.map_frame == "raw":
        if args.map_frame_transform:
            raise ValueError("--map-frame-transform is only valid with --map-frame explicit")
        return np.eye(3), {"kind": "raw", "source": "identity"}
    if args.map_frame == "explicit":
        if not args.map_frame_transform:
            raise ValueError("--map-frame explicit requires --map-frame-transform")
        return load_explicit_scene_basis(args.map_frame_transform), {
            "kind": "explicit",
            "source": str(Path(args.map_frame_transform).resolve()),
            "column_vector_convention": "x_scene = raw_to_scene @ x_raw",
        }
    raise ValueError("map_frame must be raw or explicit")


def _prepare_output_directory(
    path: Path,
    *,
    protected: Sequence[Path],
    overwrite: bool,
) -> Path:
    output = path.resolve()
    for source in protected:
        resolved = source.resolve()
        if output == resolved or output in resolved.parents:
            raise ValueError("Canonical-view output must not contain a source artifact")
    if output.exists():
        if not overwrite:
            raise FileExistsError(f"Canonical-view output directory already exists: {output}")
        if output == Path(output.anchor) or len(output.parts) < 3:
            raise ValueError(f"Refusing unsafe canonical-view output replacement: {output}")
        shutil.rmtree(output)
    output.mkdir(parents=True)
    return output


def _output_paths(
    output_dir: Path,
    scripts: CanonicalViewScripts,
    status: str,
) -> dict[str, object]:
    return {
        "manifest": str(output_dir / "canonical_views.json"),
        "chimerax_python": str(scripts.python_script),
        "chimerax_cxc": str(scripts.cxc_script),
        "chimerax_render_log": (
            str(output_dir / "logs" / "chimerax_render.log")
            if status == "rendered"
            else None
        ),
        "chimerax_render_status": (
            str(scripts.status_path) if status == "rendered" else None
        ),
        "views": str(output_dir / "views") if status == "rendered" else None,
    }


def _quote_cxc_path(path: Path) -> str:
    return '"' + str(path).replace("\\", "/").replace('"', '\\"') + '"'


def _python_script_text(config: dict[str, object]) -> str:
    config_json = json.dumps(config, ensure_ascii=False, default=_json_default)
    return f'''"""Generated by cryorole canonical-views; camera-only rendering."""
import json
import traceback
from pathlib import Path

import numpy as np
from chimerax.core.commands import run
from chimerax.geometry import Place

CONFIG = json.loads({config_json!r})
STATE = {{"frame_count": 0, "models_unchanged": False, "camera_restored": False}}


def quote_command_path(value):
    return '"' + str(value).replace("\\\\", "/").replace('"', '\\\\"') + '"'


def write_status(payload):
    path = Path(CONFIG["status_path"])
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    temporary.replace(path)


def unchanged_models(baselines):
    return all(
        np.allclose(model.scene_position.matrix, baseline)
        for model, baseline in baselines
    )


def render():
    run(session, "open " + quote_command_path(CONFIG["session_path"]))
    main_view = session.main_view
    camera = main_view.camera
    camera_baseline = camera.position
    camera_matrix = np.asarray(camera_baseline.matrix, dtype=float).copy()
    center = np.asarray(main_view.center_of_rotation, dtype=float)
    if camera_matrix.shape != (3, 4) or center.shape != (3,):
        raise RuntimeError("Unexpected ChimeraX camera or center-of-rotation shape")
    distance = float(np.linalg.norm(camera_matrix[:, 3] - center))
    if not np.isfinite(distance) or distance <= 0:
        raise RuntimeError("Cannot resolve a positive camera distance from scene center")
    model_baselines = [
        (model, np.asarray(model.scene_position.matrix, dtype=float).copy())
        for model in session.models.list()
        if hasattr(model, "scene_position")
    ]
    output_dir = Path(CONFIG["view_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Camera Place maps camera coordinates to scene coordinates. Its first,
    # second, and third columns are screen-right, screen-up, and toward-viewer.
    # Every view is absolute from the saved camera baseline and scene center.
    try:
        for view in CONFIG["views"]:
            right = np.asarray(view["screen_right"], dtype=float)
            up = np.asarray(view["screen_up"], dtype=float)
            outward = np.asarray(view["toward_viewer"], dtype=float)
            rotation = np.column_stack((right, up, outward))
            origin = center + distance * outward
            matrix = np.column_stack((rotation, origin))
            camera_place = Place(matrix=tuple(tuple(float(value) for value in row) for row in matrix))
            camera.position = camera_place
            output_path = output_dir / view["file_name"]
            run(
                session,
                "save {{}} width {{}} height {{}} supersample 1".format(
                    quote_command_path(output_path),
                    CONFIG["width"],
                    CONFIG["height"],
                ),
            )
            STATE["frame_count"] += 1
            if not unchanged_models(model_baselines):
                raise RuntimeError("A model scene transform changed during camera rendering")
        STATE["models_unchanged"] = bool(unchanged_models(model_baselines))
    finally:
        camera.position = camera_baseline
        STATE["camera_restored"] = bool(
            np.allclose(camera.position.matrix, camera_matrix)
        )
        STATE["models_unchanged"] = bool(unchanged_models(model_baselines))
    if not STATE["models_unchanged"]:
        raise RuntimeError("Model scene transforms were not preserved")
    if not STATE["camera_restored"]:
        raise RuntimeError("Camera baseline was not restored")


try:
    render()
except Exception as error:
    write_status({{
        "status": "failed",
        "error": "{{}}: {{}}".format(type(error).__name__, error),
        "traceback": traceback.format_exc(),
        **STATE,
    }})
    raise
else:
    write_status({{"status": "success", **STATE}})
'''


def _json_default(value: object):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")
