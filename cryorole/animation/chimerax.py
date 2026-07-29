"""ChimeraX script generation and optional headless execution."""

from __future__ import annotations

import csv
import json
import struct
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


TRANSFORM_COLUMNS = ("frame_index",) + tuple(
    f"m{row}{column}" for row in range(3) for column in range(4)
)


@dataclass(frozen=True)
class ChimeraXScripts:
    python_script: Path
    cxc_script: Path
    transform_csv: Path
    status_path: Path


@dataclass(frozen=True)
class StructureFrameValidation:
    frame_count: int
    dimensions: tuple[int, int]
    paths: tuple[Path, ...]


class ChimeraXExecutionError(RuntimeError):
    """A ChimeraX process failure carrying its exact return code."""

    def __init__(self, return_code: int, log_path: Path) -> None:
        self.return_code = int(return_code)
        self.log_path = log_path
        super().__init__(
            f"ChimeraX rendering failed with return code {return_code}; see {log_path}"
        )


def write_scene_transforms_csv(
    transforms: np.ndarray,
    path: str | Path,
) -> Path:
    """Persist pivoted scene deltas as contiguous column-vector 3x4 matrices."""

    values = np.asarray(transforms, dtype=float)
    if values.ndim != 3 or values.shape[1:] != (4, 4):
        raise ValueError("scene transforms must have shape (n, 4, 4)")
    if len(values) == 0 or not np.isfinite(values).all():
        raise ValueError("scene transforms must be non-empty and finite")
    if not np.allclose(values[:, 3, :], np.asarray([0.0, 0.0, 0.0, 1.0])):
        raise ValueError("scene transforms must be affine homogeneous matrices")
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRANSFORM_COLUMNS)
        writer.writeheader()
        for index, transform in enumerate(values):
            row = {"frame_index": index}
            row.update(
                {
                    f"m{row_index}{column_index}": transform[row_index, column_index]
                    for row_index in range(3)
                    for column_index in range(4)
                }
            )
            writer.writerow(row)
    return output


def generate_chimerax_scripts(
    *,
    output_dir: str | Path,
    session_path: str | Path,
    transform_csv: str | Path,
    reference_model_id: str,
    moving_model_id: str,
    structure_width: int,
    structure_height: int,
    view_name: str | None = None,
) -> ChimeraXScripts:
    """Generate auditable scripts without requiring a ChimeraX installation."""

    root = Path(output_dir)
    session = Path(session_path).resolve()
    transforms = Path(transform_csv).resolve()
    if not session.is_file():
        raise ValueError(f"ChimeraX session does not exist: {session}")
    if not transforms.is_file():
        raise ValueError(f"Scene transform artifact does not exist: {transforms}")
    reference = _normalize_model_id(reference_model_id)
    moving = _normalize_model_id(moving_model_id)
    if reference == moving:
        raise ValueError("reference and moving ChimeraX model IDs must differ")
    if structure_width <= 0 or structure_height <= 0:
        raise ValueError("structure frame dimensions must be positive")
    if view_name not in {None, "primary", "secondary"}:
        raise ValueError("view_name must be primary, secondary, or None")
    script_dir = root / "chimerax" / view_name if view_name else root / "chimerax"
    structure_dir = (
        root / "frames" / "structure" / view_name
        if view_name
        else root / "frames" / "structure"
    )
    script_dir.mkdir(parents=True, exist_ok=True)
    python_path = script_dir / "render_structure.py"
    cxc_path = script_dir / "run_render.cxc"
    status_path = root / "logs" / (
        f"chimerax_render_{view_name}_status.json"
        if view_name
        else "chimerax_render_status.json"
    )
    config = {
        "session_path": str(session),
        "transform_csv": str(transforms),
        "structure_dir": str(structure_dir.resolve()),
        "reference_model_id": reference,
        "moving_model_id": moving,
        "structure_width": int(structure_width),
        "structure_height": int(structure_height),
        "status_path": str(status_path.resolve()),
    }
    python_path.write_text(_python_script_text(config), encoding="utf-8")
    cxc_path.write_text(
        f"runscript {_quote_cxc_path(python_path.resolve())}\nexit\n",
        encoding="utf-8",
    )
    return ChimeraXScripts(
        python_script=python_path,
        cxc_script=cxc_path,
        transform_csv=transforms,
        status_path=status_path,
    )


def execute_chimerax(
    *,
    chimerax_bin: str | Path,
    cxc_script: str | Path,
    log_path: str | Path,
    platform_name: str | None = None,
) -> subprocess.CompletedProcess[str]:
    """Execute ChimeraX with an argument list and capture a structured log."""

    executable = Path(chimerax_bin)
    script = Path(cxc_script)
    if not executable.is_file():
        raise ValueError(f"ChimeraX executable does not exist: {executable}")
    if not script.is_file():
        raise ValueError(f"ChimeraX CXC script does not exist: {script}")
    platform_value = (platform_name or sys.platform).lower()
    if platform_value.startswith("linux"):
        arguments = [
            str(executable),
            "--offscreen",
            "--exit",
            "--nocolor",
            "--script",
            str(script.resolve()),
        ]
    else:
        arguments = [
            str(executable),
            "--nogui",
            "--exit",
            "--nocolor",
            "--script",
            str(script.resolve()),
        ]
    result = subprocess.run(
        arguments,
        shell=False,
        capture_output=True,
        text=True,
        check=False,
    )
    log = Path(log_path)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(
        json.dumps(
            {
                "arguments": arguments,
                "return_code": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    if result.returncode != 0:
        raise ChimeraXExecutionError(result.returncode, log)
    return result


def validate_chimerax_completion(
    path: str | Path,
    *,
    expected_count: int,
) -> dict[str, object]:
    """Require explicit renderer completion independent of the process return code."""

    status_path = Path(path)
    if not status_path.is_file():
        raise RuntimeError(f"ChimeraX renderer completion artifact is missing: {status_path}")
    try:
        payload = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"ChimeraX renderer completion artifact is invalid: {status_path}"
        ) from exc
    if not isinstance(payload, dict):
        raise RuntimeError("ChimeraX renderer completion artifact must contain an object")
    if payload.get("status") != "success":
        detail = payload.get("error") or payload.get("traceback") or "unknown renderer failure"
        raise RuntimeError(f"ChimeraX renderer reported failure: {detail}")
    if payload.get("frame_count") != expected_count:
        raise RuntimeError(
            "ChimeraX renderer completion frame count does not match trajectory: "
            f"{payload.get('frame_count')} != {expected_count}"
        )
    if payload.get("reference_unchanged") is not True:
        raise RuntimeError("ChimeraX renderer did not confirm unchanged reference transform")
    if payload.get("moving_restored") is not True:
        raise RuntimeError("ChimeraX renderer did not confirm moving-transform restoration")
    return payload


def validate_matching_session_baselines(
    primary: dict[str, object],
    secondary: dict[str, object],
    *,
    atol: float = 1e-8,
) -> dict[str, object]:
    """Require both dual-view sessions to start from identical model transforms."""

    if not np.isfinite(atol) or atol <= 0:
        raise ValueError("session baseline comparison tolerance must be positive")
    differences: dict[str, float] = {}
    for name in ("reference_baseline_matrix", "moving_baseline_matrix"):
        try:
            primary_matrix = np.asarray(primary[name], dtype=float)
            secondary_matrix = np.asarray(secondary[name], dtype=float)
        except (KeyError, TypeError, ValueError) as exc:
            raise RuntimeError(
                f"ChimeraX completion artifacts do not contain valid {name}"
            ) from exc
        if (
            primary_matrix.shape != (3, 4)
            or secondary_matrix.shape != (3, 4)
            or not np.isfinite(primary_matrix).all()
            or not np.isfinite(secondary_matrix).all()
        ):
            raise RuntimeError(
                f"ChimeraX completion artifacts do not contain valid {name}"
            )
        difference = float(np.max(np.abs(primary_matrix - secondary_matrix)))
        differences[name] = difference
        if not np.allclose(primary_matrix, secondary_matrix, rtol=0.0, atol=atol):
            raise RuntimeError(
                "Dual ChimeraX session initial scene transforms differ "
                f"for {name}"
            )
    return {
        "status": "validated",
        "matched": True,
        "absolute_tolerance": float(atol),
        "maximum_absolute_differences": differences,
    }


def validate_structure_frames(
    directory: str | Path,
    *,
    expected_count: int,
) -> StructureFrameValidation:
    """Validate contiguous, non-empty, fixed-size ChimeraX PNG output."""

    root = Path(directory)
    paths = tuple(sorted(root.glob("frame_*.png"))) if root.is_dir() else ()
    expected_names = tuple(f"frame_{index:06d}.png" for index in range(expected_count))
    names = tuple(path.name for path in paths)
    if names != expected_names:
        raise RuntimeError(
            "ChimeraX structure frame indices are not contiguous or do not match "
            f"trajectory count {expected_count}: {names}"
        )
    dimensions: set[tuple[int, int]] = set()
    for path in paths:
        if path.stat().st_size <= 0:
            raise RuntimeError(f"ChimeraX structure frame is empty: {path}")
        dimensions.add(_png_dimensions(path))
    if len(dimensions) != 1:
        raise RuntimeError("ChimeraX structure frame dimensions are inconsistent")
    return StructureFrameValidation(
        frame_count=len(paths),
        dimensions=next(iter(dimensions)),
        paths=paths,
    )


def _png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError(f"Structure frame is not a valid PNG: {path}")
    width, height = struct.unpack(">II", header[16:24])
    return int(width), int(height)


def _normalize_model_id(value: str) -> str:
    normalized = str(value).strip()
    if normalized.startswith("#"):
        normalized = normalized[1:]
    if not normalized:
        raise ValueError("ChimeraX model ID must be non-empty")
    return normalized


def _quote_cxc_path(path: Path) -> str:
    return '"' + str(path).replace("\\", "/").replace('"', '\\"') + '"'


def _python_script_text(config: dict[str, object]) -> str:
    config_json = json.dumps(config, ensure_ascii=False)
    return f'''"""Generated by cryorole animate; do not edit scientific transform formulas here."""
import csv
import json
import traceback
from pathlib import Path

import numpy as np
from chimerax.core.commands import run
from chimerax.geometry import Place

CONFIG = json.loads({config_json!r})
STATE = {{
    "frame_count": 0,
    "reference_unchanged": False,
    "moving_restored": False,
}}


def quote_command_path(value):
    return '"' + str(value).replace("\\\\", "/").replace('"', '\\\\"') + '"'


def exact_model(model_id):
    requested = str(model_id).lstrip("#")
    available = [(model.id_string, getattr(model, "name", "")) for model in session.models.list()]
    matches = [
        model for model in session.models.list()
        if model.id_string == requested
    ]
    if len(matches) != 1:
        raise RuntimeError(
            "Expected exactly one model with id #{{}}, available={{}}".format(requested, available)
        )
    return matches[0]


def read_transforms(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for expected_index, row in enumerate(csv.DictReader(handle)):
            if int(row["frame_index"]) != expected_index:
                raise RuntimeError("transform frame indices are not contiguous")
            rows.append(tuple(
                tuple(float(row["m{{}}{{}}".format(r, c)]) for c in range(4))
                for r in range(3)
            ))
    return rows


def write_status(payload):
    path = Path(CONFIG["status_path"])
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    temporary.replace(path)


def render():
    moving = None
    moving_baseline = None
    run(session, "open " + quote_command_path(CONFIG["session_path"]))
    reference = exact_model(CONFIG["reference_model_id"])
    moving = exact_model(CONFIG["moving_model_id"])
    if reference is moving:
        raise RuntimeError("reference and moving model IDs resolve to the same model")
    reference_baseline = reference.scene_position
    moving_baseline = moving.scene_position
    STATE["reference_baseline_matrix"] = np.asarray(
        reference_baseline.matrix, dtype=float
    ).tolist()
    STATE["moving_baseline_matrix"] = np.asarray(
        moving_baseline.matrix, dtype=float
    ).tolist()
    output_dir = Path(CONFIG["structure_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    # ChimeraX Place maps local coordinates to scene coordinates. Place multiplication
    # acts right-to-left, so a scene-space pivot delta left-multiplies the saved
    # local-to-scene baseline: absolute = delta_place * moving_baseline.
    try:
        for frame_index, matrix in enumerate(read_transforms(CONFIG["transform_csv"])):
            delta_place = Place(matrix=matrix)
            moving.scene_position = delta_place * moving_baseline
            output_path = output_dir / "frame_{{:06d}}.png".format(frame_index)
            run(
                session,
                "save {{}} width {{}} height {{}} supersample 1".format(
                    quote_command_path(output_path),
                    CONFIG["structure_width"],
                    CONFIG["structure_height"],
                ),
            )
            STATE["frame_count"] = frame_index + 1
            if not np.allclose(reference.scene_position.matrix, reference_baseline.matrix):
                raise RuntimeError("reference scene transform changed during rendering")
        STATE["reference_unchanged"] = bool(
            np.allclose(reference.scene_position.matrix, reference_baseline.matrix)
        )
    finally:
        if moving is not None and moving_baseline is not None:
            moving.scene_position = moving_baseline
            STATE["moving_restored"] = bool(
                np.allclose(moving.scene_position.matrix, moving_baseline.matrix)
            )
    if not STATE["reference_unchanged"]:
        raise RuntimeError("reference scene transform changed during rendering")
    if not STATE["moving_restored"]:
        raise RuntimeError("moving scene transform was not restored")


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
