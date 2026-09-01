"""ChimeraX script generation and optional headless execution."""

from __future__ import annotations

import csv
import json
import struct
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

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
    reference_model_id: str | Sequence[str],
    moving_model_id: str | Sequence[str],
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
    reference_ids, moving_ids = normalize_model_id_groups(
        reference_model_id,
        moving_model_id,
    )
    if structure_width <= 0 or structure_height <= 0:
        raise ValueError("structure frame dimensions must be positive")
    if view_name not in {None, "primary", "secondary", "tertiary"}:
        raise ValueError("view_name must be primary, secondary, tertiary, or None")
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
        "reference_model_ids": reference_ids,
        "moving_model_ids": moving_ids,
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
    expected_reference_model_ids: str | Sequence[str] | None = None,
    expected_moving_model_ids: str | Sequence[str] | None = None,
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
    if expected_reference_model_ids is not None or expected_moving_model_ids is not None:
        reference_ids, moving_ids = normalize_model_id_groups(
            expected_reference_model_ids or (),
            expected_moving_model_ids or (),
        )
        strict_groups = (
            len(reference_ids) > 1
            or len(moving_ids) > 1
            or "reference_model_ids" in payload
            or "moving_model_ids" in payload
        )
        if strict_groups:
            _validate_completion_model_group(
                payload,
                group_name="reference",
                expected_ids=reference_ids,
                result_field="reference_unchanged_by_model",
            )
            _validate_completion_model_group(
                payload,
                group_name="moving",
                expected_ids=moving_ids,
                result_field="moving_restored_by_model",
            )
            if payload.get("stationary_unchanged") is not True:
                raise RuntimeError(
                    "ChimeraX renderer did not confirm unchanged stationary transforms"
                )
            stationary_checks = payload.get("stationary_unchanged_by_model")
            if not isinstance(stationary_checks, dict) or not all(
                value is True for value in stationary_checks.values()
            ):
                raise RuntimeError(
                    "ChimeraX renderer did not confirm every stationary model transform"
                )
            if not set(reference_ids).issubset(map(str, stationary_checks)):
                raise RuntimeError(
                    "ChimeraX renderer stationary checks omit a reference model"
                )
            relative_checks = payload.get("moving_group_relative_transform_checks")
            expected_relative_checks = len(moving_ids) * (len(moving_ids) - 1) // 2
            if (
                not isinstance(relative_checks, dict)
                or len(relative_checks) != expected_relative_checks
                or not all(value is True for value in relative_checks.values())
            ):
                raise RuntimeError(
                    "ChimeraX renderer did not validate every moving-group relative transform"
                )
            if payload.get("moving_group_relative_transforms_preserved") is not True:
                raise RuntimeError(
                    "ChimeraX renderer did not preserve moving-group relative transforms"
                )
    return payload


def validate_matching_session_baselines(
    primary: dict[str, object],
    secondary: dict[str, object],
    tertiary: dict[str, object] | None = None,
    *,
    atol: float = 1e-8,
) -> dict[str, object]:
    """Require all requested sessions to start from identical model transforms."""

    if not np.isfinite(atol) or atol <= 0:
        raise ValueError("session baseline comparison tolerance must be positive")
    comparisons = [("secondary", secondary)]
    if tertiary is not None:
        comparisons.append(("tertiary", tertiary))
    primary_groups = _completion_baseline_groups(primary)
    differences: dict[str, float] = {}
    per_view_differences: dict[str, dict[str, float]] = {}
    for view_name, candidate in comparisons:
        candidate_groups = _completion_baseline_groups(candidate)
        view_differences: dict[str, float] = {}
        for group_name in ("reference", "moving"):
            primary_group = primary_groups[group_name]
            candidate_group = candidate_groups[group_name]
            if tuple(primary_group) != tuple(candidate_group):
                raise RuntimeError(
                    f"{view_name} ChimeraX session {group_name} model IDs differ"
                )
            for model_id, primary_matrix in primary_group.items():
                candidate_matrix = candidate_group[model_id]
                name = f"{group_name} model #{model_id}"
                difference = float(np.max(np.abs(primary_matrix - candidate_matrix)))
                view_differences[name] = difference
                differences[name] = max(differences.get(name, 0.0), difference)
                if not np.allclose(
                    primary_matrix,
                    candidate_matrix,
                    rtol=0.0,
                    atol=atol,
                ):
                    raise RuntimeError(
                        f"{view_name} ChimeraX session initial scene transforms differ "
                        f"for {name}"
                    )
        per_view_differences[view_name] = view_differences
    return {
        "status": "validated",
        "matched": True,
        "view_count": len(comparisons) + 1,
        "absolute_tolerance": float(atol),
        "maximum_absolute_differences": differences,
        "per_view_maximum_absolute_differences": per_view_differences,
        "reference_model_ids": list(primary_groups["reference"]),
        "moving_model_ids": list(primary_groups["moving"]),
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


def normalize_model_id_groups(
    reference_model_ids: str | Sequence[str],
    moving_model_ids: str | Sequence[str],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Normalize exact ChimeraX IDs and validate disjoint rigid groups."""

    reference = _normalize_model_id_group(reference_model_ids, "reference")
    moving = _normalize_model_id_group(moving_model_ids, "moving")
    overlap = set(reference).intersection(moving)
    if overlap:
        joined = ", ".join(f"#{model_id}" for model_id in sorted(overlap))
        raise ValueError(f"model IDs cannot be both reference and moving: {joined}")
    return reference, moving


def _normalize_model_id_group(
    values: str | Sequence[str],
    group_name: str,
) -> tuple[str, ...]:
    items = (values,) if isinstance(values, str) else tuple(values)
    if not items:
        raise ValueError(f"at least one {group_name} ChimeraX model ID is required")
    normalized = tuple(_normalize_model_id(value) for value in items)
    if len(set(normalized)) != len(normalized):
        raise ValueError(f"duplicate {group_name} ChimeraX model IDs are not allowed")
    return normalized


def _validate_completion_model_group(
    payload: dict[str, object],
    *,
    group_name: str,
    expected_ids: tuple[str, ...],
    result_field: str,
) -> None:
    recorded_ids = payload.get(f"{group_name}_model_ids")
    if recorded_ids != list(expected_ids):
        raise RuntimeError(
            f"ChimeraX renderer {group_name} model IDs do not match the request"
        )
    baselines = payload.get(f"{group_name}_baseline_matrices")
    if not isinstance(baselines, dict) or tuple(baselines) != expected_ids:
        raise RuntimeError(
            f"ChimeraX renderer did not record every {group_name} model baseline"
        )
    for model_id, matrix in baselines.items():
        values = np.asarray(matrix, dtype=float)
        if values.shape != (3, 4) or not np.isfinite(values).all():
            raise RuntimeError(
                f"ChimeraX renderer recorded an invalid {group_name} model #{model_id} baseline"
            )
    results = payload.get(result_field)
    if not isinstance(results, dict) or tuple(results) != expected_ids:
        raise RuntimeError(
            f"ChimeraX renderer did not check every {group_name} model"
        )
    failed = [model_id for model_id, value in results.items() if value is not True]
    if failed:
        raise RuntimeError(
            f"ChimeraX renderer check failed for {group_name} model #{failed[0]}"
        )


def _completion_baseline_groups(
    payload: dict[str, object],
) -> dict[str, dict[str, np.ndarray]]:
    groups: dict[str, dict[str, np.ndarray]] = {}
    for group_name in ("reference", "moving"):
        plural = payload.get(f"{group_name}_baseline_matrices")
        if plural is None:
            plural = {"legacy": payload.get(f"{group_name}_baseline_matrix")}
        if not isinstance(plural, dict) or not plural:
            raise RuntimeError(
                f"ChimeraX completion artifacts do not contain valid {group_name} baselines"
            )
        matrices: dict[str, np.ndarray] = {}
        for model_id, matrix in plural.items():
            try:
                values = np.asarray(matrix, dtype=float)
            except (TypeError, ValueError) as exc:
                raise RuntimeError(
                    f"ChimeraX completion artifacts do not contain valid {group_name} baselines"
                ) from exc
            if values.shape != (3, 4) or not np.isfinite(values).all():
                raise RuntimeError(
                    f"ChimeraX completion artifacts do not contain valid {group_name} baselines"
                )
            matrices[str(model_id)] = values
        groups[group_name] = matrices
    return groups


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
    "reference_model_ids": list(CONFIG["reference_model_ids"]),
    "moving_model_ids": list(CONFIG["moving_model_ids"]),
    "reference_baseline_matrices": {{}},
    "moving_baseline_matrices": {{}},
    "reference_unchanged_by_model": {{}},
    "stationary_unchanged_by_model": {{}},
    "moving_restored_by_model": {{}},
    "moving_group_relative_transform_checks": {{}},
    "reference_unchanged": False,
    "stationary_unchanged": False,
    "moving_restored": False,
    "moving_group_relative_transforms_preserved": False,
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


def ancestor_in_group(model, object_ids):
    parent = getattr(model, "parent", None)
    while parent is not None:
        if id(parent) in object_ids:
            return parent
        parent = getattr(parent, "parent", None)
    return None


def baseline_entries(models):
    return [
        (
            model,
            model.scene_position,
            np.asarray(model.scene_position.matrix, dtype=float).copy(),
        )
        for model in models
    ]


def unchanged_by_model(entries):
    return {{
        model.id_string: bool(np.allclose(model.scene_position.matrix, matrix))
        for model, _place, matrix in entries
    }}


def affine_matrix(place):
    matrix = np.asarray(place.matrix, dtype=float)
    affine = np.eye(4, dtype=float)
    affine[:3, :] = matrix
    return affine


def relative_matrices(entries):
    result = {{}}
    for left_index, (left, left_place, _left_matrix) in enumerate(entries):
        for right, right_place, _right_matrix in entries[left_index + 1:]:
            key = "#{{}} -> #{{}}".format(left.id_string, right.id_string)
            result[key] = np.linalg.inv(affine_matrix(left_place)) @ affine_matrix(right_place)
    return result


def relative_transform_checks(entries, baselines):
    current_entries = [
        (model, model.scene_position, matrix)
        for model, _place, matrix in entries
    ]
    current = relative_matrices(current_entries)
    return {{
        key: bool(np.allclose(current[key], baseline))
        for key, baseline in baselines.items()
    }}


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
    moving_baselines = []
    run(session, "open " + quote_command_path(CONFIG["session_path"]))
    references = [exact_model(model_id) for model_id in CONFIG["reference_model_ids"]]
    movings = [exact_model(model_id) for model_id in CONFIG["moving_model_ids"]]
    if set(map(id, references)).intersection(map(id, movings)):
        raise RuntimeError("reference and moving model IDs resolve to overlapping models")
    reference_baselines = baseline_entries(references)
    moving_baselines = baseline_entries(movings)
    moving_object_ids = set(map(id, movings))
    for moving in movings:
        ancestor = ancestor_in_group(moving, moving_object_ids)
        if ancestor is not None:
            raise RuntimeError(
                "moving model IDs cannot include both parent #{{}} and descendant #{{}}".format(
                    ancestor.id_string, moving.id_string
                )
            )
    for reference in references:
        ancestor = ancestor_in_group(reference, moving_object_ids)
        if ancestor is not None:
            raise RuntimeError(
                "reference model #{{}} cannot be a descendant of moving model #{{}}".format(
                    reference.id_string, ancestor.id_string
                )
            )
    stationary_baselines = baseline_entries([
        model for model in session.models.list()
        if (
            hasattr(model, "scene_position")
            and id(model) not in moving_object_ids
            and ancestor_in_group(model, moving_object_ids) is None
        )
    ])
    moving_relative_baselines = relative_matrices(moving_baselines)
    STATE["reference_baseline_matrices"] = {{
        model.id_string: matrix.tolist()
        for model, _place, matrix in reference_baselines
    }}
    STATE["moving_baseline_matrices"] = {{
        model.id_string: matrix.tolist()
        for model, _place, matrix in moving_baselines
    }}
    if len(reference_baselines) == 1:
        STATE["reference_baseline_matrix"] = reference_baselines[0][2].tolist()
    if len(moving_baselines) == 1:
        STATE["moving_baseline_matrix"] = moving_baselines[0][2].tolist()
    STATE["moving_group_relative_transforms_preserved"] = True
    output_dir = Path(CONFIG["structure_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    # ChimeraX Place maps local coordinates to scene coordinates. Place multiplication
    # acts right-to-left, so a scene-space pivot delta left-multiplies the saved
    # local-to-scene baseline: absolute = delta_place * moving_baseline.
    try:
        for frame_index, matrix in enumerate(read_transforms(CONFIG["transform_csv"])):
            delta_place = Place(matrix=matrix)
            for moving, moving_baseline, _matrix in moving_baselines:
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
            STATE["reference_unchanged_by_model"] = unchanged_by_model(reference_baselines)
            STATE["stationary_unchanged_by_model"] = unchanged_by_model(stationary_baselines)
            STATE["moving_group_relative_transform_checks"] = relative_transform_checks(
                moving_baselines, moving_relative_baselines
            )
            STATE["reference_unchanged"] = bool(
                all(STATE["reference_unchanged_by_model"].values())
            )
            STATE["stationary_unchanged"] = bool(
                all(STATE["stationary_unchanged_by_model"].values())
            )
            if not STATE["reference_unchanged"]:
                raise RuntimeError("reference scene transform changed during rendering")
            if not STATE["stationary_unchanged"]:
                raise RuntimeError("an undeclared stationary scene transform changed during rendering")
            if not all(STATE["moving_group_relative_transform_checks"].values()):
                STATE["moving_group_relative_transforms_preserved"] = False
                raise RuntimeError("moving-group relative transforms changed during rendering")
    finally:
        for moving, moving_baseline, _matrix in moving_baselines:
            moving.scene_position = moving_baseline
        STATE["moving_restored_by_model"] = unchanged_by_model(moving_baselines)
        STATE["moving_restored"] = bool(
            moving_baselines and all(STATE["moving_restored_by_model"].values())
        )
    if not STATE["reference_unchanged"]:
        raise RuntimeError("reference scene transform changed during rendering")
    if not STATE["stationary_unchanged"]:
        raise RuntimeError("an undeclared stationary scene transform changed during rendering")
    if not STATE["moving_restored"]:
        raise RuntimeError("moving scene transform was not restored")
    if not STATE["moving_group_relative_transforms_preserved"]:
        raise RuntimeError("moving-group relative transforms were not preserved")


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
