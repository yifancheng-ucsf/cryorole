"""Waypoint parsing and deterministic SO(3) trajectory generation."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy.spatial.transform import Rotation, Slerp

from cryorole.animation.schemas import TrajectoryFrame, Waypoint


EA_COLUMNS = ("alpha_deg", "beta_deg", "gamma_deg")
RV_COLUMNS = ("rv_x_rad", "rv_y_rad", "rv_z_rad")
TRAJECTORY_COLUMNS = (
    "frame_index",
    "time_sec",
    "segment_index",
    "from_label",
    "to_label",
    "segment_fraction",
    "quat_w",
    "quat_x",
    "quat_y",
    "quat_z",
    "rv_x_rad",
    "rv_y_rad",
    "rv_z_rad",
    "alpha_deg",
    "beta_deg",
    "gamma_deg",
    "is_waypoint",
    "waypoint_label",
    "annotation",
)


def parse_waypoint_csv(
    path: str | Path,
    *,
    path_space: str,
    scipy_euler_sequence: str,
    default_segment_frames: int,
    default_hold_frames: int,
) -> tuple[tuple[Waypoint, ...], tuple[str, ...]]:
    """Read and validate an EA or RV waypoint CSV without unit inference."""

    if path_space not in {"ea", "rv"}:
        raise ValueError("path_space must be 'ea' or 'rv'")
    _validate_integer(default_segment_frames, "frames_per_segment", minimum=2)
    _validate_integer(default_hold_frames, "hold_frames", minimum=0)
    input_path = Path(path)
    if not input_path.is_file():
        raise ValueError(f"Waypoint CSV does not exist: {input_path}")
    with input_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = tuple(reader.fieldnames or ())
        required = ("label",) + (EA_COLUMNS if path_space == "ea" else RV_COLUMNS)
        alternate = RV_COLUMNS if path_space == "ea" else EA_COLUMNS
        missing = [column for column in required if column not in fieldnames]
        if missing:
            raise ValueError(f"Waypoint CSV missing required columns: {', '.join(missing)}")
        if any(column in fieldnames for column in alternate):
            raise ValueError("Waypoint CSV must not mix EA and RV schemas")
        rows = list(reader)
    if len(rows) < 2:
        raise ValueError("Waypoint CSV must contain at least two waypoints")

    labels: set[str] = set()
    warnings: list[str] = []
    points: list[Waypoint] = []
    coordinate_columns = EA_COLUMNS if path_space == "ea" else RV_COLUMNS
    for index, row in enumerate(rows):
        label = str(row.get("label", "")).strip()
        if not label:
            raise ValueError(f"Waypoint row {index + 2} has an empty label")
        if label in labels:
            raise ValueError(f"Waypoint labels must be unique; duplicate {label!r}")
        labels.add(label)
        coordinates = np.asarray(
            [_parse_float(row.get(column), column, index + 2) for column in coordinate_columns],
            dtype=float,
        )
        if not np.isfinite(coordinates).all():
            raise ValueError(f"Waypoint {label!r} coordinates must be finite")
        rotation = (
            Rotation.from_euler(scipy_euler_sequence, coordinates, degrees=True)
            if path_space == "ea"
            else Rotation.from_rotvec(coordinates)
        )
        segment_value = _optional_integer(
            row.get("segment_frames"),
            default_segment_frames,
            "segment_frames",
            index + 2,
            minimum=2,
        )
        hold_value = _optional_integer(
            row.get("hold_frames"),
            default_hold_frames,
            "hold_frames",
            index + 2,
            minimum=0,
        )
        if index == len(rows) - 1 and _has_value(row.get("segment_frames")):
            warnings.append(
                f"Final waypoint {label!r} segment_frames is ignored because it has no outgoing segment"
            )
        points.append(
            Waypoint(
                label=label,
                rotation=rotation,
                segment_frames=segment_value,
                hold_frames=hold_value,
                annotation=str(row.get("annotation", "") or "").strip(),
            )
        )
    return tuple(points), tuple(warnings)


def reverse_waypoints(waypoints: Sequence[Waypoint]) -> tuple[Waypoint, ...]:
    """Reverse path direction while retaining each segment's frame count."""

    points = tuple(waypoints)
    reversed_points: list[Waypoint] = []
    for new_index, point in enumerate(reversed(points)):
        original_segment_index = len(points) - 2 - new_index
        segment_frames = (
            points[original_segment_index].segment_frames
            if original_segment_index >= 0
            else point.segment_frames
        )
        reversed_points.append(
            Waypoint(
                label=point.label,
                rotation=point.rotation,
                segment_frames=segment_frames,
                hold_frames=point.hold_frames,
                annotation=point.annotation,
            )
        )
    return tuple(reversed_points)


def continuous_quaternions(quaternions_xyzw: np.ndarray) -> np.ndarray:
    """Normalize SciPy-order quaternions and enforce adjacent sign continuity."""

    values = np.asarray(quaternions_xyzw, dtype=float).copy()
    if values.ndim != 2 or values.shape[1] != 4:
        raise ValueError("quaternions must have shape (n, 4)")
    if not np.isfinite(values).all():
        raise ValueError("quaternions contain non-finite values")
    norms = np.linalg.norm(values, axis=1)
    if np.any(norms <= 0):
        raise ValueError("quaternions must have non-zero norm")
    values /= norms[:, None]
    for index in range(1, len(values)):
        if float(np.dot(values[index - 1], values[index])) < 0:
            values[index] *= -1.0
    return values


def build_trajectory(
    waypoints: Sequence[Waypoint],
    *,
    fps: float,
    ping_pong: bool = False,
) -> tuple[TrajectoryFrame, ...]:
    """Interpolate validated waypoints with endpoint-inclusive SO(3) SLERP."""

    points = tuple(waypoints)
    if len(points) < 2:
        raise ValueError("At least two waypoints are required")
    if not np.isfinite(fps) or fps <= 0:
        raise ValueError("fps must be a finite positive number")
    waypoint_quaternions = continuous_quaternions(
        np.vstack([point.rotation.as_quat() for point in points])
    )
    normalized_points = tuple(
        Waypoint(
            label=point.label,
            rotation=Rotation.from_quat(quaternion),
            segment_frames=point.segment_frames,
            hold_frames=point.hold_frames,
            annotation=point.annotation,
        )
        for point, quaternion in zip(points, waypoint_quaternions)
    )

    records: list[dict[str, object]] = []
    for segment_index, (start, stop) in enumerate(
        zip(normalized_points[:-1], normalized_points[1:])
    ):
        sample_count = start.segment_frames
        _validate_integer(sample_count, "segment_frames", minimum=2)
        key_quats = continuous_quaternions(
            np.vstack([start.rotation.as_quat(), stop.rotation.as_quat()])
        )
        slerp = Slerp([0.0, 1.0], Rotation.from_quat(key_quats))
        fractions = np.linspace(0.0, 1.0, sample_count)
        rotations = slerp(fractions)
        first_sample = 0 if segment_index == 0 else 1
        for sample_index in range(first_sample, sample_count):
            fraction = float(fractions[sample_index])
            at_start = segment_index == 0 and sample_index == 0
            at_stop = sample_index == sample_count - 1
            waypoint = start if at_start else stop if at_stop else None
            records.append(
                {
                    "segment_index": segment_index,
                    "from_label": start.label,
                    "to_label": stop.label,
                    "segment_fraction": fraction,
                    "rotation": rotations[sample_index],
                    "is_waypoint": waypoint is not None,
                    "waypoint_label": waypoint.label if waypoint else "",
                    "annotation": waypoint.annotation if waypoint else "",
                }
            )
            if at_stop:
                for _ in range(stop.hold_frames):
                    records.append(
                        {
                            "segment_index": segment_index,
                            "from_label": start.label,
                            "to_label": stop.label,
                            "segment_fraction": 1.0,
                            "rotation": rotations[sample_index],
                            "is_waypoint": True,
                            "waypoint_label": stop.label,
                            "annotation": stop.annotation,
                        }
                    )
    if normalized_points[0].hold_frames:
        first = records[0]
        records[1:1] = [dict(first) for _ in range(normalized_points[0].hold_frames)]
    if ping_pong:
        records.extend(_ping_pong_tail(records))
    return tuple(
        TrajectoryFrame(
            frame_index=index,
            time_sec=index / float(fps),
            segment_index=int(record["segment_index"]),
            from_label=str(record["from_label"]),
            to_label=str(record["to_label"]),
            segment_fraction=float(record["segment_fraction"]),
            rotation=record["rotation"],
            is_waypoint=bool(record["is_waypoint"]),
            waypoint_label=str(record["waypoint_label"]),
            annotation=str(record["annotation"]),
        )
        for index, record in enumerate(records)
    )


def write_trajectory_csv(
    frames: Iterable[TrajectoryFrame],
    path: str | Path,
    *,
    scipy_euler_sequence: str,
) -> Path:
    """Write public wxyz quaternions and representations from one Rotation."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRAJECTORY_COLUMNS)
        writer.writeheader()
        for frame in frames:
            xyzw = frame.rotation.as_quat()
            rv = frame.rotation.as_rotvec()
            ea = frame.rotation.as_euler(scipy_euler_sequence, degrees=True)
            writer.writerow(
                {
                    "frame_index": frame.frame_index,
                    "time_sec": frame.time_sec,
                    "segment_index": frame.segment_index,
                    "from_label": frame.from_label,
                    "to_label": frame.to_label,
                    "segment_fraction": frame.segment_fraction,
                    "quat_w": xyzw[3],
                    "quat_x": xyzw[0],
                    "quat_y": xyzw[1],
                    "quat_z": xyzw[2],
                    "rv_x_rad": rv[0],
                    "rv_y_rad": rv[1],
                    "rv_z_rad": rv[2],
                    "alpha_deg": ea[0],
                    "beta_deg": ea[1],
                    "gamma_deg": ea[2],
                    "is_waypoint": frame.is_waypoint,
                    "waypoint_label": frame.waypoint_label,
                    "annotation": frame.annotation,
                }
            )
    return output_path


def _ping_pong_tail(records: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    if len(records) < 2:
        return []
    tail: list[dict[str, object]] = []
    for record in reversed(records[:-1]):
        copied = dict(record)
        copied["from_label"], copied["to_label"] = (
            copied["to_label"],
            copied["from_label"],
        )
        copied["segment_fraction"] = 1.0 - float(copied["segment_fraction"])
        tail.append(copied)
    return tail


def _parse_float(value: object, column: str, row_number: int) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Waypoint row {row_number} column {column!r} must be numeric"
        ) from exc
    return parsed


def _optional_integer(
    value: object,
    default: int,
    name: str,
    row_number: int,
    *,
    minimum: int,
) -> int:
    if not _has_value(value):
        return default
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Waypoint row {row_number} {name} must be an integer") from exc
    if not np.isfinite(number) or not number.is_integer() or int(number) < minimum:
        comparator = "non-negative" if minimum == 0 else f">= {minimum}"
        raise ValueError(f"Waypoint row {row_number} {name} must be a {comparator} integer")
    return int(number)


def _validate_integer(value: int, name: str, *, minimum: int) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def _has_value(value: object) -> bool:
    return value is not None and str(value).strip() != ""
