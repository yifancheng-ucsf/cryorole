"""Typed records shared by the offline animation stages."""

from __future__ import annotations

from dataclasses import dataclass

from scipy.spatial.transform import Rotation


@dataclass(frozen=True)
class Waypoint:
    """One validated SO(3) waypoint in the declared waypoint coordinate set."""

    label: str
    rotation: Rotation
    segment_frames: int
    hold_frames: int
    annotation: str = ""


@dataclass(frozen=True)
class TrajectoryFrame:
    """One deterministic frame produced by SO(3) interpolation."""

    frame_index: int
    time_sec: float
    segment_index: int
    from_label: str
    to_label: str
    segment_fraction: float
    rotation: Rotation
    is_waypoint: bool
    waypoint_label: str
    annotation: str = ""
