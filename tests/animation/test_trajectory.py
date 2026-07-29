from __future__ import annotations

import csv

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole.animation.trajectory import (
    build_trajectory,
    parse_waypoint_csv,
)


def _write_csv(path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_parse_ea_waypoints_validates_schema_and_values(tmp_path):
    path = tmp_path / "path.csv"
    _write_csv(
        path,
        ["label", "alpha_deg", "beta_deg", "gamma_deg"],
        [
            {"label": "a", "alpha_deg": 0, "beta_deg": 0, "gamma_deg": 0},
            {"label": "b", "alpha_deg": 10, "beta_deg": 20, "gamma_deg": 30},
        ],
    )
    waypoints, warnings = parse_waypoint_csv(
        path,
        path_space="ea",
        scipy_euler_sequence="zyx",
        default_segment_frames=5,
        default_hold_frames=0,
    )
    assert [point.label for point in waypoints] == ["a", "b"]
    assert warnings == ()


@pytest.mark.parametrize(
    "rows,match",
    [
        (
            [
                {"label": "a", "rv_x_rad": 0, "rv_y_rad": 0, "rv_z_rad": 0},
                {"label": "a", "rv_x_rad": 0, "rv_y_rad": 0, "rv_z_rad": 1},
            ],
            "unique",
        ),
        (
            [
                {"label": "a", "rv_x_rad": 0, "rv_y_rad": 0, "rv_z_rad": 0},
                {"label": "b", "rv_x_rad": "nan", "rv_y_rad": 0, "rv_z_rad": 1},
            ],
            "finite",
        ),
    ],
)
def test_parse_waypoints_rejects_duplicate_and_nonfinite(tmp_path, rows, match):
    path = tmp_path / "path.csv"
    _write_csv(path, ["label", "rv_x_rad", "rv_y_rad", "rv_z_rad"], rows)
    with pytest.raises(ValueError, match=match):
        parse_waypoint_csv(
            path,
            path_space="rv",
            scipy_euler_sequence="zyx",
            default_segment_frames=5,
            default_hold_frames=0,
        )


def test_parse_waypoints_rejects_mixed_schema(tmp_path):
    path = tmp_path / "path.csv"
    _write_csv(
        path,
        [
            "label",
            "alpha_deg",
            "beta_deg",
            "gamma_deg",
            "rv_x_rad",
            "rv_y_rad",
            "rv_z_rad",
        ],
        [
            {
                "label": "a",
                "alpha_deg": 0,
                "beta_deg": 0,
                "gamma_deg": 0,
                "rv_x_rad": 0,
                "rv_y_rad": 0,
                "rv_z_rad": 0,
            },
            {
                "label": "b",
                "alpha_deg": 1,
                "beta_deg": 2,
                "gamma_deg": 3,
                "rv_x_rad": 0,
                "rv_y_rad": 0,
                "rv_z_rad": 1,
            },
        ],
    )
    with pytest.raises(ValueError, match="mix"):
        parse_waypoint_csv(
            path,
            path_space="ea",
            scipy_euler_sequence="zyx",
            default_segment_frames=5,
            default_hold_frames=0,
        )


def test_slerp_multisegment_holds_reverse_and_ping_pong(tmp_path):
    path = tmp_path / "path.csv"
    _write_csv(
        path,
        ["label", "rv_x_rad", "rv_y_rad", "rv_z_rad", "segment_frames", "hold_frames"],
        [
            {
                "label": "a",
                "rv_x_rad": 0,
                "rv_y_rad": 0,
                "rv_z_rad": 0,
                "segment_frames": 3,
                "hold_frames": 0,
            },
            {
                "label": "b",
                "rv_x_rad": 0,
                "rv_y_rad": 0,
                "rv_z_rad": np.pi / 2,
                "segment_frames": 4,
                "hold_frames": 2,
            },
            {
                "label": "c",
                "rv_x_rad": 0,
                "rv_y_rad": np.pi / 2,
                "rv_z_rad": 0,
                "segment_frames": 7,
                "hold_frames": 0,
            },
        ],
    )
    waypoints, warnings = parse_waypoint_csv(
        path,
        path_space="rv",
        scipy_euler_sequence="zyx",
        default_segment_frames=7,
        default_hold_frames=0,
    )
    assert warnings
    forward = build_trajectory(waypoints, fps=10)
    assert len(forward) == 3 + 4 - 1 + 2
    assert np.allclose(forward[0].rotation.as_matrix(), waypoints[0].rotation.as_matrix())
    assert np.allclose(forward[-1].rotation.as_matrix(), waypoints[-1].rotation.as_matrix())
    assert sum(frame.waypoint_label == "b" for frame in forward) == 3
    assert all(np.isclose(np.linalg.norm(frame.rotation.as_quat()), 1.0) for frame in forward)

    reverse = build_trajectory(tuple(reversed(waypoints)), fps=10)
    assert np.allclose(reverse[0].rotation.as_matrix(), waypoints[-1].rotation.as_matrix())
    ping_pong = build_trajectory(waypoints, fps=10, ping_pong=True)
    assert np.allclose(ping_pong[-1].rotation.as_matrix(), waypoints[0].rotation.as_matrix())
    assert not np.allclose(
        ping_pong[len(forward) - 1].rotation.as_matrix(),
        ping_pong[len(forward)].rotation.as_matrix(),
    )
    assert [frame.frame_index for frame in ping_pong] == list(range(len(ping_pong)))


def test_quaternion_sign_continuity():
    from cryorole.animation.trajectory import continuous_quaternions

    values = np.asarray(
        [
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, -0.5, -np.sqrt(0.75)],
        ]
    )
    continuous = continuous_quaternions(values)
    assert np.dot(continuous[0], continuous[1]) >= 0


def test_ea_rv_round_trip():
    rotation = Rotation.from_euler("zyx", [25.0, -15.0, 40.0], degrees=True)
    rebuilt = Rotation.from_rotvec(rotation.as_rotvec())
    assert np.allclose(rotation.as_matrix(), rebuilt.as_matrix(), atol=1e-12)
