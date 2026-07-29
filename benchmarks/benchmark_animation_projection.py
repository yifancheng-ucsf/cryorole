"""Manual production-scale benchmark for static animation backgrounds.

Example:
    python benchmarks/benchmark_animation_projection.py --rows 1000000 --frames 10
"""

from __future__ import annotations

import argparse
import tempfile
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from cryorole.animation.projection import ProjectionRenderer, resolve_axis_limits
from cryorole.animation.schemas import TrajectoryFrame


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=int, default=1_000_000)
    parser.add_argument("--frames", type=int, default=10)
    args = parser.parse_args()
    rng = np.random.default_rng(0)
    landscape = rng.normal(size=(args.rows, 3)).astype(np.float32) * 25
    density = rng.lognormal(size=args.rows).astype(np.float32)
    rotations = Rotation.from_euler(
        "zyx",
        np.column_stack(
            [
                np.linspace(-15, 15, args.frames),
                np.linspace(5, 20, args.frames),
                np.linspace(-10, 25, args.frames),
            ]
        ),
        degrees=True,
    )
    frames = tuple(
        TrajectoryFrame(
            frame_index=index,
            time_sec=index / 30,
            segment_index=0,
            from_label="start",
            to_label="stop",
            segment_fraction=index / max(args.frames - 1, 1),
            rotation=rotation,
            is_waypoint=index in {0, args.frames - 1},
            waypoint_label="start" if index == 0 else "stop" if index == args.frames - 1 else "",
        )
        for index, rotation in enumerate(rotations)
    )
    trajectory = np.vstack([frame.rotation.as_euler("zyx", degrees=True) for frame in frames])
    limits, _ = resolve_axis_limits(
        landscape,
        trajectory,
        explicit=None,
        margin_degrees=2,
    )
    started = time.perf_counter()
    renderer = ProjectionRenderer(
        landscape_euler=landscape,
        sld_display=density,
        trajectory_frames=frames,
        scipy_euler_sequence="zyx",
        axis_limits=limits,
        canvas_width=1800,
        canvas_height=600,
        point_size=0.2,
        landscape_alpha=0.4,
        trail_frames=5,
        show_current_values=False,
    )
    with tempfile.TemporaryDirectory(prefix="cryorole_animation_benchmark_") as directory:
        renderer.render(Path(directory))
    elapsed = time.perf_counter() - started
    if renderer.static_scatter_count != 3:
        raise RuntimeError("Renderer rebuilt static landscape scatter artists")
    if renderer.static_background_draw_count != 1:
        raise RuntimeError("Renderer redrew the static landscape background")
    print(
        f"rows={args.rows} frames={args.frames} "
        f"static_scatters={renderer.static_scatter_count} "
        f"static_background_draws={renderer.static_background_draw_count} "
        f"wall_seconds={elapsed:.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
