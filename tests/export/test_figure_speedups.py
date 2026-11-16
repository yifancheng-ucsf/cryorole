"""Phase 2 figure speed-ups: parallel figure processes and the "density" point rendering."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from PIL import Image

from cryorole.cli.main import build_parser
from cryorole.export.figure_jobs import PARALLEL_FIGURES_MIN_POINTS, resolve_figure_jobs

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run"))
from test_core_productionization import _write_cs  # noqa: E402

N = 3000


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    rng = np.random.default_rng(3)
    ref_rv = rng.normal(0.0, 0.5, (N, 3))
    mov_rv = ref_rv + rng.normal(0.0, 0.1, (N, 3))
    uids = list(range(1, N + 1))
    _write_cs(tmp_path / "ref.cs", uids, ref_rv)
    _write_cs(tmp_path / "mov.cs", uids, mov_rv)
    return tmp_path / "ref.cs", tmp_path / "mov.cs"


def _run(argv: list[str]) -> int:
    args = build_parser().parse_args(argv)
    return args.handler(args)


def _pixels(path: Path) -> np.ndarray:
    return np.asarray(Image.open(path).convert("RGB"), dtype=int)


def _same_pngs(first: Path, second: Path) -> int:
    pngs = sorted(first.rglob("*.png"))
    assert pngs
    for png in pngs:
        other = second / png.relative_to(first)
        assert np.array_equal(_pixels(png), _pixels(other)), png.name
    return len(pngs)


def test_auto_jobs_stay_serial_for_small_inputs_and_explicit_jobs_are_capped() -> None:
    assert resolve_figure_jobs("auto", n_tasks=3, n_points=PARALLEL_FIGURES_MIN_POINTS - 1, worker_bytes_per_point=1) == 1
    assert resolve_figure_jobs(None, n_tasks=3, n_points=10, worker_bytes_per_point=1) == 1
    assert resolve_figure_jobs(8, n_tasks=3, n_points=10, worker_bytes_per_point=1) == 3
    assert resolve_figure_jobs(2, n_tasks=3, n_points=10, worker_bytes_per_point=1, parallel_capable=False) == 1
    with pytest.raises(ValueError):
        resolve_figure_jobs(0, n_tasks=3, n_points=10, worker_bytes_per_point=1)


def test_jobs_option_rejects_non_positive_values() -> None:
    with pytest.raises(SystemExit):
        build_parser().parse_args(["canonicalize", "--jobs", "0"])
    assert build_parser().parse_args(["run", "--ref", "a", "--mov", "b", "--jobs", "auto"]).jobs == "auto"


def test_parallel_figures_are_pixel_identical_to_serial(tmp_path) -> None:
    ref, mov = _inputs(tmp_path)
    for jobs in ("1", "3"):
        assert _run(["run", "--ref", str(ref), "--mov", str(mov), "--output-dir", str(tmp_path / f"run{jobs}"),
                     "--jobs", jobs, "--quiet"]) == 0
        assert _run(["canonicalize", "--run-dir", str(tmp_path / f"run{jobs}"), "--jobs", jobs]) == 0

    assert _same_pngs(tmp_path / "run1" / "visualizations", tmp_path / "run3" / "visualizations") == 7 + 12
    summary = json.loads((tmp_path / "run3" / "canonical" / "default" / "canonicalize_summary.json").read_text())
    assert summary["preview_jobs"] == 3
    serial = json.loads((tmp_path / "run1" / "canonical" / "default" / "canonicalize_summary.json").read_text())
    assert serial["preview_jobs"] == 1


def test_density_point_style_matches_the_scatter_picture(tmp_path) -> None:
    ref, mov = _inputs(tmp_path)
    for style in ("scatter", "density"):
        assert _run(["run", "--ref", str(ref), "--mov", str(mov), "--output-dir", str(tmp_path / style),
                     "--point-style", style, "--jobs", "1", "--quiet"]) == 0

    for name in ("all_rotvec_3view_projection.png", "all_euler_3view_projection.png"):
        scatter = _pixels(tmp_path / "scatter" / "visualizations" / "quicklook" / name)
        density = _pixels(tmp_path / "density" / "visualizations" / "quicklook" / name)
        assert scatter.shape == density.shape  # same layout, axes and colour bar
        coloured_scatter = scatter.min(-1) <= 245
        coloured_density = density.min(-1) <= 245
        # Every density cell is where the scatter also drew; the scatter's markers are
        # ~1.4 px wide, so it colours some extra edge pixels.
        assert (coloured_density & ~coloured_scatter).sum() <= 0.01 * coloured_density.sum()
        assert coloured_density.sum() >= 0.6 * coloured_scatter.sum()
        both = coloured_scatter & coloured_density
        agree = np.abs(scatter - density).max(-1)[both] <= 32
        assert agree.mean() >= 0.95
