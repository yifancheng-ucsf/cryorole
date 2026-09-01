"""100k preflight, visualization, and selection benchmark for Workflow UX."""

from __future__ import annotations

import argparse
import ctypes
from ctypes import wintypes
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from time import perf_counter

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryorole.interactive import ExploreSession  # noqa: E402
from cryorole.preflight import PreflightRequest, run_preflight  # noqa: E402
from cryorole.select import SelectRequest, create_selection  # noqa: E402
from cryorole.visualize import VisualizationRequest, visualize  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=int, default=100_000)
    parser.add_argument("--max-display-points", type=int, default=50_000)
    parser.add_argument(
        "--worker",
        choices=("preflight", "interactive", "visualize", "select"),
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--root", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(_worker(args), sort_keys=True))
        return 0
    temporary = Path(tempfile.mkdtemp(prefix="cryorole-workflow-ux-benchmark-"))
    try:
        _write_inputs_and_bundle(temporary, args.rows)
        results = [
            _run_worker(temporary, args, "preflight"),
            _run_worker(temporary, args, "interactive"),
            _run_worker(temporary, args, "visualize"),
            _run_worker(temporary, args, "select"),
        ]
        print(json.dumps({
            "artifact_type": "cryorole_workflow_ux_benchmark",
            "rows": args.rows,
            "max_display_points": args.max_display_points,
            "results": results,
            "notes": [
                "Workers run in separate subprocesses.",
                "Interactive worker loads NPZ once; no raw CSV is present.",
                "Radius evaluation is exact on the full parent landscape.",
                "Visualization filters compact NPZ arrays before plotting-table materialization.",
                "Scientific select evaluates the full parent NPZ, never the display sample.",
            ],
        }, indent=2, sort_keys=True))
    finally:
        shutil.rmtree(temporary)
    return 0


def _write_inputs_and_bundle(root: Path, rows: int) -> None:
    rng = np.random.default_rng(4)
    dtype = np.dtype([("uid", "<u8"), ("alignments3D/pose", "<f4", (3,))])
    ref = np.zeros(rows, dtype=dtype)
    mov = np.zeros(rows, dtype=dtype)
    ref["uid"] = np.arange(rows, dtype=np.uint64)
    mov["uid"] = ref["uid"]
    ref["alignments3D/pose"] = rng.normal(0, 0.2, size=(rows, 3))
    mov["alignments3D/pose"] = ref["alignments3D/pose"] + rng.normal(0, 0.05, size=(rows, 3))
    for name, values in (("ref.cs", ref), ("mov.cs", mov)):
        with (root / name).open("wb") as handle:
            np.save(handle, values)
    run = root / "run"
    (run / "data").mkdir(parents=True)
    coordinates = rng.normal(0, 0.25, size=(rows, 3))
    np.savez_compressed(
        run / "data" / "raw_landscape.npz",
        artifact_type=np.asarray("raw_landscape"), schema_version=np.asarray("1"),
        particle_key=np.asarray([str(index) for index in range(rows)]),
        coordinates_analysis=coordinates, coordinates_display=coordinates,
        sld_unfloored=np.ones(rows), sld_raw=np.ones(rows), sld_display=np.ones(rows),
        sld_display_is_outlier=np.zeros(rows, bool), sld_was_floored=np.zeros(rows, bool),
        sld_local_k_mean=np.ones(rows), sld_effective_local_k_mean=np.ones(rows),
        sld_distance_floor=np.full(rows, 1e-4), ref_source_row_id=np.arange(rows),
        mov_source_row_id=np.arange(rows),
    )
    (run / "run_manifest.json").write_text(
        json.dumps({"run_id": "benchmark", "bundle_transaction": {"state": "completed"}}), encoding="utf-8"
    )
    (run / "run_summary.json").write_text(
        json.dumps({"run_id": "benchmark", "euler_convention": "extrinsic_zyx", "scipy_euler_sequence": "zyx"}), encoding="utf-8"
    )
    (run / ".cryorole_bundle_complete").write_text("{}", encoding="utf-8")


def _run_worker(root: Path, args, worker: str) -> dict[str, object]:
    completed = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), "--worker", worker, "--root", str(root),
         "--rows", str(args.rows), "--max-display-points", str(args.max_display_points)],
        capture_output=True, text=True, check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"{worker} benchmark worker failed:\n{completed.stderr}")
    return json.loads(completed.stdout.splitlines()[-1])


def _worker(args) -> dict[str, object]:
    root = Path(args.root)
    if args.worker == "preflight":
        started = perf_counter()
        result = run_preflight(PreflightRequest(ref=root / "ref.cs", mov=root / "mov.cs", visualize=False))
        elapsed = perf_counter() - started
        return {
            "stage": "preflight", "readiness": result.report["readiness"],
            "matched_count": result.report["matching"]["matched_count"],
            "wall_time_sec": elapsed, "peak_rss_mib": _peak_rss_bytes() / (1024 * 1024),
        }
    if args.worker == "visualize":
        started = perf_counter()
        result = visualize(
            VisualizationRequest(
                run_dir=str(root / "run"),
                visual_id="benchmark",
                representation="rotvec",
                formats="png",
                overwrite=True,
            )
        )
        elapsed = perf_counter() - started
        return {
            "stage": "visualize",
            "wall_time_sec": elapsed,
            "peak_rss_mib": _peak_rss_bytes() / (1024 * 1024),
            "full_parent_count": result.report.get("n_points_full_parent"),
            "input_count": result.report["n_points_input"],
            "filtered_count": result.report["n_points_after_display_filter"],
            "count_2d": result.report["n_points_2d"],
            "count_3d": result.report["n_points_3d"],
            "artifact_size_bytes": sum(
                path.stat().st_size for path in result.output_dir.rglob("*") if path.is_file()
            ),
        }
    if args.worker == "select":
        started = perf_counter()
        result = create_selection(
            SelectRequest(
                run_dir=str(root / "run"),
                selection_id="benchmark_radius",
                center=(0.0, 0.0, 0.0),
                center_representation="rotvec",
                radius_rad=0.1,
                overwrite=True,
            )
        )
        elapsed = perf_counter() - started
        selection_payload = json.loads(
            (result.output_dir / "selection.json").read_text(encoding="utf-8")
        )
        return {
            "stage": "select",
            "wall_time_sec": elapsed,
            "peak_rss_mib": _peak_rss_bytes() / (1024 * 1024),
            "candidate_count": selection_payload["total_count"],
            "selected_count": selection_payload["selected_count"],
            "artifact_size_bytes": sum(
                path.stat().st_size for path in result.output_dir.rglob("*") if path.is_file()
            ),
        }
    started = perf_counter()
    session = ExploreSession(root / "run", max_display_points=args.max_display_points)
    load_time = perf_counter() - started
    evaluation_started = perf_counter()
    draft = session.evaluate(center=[0, 0, 0], representation="rotvec", radius_deg=6)
    evaluation_time = perf_counter() - evaluation_started
    return {
        "stage": "interactive", "row_count": session.arrays.n_points,
        "displayed_point_count": len(session.display_indices),
        "initial_load_time_sec": load_time, "full_radius_evaluation_time_sec": evaluation_time,
        "selected_count": draft["selected_count"], "exact": draft["exact"],
        "raw_csv_present": (root / "run" / "data" / "raw_landscape.csv").exists(),
        "peak_rss_mib": _peak_rss_bytes() / (1024 * 1024),
    }


def _peak_rss_bytes() -> int:
    if os.name != "nt":
        import resource
        value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        return value if sys.platform == "darwin" else value * 1024
    class Counters(ctypes.Structure):
        _fields_ = [
            ("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
            ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t),
        ]
    counters = Counters()
    counters.cb = ctypes.sizeof(counters)
    get_process_memory_info = ctypes.windll.psapi.GetProcessMemoryInfo
    get_process_memory_info.argtypes = [
        wintypes.HANDLE,
        ctypes.POINTER(Counters),
        wintypes.DWORD,
    ]
    get_process_memory_info.restype = wintypes.BOOL
    ok = get_process_memory_info(
        ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb
    )
    if not ok:
        raise RuntimeError("Could not read Windows peak RSS")
    return int(counters.PeakWorkingSetSize)


if __name__ == "__main__":
    raise SystemExit(main())
