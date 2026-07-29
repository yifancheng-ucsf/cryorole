"""Compare production SLD metrics outside the normal pytest suite.

Run from the repository root:

    python benchmarks/benchmark_sld_metrics.py
"""

from __future__ import annotations

import argparse
import ctypes
from ctypes import wintypes
import json
import os
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import numpy as np
from scipy.spatial.transform import Rotation


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from cryorole.core.density import SLD_METRICS, compute_sld_values


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Benchmark cryoROLE SLD metrics.")
    parser.add_argument("--rows", type=int, default=100_000)
    parser.add_argument("--k-neighbors", type=int, default=50)
    parser.add_argument("--query-batch-size", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--metric", choices=("both", *SLD_METRICS), default="both")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.rows < 2:
        raise ValueError("--rows must be at least 2")
    if args.k_neighbors < 1:
        raise ValueError("--k-neighbors must be at least 1")
    if args.query_batch_size < 1:
        raise ValueError("--query-batch-size must be at least 1")

    if args.worker:
        print(json.dumps(_run_worker(args), sort_keys=True))
        return 0

    metrics = SLD_METRICS if args.metric == "both" else (args.metric,)
    results = [_run_metric_subprocess(args, metric) for metric in metrics]
    print(
        json.dumps(
            {
                "artifact_type": "sld_metric_benchmark",
                "rows": args.rows,
                "k_neighbors": args.k_neighbors,
                "query_batch_size": args.query_batch_size,
                "seed": args.seed,
                "results": results,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


def _run_metric_subprocess(args, metric: str) -> dict[str, object]:
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--metric",
        metric,
        "--rows",
        str(args.rows),
        "--k-neighbors",
        str(args.k_neighbors),
        "--query-batch-size",
        str(args.query_batch_size),
        "--seed",
        str(args.seed),
    ]
    completed = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(completed.stdout)


def _run_worker(args) -> dict[str, object]:
    rng = np.random.default_rng(args.seed)
    coordinates = Rotation.random(args.rows, random_state=rng).as_rotvec()
    start = perf_counter()
    result = compute_sld_values(
        coordinates,
        k_neighbors=args.k_neighbors,
        distance_floor_fraction=1e-4,
        sld_metric=args.metric,
        query_batch_size=args.query_batch_size,
    )
    wall_time_sec = perf_counter() - start
    sld_raw = np.asarray(result["sld_raw"], dtype=float)
    peak_rss_bytes = _peak_rss_bytes()
    return {
        "metric": args.metric,
        "row_count": args.rows,
        "wall_time_sec": wall_time_sec,
        "peak_rss_bytes": peak_rss_bytes,
        "peak_rss_mib": peak_rss_bytes / (1024 * 1024),
        "finite_sld_count": int(np.isfinite(sld_raw).sum()),
        "sld_checksum": float(np.sum(sld_raw)),
    }


def _peak_rss_bytes() -> int:
    if os.name == "nt":
        return _peak_rss_bytes_windows()
    import resource

    peak = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return peak if sys.platform == "darwin" else peak * 1024


def _peak_rss_bytes_windows() -> int:
    class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
        _fields_ = [
            ("cb", wintypes.DWORD),
            ("PageFaultCount", wintypes.DWORD),
            ("PeakWorkingSetSize", ctypes.c_size_t),
            ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t),
            ("PeakPagefileUsage", ctypes.c_size_t),
        ]

    counters = PROCESS_MEMORY_COUNTERS()
    counters.cb = ctypes.sizeof(PROCESS_MEMORY_COUNTERS)
    get_process_memory_info = ctypes.windll.psapi.GetProcessMemoryInfo
    get_process_memory_info.argtypes = [
        wintypes.HANDLE,
        ctypes.POINTER(PROCESS_MEMORY_COUNTERS),
        wintypes.DWORD,
    ]
    get_process_memory_info.restype = wintypes.BOOL
    process = ctypes.windll.kernel32.GetCurrentProcess()
    ok = get_process_memory_info(
        process,
        ctypes.byref(counters),
        counters.cb,
    )
    if not ok:
        raise RuntimeError("Could not read peak RSS from Windows Process Status API")
    return int(counters.PeakWorkingSetSize)


if __name__ == "__main__":
    raise SystemExit(main())
