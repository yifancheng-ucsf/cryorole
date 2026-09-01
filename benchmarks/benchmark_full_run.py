"""Full ``cryorole run`` benchmark for array-native and compatibility backends.

The normal benchmark is 100k rows. Use ``--rows 500000`` or ``--rows 1000000``
as explicit slow runs. Each backend runs in its own subprocess so peak RSS is
comparable and includes read, match, normalization, RO, SLD, NPZ, CSV, reports,
and transactional commit.
"""

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


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from cryorole.cli.main import build_parser, run_command


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Benchmark complete cryoROLE run bundles.")
    parser.add_argument("--rows", type=int, default=100_000)
    parser.add_argument("--k-neighbors", type=int, default=50)
    parser.add_argument("--query-batch-size", type=int, default=25_000)
    parser.add_argument("--raw-csv-chunk-size", type=int, default=25_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--backend",
        choices=("both", "array_native", "dataframe_compat"),
        default="both",
    )
    parser.add_argument("--keep-temp", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--ref", help=argparse.SUPPRESS)
    parser.add_argument("--mov", help=argparse.SUPPRESS)
    parser.add_argument("--output", help=argparse.SUPPRESS)
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.worker:
        print(json.dumps(_worker(args), sort_keys=True))
        return 0
    if args.rows < 2:
        raise ValueError("--rows must be at least 2")
    if min(args.k_neighbors, args.query_batch_size, args.raw_csv_chunk_size) < 1:
        raise ValueError("k, query batch, and CSV chunk sizes must be positive")
    root = Path(tempfile.mkdtemp(prefix="cryorole-full-run-benchmark-"))
    try:
        ref = root / "ref.cs"
        mov = root / "mov.cs"
        _write_inputs(ref, mov, rows=args.rows, seed=args.seed)
        backends = (
            ("array_native", "dataframe_compat")
            if args.backend == "both"
            else (args.backend,)
        )
        results = [
            _subprocess_result(args, backend, ref, mov, root / backend)
            for backend in backends
        ]
        print(
            json.dumps(
                {
                    "artifact_type": "cryorole_full_run_benchmark",
                    "rows": args.rows,
                    "k_neighbors": args.k_neighbors,
                    "query_batch_size": args.query_batch_size,
                    "raw_csv_chunk_size": args.raw_csv_chunk_size,
                    "memory_budget_mib": {
                        "array_native_100k_target": 1024,
                        "note": "Peak RSS is measured per subprocess; scale validation is empirical.",
                    },
                    "results": results,
                    "temporary_root": str(root) if args.keep_temp else None,
                },
                indent=2,
                sort_keys=True,
            )
        )
    finally:
        if not args.keep_temp:
            shutil.rmtree(root)
    return 0


def _write_inputs(ref: Path, mov: Path, *, rows: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    dtype = np.dtype([("uid", "<u8"), ("alignments3D/pose", "<f4", (3,))])
    reference = np.zeros(rows, dtype=dtype)
    moving = np.zeros(rows, dtype=dtype)
    reference["uid"] = np.arange(rows, dtype=np.uint64)
    moving["uid"] = reference["uid"]
    reference["alignments3D/pose"] = rng.normal(0.0, 0.2, size=(rows, 3))
    moving["alignments3D/pose"] = (
        reference["alignments3D/pose"] + rng.normal(0.0, 0.08, size=(rows, 3))
    )
    for path, values in ((ref, reference), (mov, moving)):
        with path.open("wb") as handle:
            np.save(handle, values)


def _subprocess_result(args, backend: str, ref: Path, mov: Path, output: Path) -> dict:
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--backend",
        backend,
        "--ref",
        str(ref),
        "--mov",
        str(mov),
        "--output",
        str(output),
        "--rows",
        str(args.rows),
        "--k-neighbors",
        str(args.k_neighbors),
        "--query-batch-size",
        str(args.query_batch_size),
        "--raw-csv-chunk-size",
        str(args.raw_csv_chunk_size),
    ]
    completed = subprocess.run(command, check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"full-run benchmark worker failed for {backend}:\n{completed.stderr}"
        )
    return json.loads(completed.stdout.splitlines()[-1])


def _worker(args) -> dict[str, object]:
    cli_args = build_parser().parse_args(
        [
            "run",
            "--ref",
            args.ref,
            "--mov",
            args.mov,
            "--output-dir",
            args.output,
            "--run-backend",
            args.backend,
            "--k-neighbors",
            str(args.k_neighbors),
            "--density-query-batch-size",
            str(args.query_batch_size),
            "--raw-csv-chunk-size",
            str(args.raw_csv_chunk_size),
            "--no-visualize",
            "--quiet",
        ]
    )
    start = perf_counter()
    run_command(cli_args)
    wall = perf_counter() - start
    output = Path(args.output)
    artifact_sizes = {
        str(path.relative_to(output)): path.stat().st_size
        for path in output.rglob("*")
        if path.is_file()
    }
    peak = _peak_rss_bytes()
    return {
        "backend": args.backend,
        "row_count": args.rows,
        "wall_time_sec": wall,
        "peak_rss_bytes": peak,
        "peak_rss_mib": peak / (1024 * 1024),
        "artifact_size_bytes": int(sum(artifact_sizes.values())),
        "artifact_sizes": artifact_sizes,
    }


def _peak_rss_bytes() -> int:
    if os.name != "nt":
        import resource

        peak = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        return peak if sys.platform == "darwin" else peak * 1024

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
    counters.cb = ctypes.sizeof(counters)
    get_process_memory_info = ctypes.windll.psapi.GetProcessMemoryInfo
    get_process_memory_info.argtypes = [
        wintypes.HANDLE,
        ctypes.POINTER(PROCESS_MEMORY_COUNTERS),
        wintypes.DWORD,
    ]
    get_process_memory_info.restype = wintypes.BOOL
    ok = get_process_memory_info(
        ctypes.windll.kernel32.GetCurrentProcess(),
        ctypes.byref(counters),
        counters.cb,
    )
    if not ok:
        raise RuntimeError("Could not read Windows peak RSS")
    return int(counters.PeakWorkingSetSize)


if __name__ == "__main__":
    raise SystemExit(main())
