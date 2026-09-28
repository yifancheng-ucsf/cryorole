"""Benchmark ``cryorole align`` on production-size RELION STAR files.

Target (align plan v3, A1): two 1.2M-particle, 23-column RELION 3.1 STAR files
(~150 MB each, as in job043/job055) aligned in ≤ 60 s wall time with ≤ 2 GB
peak RSS. The report records the machine (CPU, cores, RAM, OS, disk of the
work directory) next to the timings, because the target is hardware-specific.

Example:
    python benchmarks/benchmark_align.py --rows 1200000 --workdir /path/on/local/disk
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

COLUMNS = [
    "_rlnImageName", "_rlnMicrographName", "_rlnCoordinateX", "_rlnCoordinateY", "_rlnAngleRot", "_rlnAngleTilt",
    "_rlnAnglePsi", "_rlnOriginXAngst", "_rlnOriginYAngst", "_rlnDefocusU", "_rlnDefocusV", "_rlnDefocusAngle",
    "_rlnPhaseShift", "_rlnCtfBfactor", "_rlnOpticsGroup", "_rlnClassNumber", "_rlnNormCorrection",
    "_rlnRandomSubset", "_rlnLogLikeliContribution", "_rlnMaxValueProbDistribution", "_rlnNrOfSignificantSamples",
    "_rlnImageOriginalName", "_rlnGroupNumber",
]


def _write_star(path: Path, rows: int, *, image_prefix: str, original_prefix: str, order: np.ndarray, rng) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\n# version 30001\n\ndata_optics\n\nloop_\n_rlnOpticsGroupName #1\n_rlnOpticsGroup #2\n"
                     "_rlnMicrographOriginalPixelSize #3\n_rlnImagePixelSize #4\n_rlnImageSize #5\n"
                     "_rlnImageDimensionality #6\nopticsGroup1 1 0.835 0.835 448 2\n\n\n# version 30001\n\n"
                     "data_particles\n\nloop_\n")
        for index, column in enumerate(COLUMNS, start=1):
            handle.write(f"{column} #{index}\n")
        chunk = 100_000
        for start in range(0, rows, chunk):
            ids = order[start:start + chunk]
            n = len(ids)
            angles = rng.uniform(-180, 180, (n, 3))
            lines = []
            for k, particle in enumerate(ids):
                mic = particle // 8
                lines.append(
                    f"{particle + 1:06d}@{image_prefix}/mic_{mic:06d}.mrcs mics/mic_{mic:06d}_DW.mrc "
                    f"{1000 + particle % 3000:.6f} {2000 + particle % 2000:.6f} {angles[k, 0]:.6f} "
                    f"{abs(angles[k, 1]):.6f} {angles[k, 2]:.6f} 1.234567 -2.345678 18540.552734 17598.398438 "
                    f"-52.03189 0.000000 0.000000 1 1 0.516166 {1 + particle % 2} 8.816460e+05 0.443101 5 "
                    f"{particle + 1:06d}@{original_prefix}/mic_{mic:06d}.mrcs 1\n"
                )
            handle.write("".join(lines))


def _machine() -> dict:
    info = {
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "processor": platform.processor() or platform.machine(),
        "cpu_count": os.cpu_count(),
    }
    try:
        info["ram_gib"] = round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30, 1)
    except (ValueError, OSError, AttributeError):
        info["ram_gib"] = None
    if sys.platform == "darwin":
        try:
            info["cpu_model"] = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True,
                                               text=True, check=False).stdout.strip()
        except OSError:
            pass
    elif Path("/proc/cpuinfo").is_file():
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                info["cpu_model"] = line.split(":", 1)[1].strip()
                break
    return info


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=int, default=1_200_000)
    parser.add_argument("--workdir", help="Directory for the synthetic STAR files (use a local disk).")
    parser.add_argument("--report", help="Write the JSON report here.")
    args = parser.parse_args()

    from cryorole.align.star_align import align_star_files

    workdir = Path(args.workdir) if args.workdir else Path(tempfile.mkdtemp(prefix="cryorole_align_bench_"))
    workdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    ref, mov = workdir / "consensus.star", workdir / "body.star"
    started = time.perf_counter()
    _write_star(ref, args.rows, image_prefix="Extract/job018", original_prefix="Import/job001",
                order=np.arange(args.rows), rng=rng)
    _write_star(mov, args.rows, image_prefix="Subtract/job041", original_prefix="Extract/job018",
                order=rng.permutation(args.rows), rng=rng)
    generation_s = time.perf_counter() - started

    started = time.perf_counter()
    report = align_star_files(ref=ref, mov=mov, key_pairs=["_rlnImageName=_rlnImageOriginalName"],
                              output_dir=workdir / "out", overwrite=True)
    wall_s = time.perf_counter() - started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_gib = peak / 2**30 if sys.platform == "darwin" else peak / 2**20  # bytes on macOS, KiB on Linux
    result = {
        "artifact_type": "cryorole_align_benchmark",
        "rows": args.rows,
        "columns": len(COLUMNS),
        "ref_bytes": ref.stat().st_size,
        "mov_bytes": mov.stat().st_size,
        "workdir": str(workdir),
        "generation_seconds": round(generation_s, 2),
        "align_wall_seconds": round(wall_s, 2),
        "peak_rss_gib": round(peak_gib, 3),
        "matched": report["matched_count"],
        "target": {"wall_seconds": 60, "peak_rss_gib": 2.0},
        "meets_target": wall_s <= 60 and peak_gib <= 2.0,
        "machine": _machine(),
    }
    text = json.dumps(result, indent=2)
    print(text)
    if args.report:
        Path(args.report).write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
