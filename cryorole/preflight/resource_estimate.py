"""Transparent deterministic resource estimates for production run preflight."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
import shutil


MIB = 1024 * 1024

# Peak-RSS model, fitted on 2026-09-28 to measured `run` peaks (os.wait4 ru_maxrss)
# from 20k to 1.16M particles: synthetic CryoSPARC inputs, the J75/J80 CryoSPARC
# pair and a RELION job018/job022 pair (Linux, NumPy 2.2, Matplotlib 3.10). Every
# measured case lies within -5 % ... +15 % of the model (it mostly errs high).
FIXED_RUNTIME_BYTES = 110 * MIB
COMPACT_ARRAY_BYTES_PER_PARTICLE = 740
# .cs inputs are memory-mapped; the touched pages count toward RSS.
CS_INPUT_RESIDENT_FRACTION = 0.75
# .star inputs are parsed from text; the parse holds about twice the file size.
STAR_INPUT_RESIDENT_FACTOR = 2.0
# Raw quick-look figures (run_service.RUN_QUICKLOOK_MAX_POINTS caps their points).
PREVIEW_FIXED_BYTES = 60 * MIB
PREVIEW_BYTES_PER_PLOTTED_PARTICLE = 1300
PREVIEW_MAX_PLOTTED_PARTICLES = 500_000


def estimate_run_resources(
    *,
    matched_count: int,
    k_neighbors: int,
    query_batch_size: int,
    raw_csv: bool,
    visualize: bool,
    output_dir: str | Path,
    input_paths: Sequence[str | Path] = (),
) -> dict[str, object]:
    """Return a conservative estimate with every formula input exposed."""

    n = max(int(matched_count), 0)
    k = max(int(k_neighbors), 1)
    batch = max(int(query_batch_size), 1)
    bounded_batch = min(batch, max(1, 2_000_000 // (k + 1)))
    fixed_bytes = FIXED_RUNTIME_BYTES
    compact_arrays_bytes = n * COMPACT_ARRAY_BYTES_PER_PARTICLE
    knn_query_bytes = min(bounded_batch, n) * (k + 1) * 16
    cs_input_bytes = 0
    star_input_bytes = 0
    for input_path in input_paths:
        path = Path(input_path)
        try:
            size = path.stat().st_size
        except OSError:
            continue
        if path.suffix.lower() == ".cs":
            cs_input_bytes += size
        else:
            star_input_bytes += size
    input_resident_bytes = int(
        cs_input_bytes * CS_INPUT_RESIDENT_FRACTION + star_input_bytes * STAR_INPUT_RESIDENT_FACTOR
    )
    plotted = min(n, PREVIEW_MAX_PLOTTED_PARTICLES) if visualize else 0
    preview_memory_bytes = (PREVIEW_FIXED_BYTES + plotted * PREVIEW_BYTES_PER_PLOTTED_PARTICLE) if visualize else 0
    estimated_peak = (
        fixed_bytes + compact_arrays_bytes + knn_query_bytes + input_resident_bytes + preview_memory_bytes
    )
    npz_bytes = n * 96 + 64 * 1024
    csv_bytes = n * 420 if raw_csv else 0
    preview_rows = min(n, 50_000) if visualize else 0
    preview_bytes = (8 * MIB + preview_rows * 280) if visualize else 0
    report_bytes = 2 * MIB
    disk_bytes = npz_bytes + csv_bytes + preview_bytes + report_bytes

    target = Path(output_dir).expanduser().resolve()
    probe = target.parent
    while not probe.exists() and probe.parent != probe:
        probe = probe.parent
    free_bytes = int(shutil.disk_usage(probe).free)
    disk_sufficient = free_bytes >= int(disk_bytes * 1.2)
    recommend_no_visualize = bool(n >= 500_000 or preview_bytes >= 128 * MIB)
    return {
        "estimated_peak_memory_bytes": int(estimated_peak),
        "estimated_peak_memory_mib": round(estimated_peak / MIB, 2),
        "estimated_bundle_bytes": int(disk_bytes),
        "estimated_bundle_mib": round(disk_bytes / MIB, 2),
        "available_disk_bytes": free_bytes,
        "disk_sufficient_with_20pct_margin": disk_sufficient,
        "recommend_no_visualize": recommend_no_visualize,
        "assumptions": {
            "matched_count": n,
            "k_neighbors": k,
            "requested_query_batch_size": batch,
            "memory_bounded_query_batch_size": bounded_batch,
            "fixed_runtime_bytes": fixed_bytes,
            "compact_array_bytes_per_particle": COMPACT_ARRAY_BYTES_PER_PARTICLE,
            "knn_distance_and_index_bytes_per_element": 16,
            "cs_input_bytes": cs_input_bytes,
            "cs_input_resident_fraction": CS_INPUT_RESIDENT_FRACTION,
            "star_input_bytes": star_input_bytes,
            "star_input_resident_factor": STAR_INPUT_RESIDENT_FACTOR,
            "input_resident_bytes": input_resident_bytes,
            "preview_memory_fixed_bytes": PREVIEW_FIXED_BYTES if visualize else 0,
            "preview_memory_bytes_per_plotted_particle": PREVIEW_BYTES_PER_PLOTTED_PARTICLE if visualize else 0,
            "preview_max_plotted_particles": PREVIEW_MAX_PLOTTED_PARTICLES if visualize else 0,
            "preview_memory_bytes": preview_memory_bytes,
            "model_fit": "measured run peak RSS, 20k-1.16M particles, CryoSPARC and RELION inputs, 2026-09-28",
            "npz_bytes_per_particle": 96,
            "raw_csv_bytes_per_particle": 420 if raw_csv else 0,
            "preview_max_rows": 50_000 if visualize else 0,
            "preview_fixed_bytes": 8 * MIB if visualize else 0,
            "preview_bytes_per_displayed_particle": 280 if visualize else 0,
            "report_and_manifest_bytes": report_bytes,
            "disk_safety_margin_fraction": 0.20,
        },
        "formulas": {
            "peak_memory": (
                "fixed_runtime_bytes + matched_count * compact_array_bytes_per_particle + "
                "min(bounded_query_batch, matched_count) * (k_neighbors + 1) * 16 + "
                "cs_input_bytes * cs_input_resident_fraction + star_input_bytes * star_input_resident_factor + "
                "[visualize] (preview_memory_fixed_bytes + min(matched_count, preview_max_plotted_particles) * "
                "preview_memory_bytes_per_plotted_particle)"
            ),
            "bundle_disk": (
                "matched_count * npz_bytes_per_particle + matched_count * "
                "raw_csv_bytes_per_particle + preview estimate + reports"
            ),
        },
    }
