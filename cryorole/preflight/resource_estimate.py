"""Transparent deterministic resource estimates for production run preflight."""

from __future__ import annotations

from pathlib import Path
import shutil


MIB = 1024 * 1024


def estimate_run_resources(
    *,
    matched_count: int,
    k_neighbors: int,
    query_batch_size: int,
    raw_csv: bool,
    visualize: bool,
    output_dir: str | Path,
) -> dict[str, object]:
    """Return a conservative estimate with every formula input exposed."""

    n = max(int(matched_count), 0)
    k = max(int(k_neighbors), 1)
    batch = max(int(query_batch_size), 1)
    bounded_batch = min(batch, max(1, 2_000_000 // (k + 1)))
    fixed_bytes = 96 * MIB
    compact_arrays_bytes = n * 360
    knn_query_bytes = bounded_batch * (k + 1) * 16
    estimated_peak = fixed_bytes + compact_arrays_bytes + knn_query_bytes
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
            "compact_array_bytes_per_particle": 360,
            "knn_distance_and_index_bytes_per_element": 16,
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
                "bounded_query_batch * (k_neighbors + 1) * 16"
            ),
            "bundle_disk": (
                "matched_count * npz_bytes_per_particle + matched_count * "
                "raw_csv_bytes_per_particle + preview estimate + reports"
            ),
        },
    }
