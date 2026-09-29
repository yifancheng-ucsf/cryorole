"""Process-level parallelism for independent static figures.

Matplotlib's Agg renderer holds the GIL, so figures are rendered in separate
processes. The "spawn" start method is used everywhere: it is the only one that
is safe on macOS, and it gives the same behaviour on Linux.
"""

from __future__ import annotations

import os
from typing import Any

# Below this many particles, starting worker processes (each imports NumPy,
# pandas and Matplotlib, about 1-2 s) costs more than rendering saves.
PARALLEL_FIGURES_MIN_POINTS = 100_000
# A worker's own footprint before data: interpreter, imports, figure buffers.
WORKER_FIXED_BYTES = 200 * 1024 * 1024
# Keep this share of the available memory free for the parent and the system.
MEMORY_RESERVE_FRACTION = 1 / 3


def resolve_figure_jobs(
    requested: int | str | None,
    *,
    n_tasks: int,
    n_points: int,
    worker_bytes_per_point: float,
    parallel_capable: bool = True,
) -> int:
    """Return how many worker processes to use for ``n_tasks`` independent figure tasks.

    ``requested`` is a positive integer or ``None``/``"auto"``. Auto uses up to
    ``n_tasks`` processes, limited by CPU cores and by the memory available for
    workers of ``WORKER_FIXED_BYTES + n_points * worker_bytes_per_point`` each,
    and 1 below ``PARALLEL_FIGURES_MIN_POINTS`` particles.
    """

    if not parallel_capable or n_tasks <= 1:
        return 1
    if requested not in (None, "auto"):
        jobs = int(requested)
        if jobs < 1:
            raise ValueError("jobs must be a positive integer or 'auto'")
        return min(jobs, n_tasks)
    if n_points < PARALLEL_FIGURES_MIN_POINTS:
        return 1
    jobs = min(n_tasks, os.cpu_count() or 1)
    from cryorole.preflight.environment import _available_memory_bytes

    available, _backend = _available_memory_bytes()
    if available is not None:
        per_worker = WORKER_FIXED_BYTES + n_points * worker_bytes_per_point
        usable = available * (1 - MEMORY_RESERVE_FRACTION)
        jobs = min(jobs, max(1, int(usable // per_worker)))
    return max(1, jobs)


def spawn_pool(max_workers: int) -> Any:
    """A ``ProcessPoolExecutor`` using the "spawn" start method."""

    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor

    return ProcessPoolExecutor(max_workers=max_workers, mp_context=multiprocessing.get_context("spawn"))
