"""Side-effect-free checks of the machine a run would use (preflight only).

Reports package versions, whether the output location can be written, whether
it already exists, and available memory next to the run estimate. Nothing is
created: writability is checked with ``os.access`` on the nearest existing
parent directory.
"""

from __future__ import annotations

import os
import platform
import sys
from pathlib import Path
from typing import Any

_PACKAGES = ("numpy", "scipy", "pandas", "matplotlib")


def _package_versions() -> dict[str, str | None]:
    from importlib.metadata import PackageNotFoundError, version

    versions: dict[str, str | None] = {}
    for name in ("cryorole", *_PACKAGES):
        try:
            versions[name] = version(name)
        except PackageNotFoundError:
            versions[name] = None
    if versions["cryorole"] is None:
        from cryorole import __version__

        versions["cryorole"] = __version__
    return versions


def _available_memory_bytes() -> tuple[int | None, str]:
    try:
        import psutil  # type: ignore

        return int(psutil.virtual_memory().available), "psutil"
    except Exception:  # noqa: BLE001 - optional dependency
        pass
    meminfo = Path("/proc/meminfo")
    if meminfo.is_file():
        for line in meminfo.read_text(encoding="utf-8").splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024, "linux_proc_meminfo"
    return None, "unavailable"


def _nearest_existing(path: Path) -> Path:
    probe = path
    while not probe.exists() and probe.parent != probe:
        probe = probe.parent
    return probe


def check_environment(output_dir: str | Path, *, estimated_peak_memory_bytes: int | None = None) -> dict[str, Any]:
    """Return an ``environment`` record plus ``errors``/``warnings`` for preflight."""

    target = Path(output_dir).expanduser().resolve()
    probe = _nearest_existing(target if not target.exists() else target.parent)
    writable = probe.is_dir() and os.access(probe, os.W_OK | os.X_OK)
    output_exists = target.exists()
    memory, memory_backend = _available_memory_bytes()
    errors: list[str] = []
    warnings: list[str] = []
    if not writable:
        errors.append(
            f"[OUTPUT_NOT_WRITABLE] The output location cannot be written ({probe}); "
            "choose another --output-dir."
        )
    if output_exists:
        warnings.append(
            f"[OUTPUT_EXISTS] {target} already exists; `run` will refuse it unless --overwrite is given "
            "(or choose another --output-dir)."
        )
    if memory is not None and estimated_peak_memory_bytes and memory < estimated_peak_memory_bytes:
        warnings.append(
            f"[LOW_MEMORY] Available memory ({memory / 2**20:.0f} MiB) is below the estimated peak "
            f"({estimated_peak_memory_bytes / 2**20:.0f} MiB); consider a smaller --density-query-batch-size "
            "or another machine."
        )
    return {
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "packages": _package_versions(),
            "output_dir": str(target),
            "output_exists": output_exists,
            "output_parent_checked": str(probe),
            "output_writable": writable,
            "available_memory_bytes": memory,
            "available_memory_backend": memory_backend,
        },
        "errors": errors,
        "warnings": warnings,
    }
