"""CLI progress, timing, and memory diagnostics helpers."""

from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import ctypes
import os
from pathlib import Path
import sys
from time import perf_counter
from typing import Any, Iterator, TextIO


class ProgressReporter:
    """Small stderr reporter with optional timing and RSS profiling."""

    def __init__(
        self,
        *,
        command: str = "run",
        quiet: bool = False,
        verbose: bool = False,
        profile_time: bool = False,
        profile_memory: bool = False,
        stream: TextIO | None = None,
    ) -> None:
        self.command = command
        self.quiet = quiet
        self.verbose = verbose
        self.profile_time = profile_time
        self.profile_memory = profile_memory
        self.stream = stream or sys.stderr
        self._start_time = perf_counter()
        self._stages: list[dict[str, Any]] = []
        self._memory_samples: list[dict[str, Any]] = []
        self._memory_backend: str | None = None

    def stage_start(
        self,
        label: str,
        *,
        stage_number: int | None = None,
        stage_total: int | None = None,
        detail: str | None = None,
    ) -> None:
        """Emit a stage-start progress line."""

        if self.quiet:
            return
        prefix = (
            f"stage {stage_number}/{stage_total}: "
            if stage_number is not None and stage_total is not None
            else ""
        )
        suffix = f": {detail}" if detail else ""
        self._emit(f"{prefix}{label}{suffix}")

    def stage_done(self, label: str, *, elapsed_sec: float, detail: str | None = None) -> None:
        """Emit a stage completion line."""

        if self.quiet:
            return
        suffix = f": {detail}" if detail else ""
        self._emit(f"{label} completed in {elapsed_sec:.3f} s{suffix}")

    @contextmanager
    def timed_stage(
        self,
        stage: str,
        *,
        label: str,
        stage_number: int | None = None,
        stage_total: int | None = None,
        detail: str | None = None,
        done_label: str | None = None,
        **metadata: Any,
    ) -> Iterator[dict[str, Any]]:
        """Time a stage and record stable metadata."""

        self.stage_start(
            label,
            stage_number=stage_number,
            stage_total=stage_total,
            detail=detail,
        )
        start = perf_counter()
        stage_metadata: dict[str, Any] = dict(metadata)
        try:
            yield stage_metadata
        finally:
            elapsed = perf_counter() - start
            record = {"stage": stage, "elapsed_sec": elapsed}
            record.update(stage_metadata)
            self._stages.append(record)
            self.stage_done(done_label or label, elapsed_sec=elapsed, detail=_detail(stage_metadata))

    def batch_update(self, name: str, processed: int, total: int) -> None:
        """Emit an approximate batch update when measurable."""

        if self.quiet or total <= 0:
            return
        elapsed = max(perf_counter() - self._start_time, 1e-9)
        rate = processed / elapsed
        eta = (total - processed) / rate if rate > 0 and processed < total else None
        percent = 100.0 * processed / total
        eta_text = f", ETA ~{eta:.0f}s" if eta is not None else ""
        self._emit(f"{name}: {processed} / {total} rows, {percent:.1f}%{eta_text}")

    def info(self, message: str, *, verbose_only: bool = False) -> None:
        """Emit an informational progress line."""

        if self.quiet or (verbose_only and not self.verbose):
            return
        self._emit(message)

    def warning(self, message: str) -> None:
        """Emit a warning even in quiet mode."""

        self._emit(f"warning: {message}")

    def sample_memory(self, stage: str) -> None:
        """Record and optionally emit an RSS memory sample."""

        if not self.profile_memory and not self.verbose:
            return
        rss_bytes, backend = _current_rss_bytes()
        if self._memory_backend is None:
            self._memory_backend = backend
        sample = {
            "stage": stage,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "rss_bytes": rss_bytes,
            "rss_mib": rss_bytes / (1024 * 1024),
        }
        if self.profile_memory:
            self._memory_samples.append(sample)
        if self.verbose and not self.quiet:
            self._emit(f"rss after {stage}: {sample['rss_mib']:.1f} MiB")

    def timing_profile_payload(self) -> dict[str, Any]:
        """Return the JSON-safe timing profile payload."""

        return {
            "artifact_type": "run_timing_profile",
            "schema_version": "1",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "total_elapsed_sec": perf_counter() - self._start_time,
            "stages": self._stages,
        }

    def memory_profile_payload(self) -> dict[str, Any]:
        """Return the JSON-safe memory profile payload."""

        peak = max((sample["rss_bytes"] for sample in self._memory_samples), default=0)
        return {
            "artifact_type": "run_memory_profile",
            "schema_version": "1",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "metric": "rss",
            "units": "bytes",
            "backend": self._memory_backend,
            "peak_rss_bytes": peak,
            "peak_rss_mib": peak / (1024 * 1024),
            "samples": self._memory_samples,
        }

    def _emit(self, message: str) -> None:
        self.stream.write(f"[{self.command}] {message}\n")
        self.stream.flush()


def _detail(metadata: dict[str, Any]) -> str | None:
    parts = []
    for key in ("row_count", "matched_count", "path"):
        if key in metadata:
            parts.append(f"{key}={metadata[key]}")
    return ", ".join(parts) if parts else None


def _current_rss_bytes() -> tuple[int, str]:
    try:
        import psutil  # type: ignore

        return int(psutil.Process().memory_info().rss), "psutil"
    except Exception:
        pass
    try:
        if os.name == "nt":
            return _current_rss_bytes_windows(), "windows_psapi"
        if sys.platform.startswith("linux"):
            return _current_rss_bytes_linux_proc(), "linux_proc_statm"
    except Exception:
        pass
    return 0, "unavailable"


def _current_rss_bytes_linux_proc() -> int:
    with Path("/proc/self/statm").open("r", encoding="utf-8") as handle:
        fields = handle.read().split()
    if len(fields) < 2:
        raise RuntimeError("Could not read RSS from /proc/self/statm")
    return int(fields[1]) * int(os.sysconf("SC_PAGE_SIZE"))


def _current_rss_bytes_windows() -> int:
    class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
        _fields_ = [
            ("cb", ctypes.c_ulong),
            ("PageFaultCount", ctypes.c_ulong),
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
    process = ctypes.windll.kernel32.GetCurrentProcess()
    ok = ctypes.windll.psapi.GetProcessMemoryInfo(
        process,
        ctypes.byref(counters),
        counters.cb,
    )
    if not ok:
        raise RuntimeError("Could not read RSS from Windows Process Status API")
    return int(counters.WorkingSetSize)
