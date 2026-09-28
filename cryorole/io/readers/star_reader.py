"""Loop-oriented RELION STAR reader.

This module owns the single STAR structure parser used by every cryoROLE path
that needs to know which rows belong to the particle table: the production
run reader, the metadata-selection reader, and the source-preserving STAR
subset writer. Sharing one parser guarantees that "source row N" means the
same particle everywhere.

Particle-table resolution rule (all callers):

* the particle loop is the loop inside the block named ``data_particles``;
* a file with more than one ``data_particles`` loop (a repeated block or two
  loops inside the block) is rejected, never silently merged or truncated;
* a file without a ``data_particles`` block has no particle table.

The reader is raw-only: it does not interpret Euler angles, resolve particle
identity, or assume row-order correspondence. Parsed particle row indices are
named ``source_row_id`` so downstream normalization/debugging can trace rows
back to the source table without using row order as an identity key.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

import pandas as pd

from cryorole.io.import_report import ImportReport

PARTICLES_BLOCK = "data_particles"
OPTICS_BLOCK = "data_optics"


@dataclass(frozen=True)
class RelionStarData:
    """Raw RELION STAR data blocks."""

    path: str
    optics: pd.DataFrame
    particles: pd.DataFrame
    report: ImportReport


@dataclass(frozen=True)
class StarLoopLocation:
    """Line-level location of one STAR loop.

    ``data_line_indices`` are 0-based indices into the scanned line sequence of
    every non-empty, non-comment data row in the loop, in file order. Row ``i``
    of the loop is ``lines[data_line_indices[i]]``.
    """

    block: str | None
    headers: tuple[str, ...]
    data_start: int
    data_end: int
    data_line_indices: tuple[int, ...]

    @property
    def row_count(self) -> int:
        return len(self.data_line_indices)


class _LoopSink:
    """Callback interface used by :func:`_scan_star` for one loop."""

    def row(self, fields: list[str]) -> None:  # pragma: no cover - interface
        raise NotImplementedError


def _scan_star(
    lines: Iterable[str],
    *,
    on_loop: Callable[[str | None, tuple[str, ...], int], "_LoopSink | None"],
    on_loop_end: Callable[[str | None, tuple[str, ...], int, int, list[int]], None] | None = None,
    record_line_indices: bool = False,
) -> None:
    """Scan STAR text once, dispatching loop rows to per-loop sinks.

    ``on_loop(block, headers, data_start)`` is called once the header section
    of a loop is complete; it returns a sink for the loop's rows, or ``None``
    to validate the rows without retaining them. ``on_loop_end`` receives the
    loop's line span and (optionally) its data line indices.

    A loop ends at the next ``loop_``, ``data_`` block, or a ``_name value``
    line that follows data rows. Every data row is width-checked against the
    loop header regardless of whether the loop is retained.
    """

    block: str | None = None
    state = "none"  # none | headers | rows
    headers: list[str] = []
    sink: _LoopSink | None = None
    data_start = 0
    last_data_end = 0
    line_indices: list[int] = []

    def begin_rows(index: int) -> None:
        nonlocal state, sink, data_start
        state = "rows"
        data_start = index
        sink = on_loop(block, tuple(headers), data_start)

    def end_loop(index: int) -> None:
        nonlocal state, headers, sink, line_indices
        if state == "headers" and headers:
            # Header-only loop (zero rows).
            begin_rows(index)
        if state == "rows" and on_loop_end is not None:
            on_loop_end(block, tuple(headers), data_start, max(last_data_end, data_start), line_indices)
        state = "none"
        headers = []
        sink = None
        line_indices = []

    index = -1
    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("data_"):
            end_loop(index)
            block = stripped.split()[0]
            continue
        if stripped == "loop_":
            end_loop(index)
            state = "headers"
            headers = []
            continue
        if stripped.startswith("_"):
            if state == "headers":
                headers.append(stripped.split()[0])
                continue
            # A key/value pair after data rows (or outside any loop) ends the loop.
            end_loop(index)
            continue
        if state == "none":
            continue
        if state == "headers":
            if not headers:
                raise ValueError(f"STAR row before loop headers at line {index + 1}")
            begin_rows(index)
        fields = stripped.split()
        if len(fields) != len(headers):
            raise ValueError(
                f"STAR row at line {index + 1} has {len(fields)} fields; "
                f"expected {len(headers)}"
            )
        if sink is not None:
            sink.row(fields)
        if record_line_indices:
            line_indices.append(index)
        last_data_end = index + 1
    end_loop(index + 1)


class _ColumnSink(_LoopSink):
    def __init__(self, headers: tuple[str, ...], requested: tuple[str, ...] | None) -> None:
        wanted = headers if requested is None else tuple(h for h in headers if h in set(requested))
        self.positions = tuple((headers.index(h), h) for h in wanted)
        self.values: dict[str, list[str]] = {h: [] for _, h in self.positions}
        self.count = 0

    def row(self, fields: list[str]) -> None:
        for position, header in self.positions:
            self.values[header].append(fields[position])
        self.count += 1


class _TableCollector:
    """Collect particle (and optionally optics) loops with duplicate guards."""

    def __init__(self, *, particle_columns: tuple[str, ...] | None, keep_optics: bool) -> None:
        self.particle_columns = particle_columns
        self.keep_optics = keep_optics
        self.particles: _ColumnSink | None = None
        self.particle_headers: tuple[str, ...] = ()
        self.optics: _ColumnSink | None = None
        self.optics_headers: tuple[str, ...] = ()

    def on_loop(self, block: str | None, headers: tuple[str, ...], _start: int) -> _LoopSink | None:
        if block == PARTICLES_BLOCK:
            if self.particles is not None:
                raise ValueError(
                    "RELION STAR contains multiple data_particles loops; cryoROLE cannot tell "
                    "which one is the particle table. Split or merge the file with RELION first."
                )
            self.particle_headers = headers
            self.particles = _ColumnSink(headers, self.particle_columns)
            return self.particles
        if block == OPTICS_BLOCK and self.keep_optics and self.optics is None:
            self.optics_headers = headers
            self.optics = _ColumnSink(headers, None)
            return self.optics
        return None


def _frame(sink: _ColumnSink | None, order: Sequence[str] | None = None) -> pd.DataFrame:
    if sink is None:
        return pd.DataFrame(index=pd.RangeIndex(0, name="source_row_id"))
    names = [name for name in (order or [h for _, h in sink.positions]) if name in sink.values]
    frame = pd.DataFrame(
        {name: sink.values[name] for name in names},
        index=pd.RangeIndex(sink.count, name="source_row_id"),
    )
    return frame


def _read_tables(
    path: str | Path,
    *,
    particle_columns: tuple[str, ...] | None,
    keep_optics: bool,
) -> tuple[Path, _TableCollector]:
    star_path = Path(path)
    if not star_path.is_file():
        raise FileNotFoundError(f"RELION STAR file does not exist: {star_path}")
    collector = _TableCollector(particle_columns=particle_columns, keep_optics=keep_optics)
    with star_path.open("r", encoding="utf-8") as handle:
        _scan_star(handle, on_loop=collector.on_loop)
    return star_path, collector


def read_relion_star(path: str | Path) -> RelionStarData:
    """Parse a RELION STAR file into raw optics and particle tables.

    Uses the same structure parser and particle-table rule as
    :func:`read_relion_star_particle_columns`, so both give identical
    ``source_row_id`` values for the same file.
    """

    star_path, collector = _read_tables(path, particle_columns=None, keep_optics=True)
    particles = _frame(collector.particles)
    optics = _frame(collector.optics)
    report = ImportReport(
        path=str(star_path),
        source_type="relion",
        row_count=len(particles),
        columns=tuple(particles.columns),
    )
    return RelionStarData(path=str(star_path), optics=optics, particles=particles, report=report)


def read_relion_star_particle_columns(
    path: str | Path,
    *,
    columns: tuple[str, ...],
) -> RelionStarData:
    """Stream the particle loop while retaining only explicitly requested columns.

    The complete particle header and row width are still validated and reported,
    but unrequested cell values are never accumulated in Python containers. This
    is the production preflight/run reader; :func:`read_relion_star` remains the
    compatibility reader for callers that intentionally need all source fields.
    """

    requested = tuple(dict.fromkeys(columns))
    star_path, collector = _read_tables(path, particle_columns=requested, keep_optics=False)
    particles = _frame(collector.particles, order=requested)
    report = ImportReport(
        path=str(star_path),
        source_type="relion",
        row_count=collector.particles.count if collector.particles is not None else 0,
        columns=collector.particle_headers,
    )
    return RelionStarData(path=str(star_path), optics=pd.DataFrame(), particles=particles, report=report)


def locate_star_loops(lines: Sequence[str]) -> list[StarLoopLocation]:
    """Return the location of every loop in already-read STAR ``lines``."""

    loops: list[StarLoopLocation] = []

    def on_loop(_block: str | None, _headers: tuple[str, ...], _start: int) -> None:
        return None

    def on_loop_end(
        block: str | None,
        headers: tuple[str, ...],
        data_start: int,
        data_end: int,
        indices: list[int],
    ) -> None:
        loops.append(
            StarLoopLocation(
                block=block,
                headers=headers,
                data_start=data_start,
                data_end=data_end,
                data_line_indices=tuple(indices),
            )
        )

    _scan_star(lines, on_loop=on_loop, on_loop_end=on_loop_end, record_line_indices=True)
    return loops


def locate_particle_loop(lines: Sequence[str]) -> StarLoopLocation:
    """Resolve the particle loop in STAR ``lines`` using the shared block rule."""

    particle_loops = [loop for loop in locate_star_loops(lines) if loop.block == PARTICLES_BLOCK]
    if len(particle_loops) > 1:
        raise ValueError(
            "RELION STAR contains multiple data_particles loops; cryoROLE cannot tell "
            "which one is the particle table. Split or merge the file with RELION first."
        )
    if not particle_loops:
        raise ValueError(
            "STAR particle loop cannot be found: the file has no 'data_particles' block. "
            "cryoROLE supports RELION 3.1+ STAR files with a data_particles table."
        )
    return particle_loops[0]
