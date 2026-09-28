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

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np
import pandas as pd

from cryorole.io.import_report import ImportReport

PARTICLES_BLOCK = "data_particles"
OPTICS_BLOCK = "data_optics"

_QUOTED_FIELDS = re.compile(r'"[^"]*"|\'[^\']*\'|\S+')


def split_star_fields(stripped: str) -> list[str]:
    """Split one STAR data row into fields.

    Plain rows split on whitespace. Rows containing quotes keep a quoted value
    (``"a b"`` or ``'a b'``) as one field, quotes included, so that values are
    reproduced verbatim. Every cryoROLE STAR reader uses this function.
    """

    if '"' not in stripped and "'" not in stripped:
        return stripped.split()
    return _QUOTED_FIELDS.findall(stripped)


def unquote_star_value(value: str) -> str:
    """Remove one pair of matching surrounding quotes, if present."""

    text = value.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}:
        return text[1:-1]
    return text


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
        fields = split_star_fields(stripped)
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


@dataclass(frozen=True)
class StarParticleIndex:
    """Byte-level index of the ``data_particles`` loop of one STAR file.

    Row ``i`` of the particle table is ``bytes[row_start[i]:row_end[i]]`` of the
    file (including its line terminator). ``prefix_end`` is the byte offset where
    the first particle row starts (everything before it, including all headers
    and other blocks, is copied verbatim by subset writers) and ``suffix_start``
    where the text after the last particle row starts.
    """

    path: str
    headers: tuple[str, ...]
    columns: dict[str, list[str]]
    row_start: np.ndarray
    row_end: np.ndarray
    prefix_end: int
    suffix_start: int
    optics: dict[str, list[str]] = field(default_factory=dict)

    @property
    def row_count(self) -> int:
        return int(len(self.row_start))

    def column(self, name: str) -> list[str]:
        if name not in self.columns:
            raise KeyError(name)
        return self.columns[name]


def index_star_particles(
    path: str | Path,
    *,
    columns: Iterable[str] | None = None,
) -> StarParticleIndex:
    """Index the particle loop by byte offset, retaining only ``columns``.

    Uses the shared structure parser (same ``data_particles`` rule and
    multiple-loop rejection as every other reader). ``columns=None`` retains
    none; unknown requested columns are simply absent from ``columns``.
    """

    star_path = Path(path)
    if not star_path.is_file():
        raise FileNotFoundError(f"RELION STAR file does not exist: {star_path}")
    requested = tuple(dict.fromkeys(columns or ()))
    offsets: list[int] = []

    def decoded_lines(handle):
        position = 0
        for raw in handle:
            offsets.append(position)
            position += len(raw)
            yield raw.decode("utf-8")
        offsets.append(position)

    collector = _TableCollector(particle_columns=requested, keep_optics=True)
    particle_rows: list[list[int]] = []

    def on_loop_end(block, _headers, _start, _end, indices):
        if block == PARTICLES_BLOCK:
            particle_rows.append(list(indices))

    with star_path.open("rb") as handle:
        _scan_star(
            decoded_lines(handle),
            on_loop=collector.on_loop,
            on_loop_end=on_loop_end,
            record_line_indices=True,
        )
    if collector.particles is None:
        raise ValueError(
            f"STAR particle loop cannot be found in {star_path}: the file has no 'data_particles' block. "
            "cryoROLE supports RELION 3.1+ STAR files with a data_particles table."
        )
    line_offsets = np.asarray(offsets, dtype=np.int64)
    rows = np.asarray(particle_rows[0] if particle_rows else [], dtype=np.int64)
    if rows.size:
        row_start = line_offsets[rows]
        row_end = line_offsets[rows + 1]
        prefix_end = int(row_start[0])
        suffix_start = int(row_end[-1])
    else:
        row_start = np.empty(0, dtype=np.int64)
        row_end = np.empty(0, dtype=np.int64)
        prefix_end = suffix_start = int(line_offsets[-1])
    optics = {}
    if collector.optics is not None:
        optics = {name: list(values) for name, values in collector.optics.values.items()}
    return StarParticleIndex(
        path=str(star_path),
        headers=collector.particle_headers,
        columns={name: collector.particles.values[name] for name in requested if name in collector.particles.values},
        row_start=row_start,
        row_end=row_end,
        prefix_end=prefix_end,
        suffix_start=suffix_start,
        optics=optics,
    )


def _replace_fields_in_place(line: str, replacements: dict[int, str], *, expected_fields: int) -> str:
    """Replace whole fields by position, keeping every other byte of the line (spacing included)."""

    spans = [match.span() for match in _QUOTED_FIELDS.finditer(line)]
    if len(spans) != expected_fields:
        raise ValueError(f"STAR row has {len(spans)} fields; expected {expected_fields}: {line.strip()[:80]}")
    pieces = []
    cursor = 0
    for position, (start, end) in enumerate(spans):
        if position in replacements:
            pieces.append(line[cursor:start])
            pieces.append(replacements[position])
            cursor = end
    pieces.append(line[cursor:])
    return "".join(pieces)


def write_star_particle_subset(
    index: StarParticleIndex,
    row_ids: Sequence[int],
    output_path: str | Path,
    *,
    header_comment: str | None = None,
    replace_columns: dict[str, Sequence[str]] | None = None,
) -> Path:
    """Write ``row_ids`` of the indexed particle loop, everything else verbatim.

    Rows are copied byte for byte unless ``replace_columns`` gives new values
    (one per output row, in ``row_ids`` order) for named particle columns; only
    those fields are rewritten. Comment/blank lines interleaved with particle
    rows are not copied.
    """

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    replace_columns = replace_columns or {}
    positions = {name: index.headers.index(name) for name in replace_columns}
    for name, values in replace_columns.items():
        if len(values) != len(row_ids):
            raise ValueError(f"replacement values for {name} must match the selected row count")
    with open(index.path, "rb") as source, output.open("wb") as target:
        if header_comment:
            for comment_line in header_comment.splitlines():
                target.write(f"# {comment_line}\n".encode("utf-8"))
        target.write(source.read(index.prefix_end))
        for out_row, row_id in enumerate(row_ids):
            start = int(index.row_start[row_id])
            end = int(index.row_end[row_id])
            source.seek(start)
            raw = source.read(end - start)
            if positions:
                raw = _replace_fields_in_place(
                    raw.decode("utf-8"),
                    {column_index: str(replace_columns[name][out_row]) for name, column_index in positions.items()},
                    expected_fields=len(index.headers),
                ).encode("utf-8")
            if not raw.endswith((b"\n", b"\r")):
                raw += b"\n"
            target.write(raw)
        source.seek(index.suffix_start)
        target.write(source.read())
    return output


class _StopScan(Exception):
    pass


def read_star_particle_headers(path: str | Path) -> tuple[str, ...]:
    """Return the ``data_particles`` loop headers without reading the rows."""

    star_path = Path(path)
    if not star_path.is_file():
        raise FileNotFoundError(f"RELION STAR file does not exist: {star_path}")
    found: list[tuple[str, ...]] = []

    def on_loop(block, headers, _start):
        if block == PARTICLES_BLOCK:
            found.append(headers)
            raise _StopScan
        return None

    with star_path.open("r", encoding="utf-8") as handle:
        try:
            _scan_star(handle, on_loop=on_loop)
        except _StopScan:
            pass
    if not found:
        raise ValueError(
            f"STAR particle loop cannot be found in {star_path}: the file has no 'data_particles' block."
        )
    return found[0]
