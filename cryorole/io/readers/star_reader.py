"""Minimal loop-oriented RELION STAR reader.

This parser supports the simple loop-based RELION STAR layout used by Phase 1
and Phase 2 tests. It is raw-only: it does not interpret Euler angles, resolve
particle identity, or assume row-order correspondence. Parsed particle row
indices are named ``source_row_id`` so downstream normalization/debugging can
trace rows back to the source table without using row order as an identity key.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from cryorole.io.import_report import ImportReport


@dataclass(frozen=True)
class RelionStarData:
    """Raw RELION STAR data blocks."""

    path: str
    optics: pd.DataFrame
    particles: pd.DataFrame
    report: ImportReport


def _finalize_loop(headers: list[str], rows: list[list[str]]) -> pd.DataFrame:
    if not headers:
        return pd.DataFrame()
    return pd.DataFrame(rows, columns=headers)


def read_relion_star(path: str | Path) -> RelionStarData:
    """Parse a RELION STAR file into raw optics and particle tables."""

    star_path = Path(path)
    if not star_path.is_file():
        raise FileNotFoundError(f"RELION STAR file does not exist: {star_path}")

    tables: dict[str, pd.DataFrame] = {}
    current_block: str | None = None
    in_loop = False
    headers: list[str] = []
    rows: list[list[str]] = []

    def flush() -> None:
        nonlocal headers, rows, in_loop, current_block
        if current_block and headers:
            tables[current_block] = _finalize_loop(headers, rows)
        headers = []
        rows = []
        in_loop = False

    with star_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            if stripped.startswith("data_"):
                flush()
                current_block = stripped
                continue

            if stripped == "loop_":
                in_loop = True
                headers = []
                rows = []
                continue

            if in_loop and stripped.startswith("_"):
                headers.append(stripped.split()[0])
                continue

            if in_loop:
                if not headers:
                    raise ValueError(f"STAR row before loop headers at line {line_number}")
                row = stripped.split()
                if len(row) != len(headers):
                    raise ValueError(
                        f"STAR row at line {line_number} has {len(row)} fields; "
                        f"expected {len(headers)}"
                    )
                rows.append(row)
                continue

    flush()

    particles = tables.get("data_particles", pd.DataFrame())
    optics = tables.get("data_optics", pd.DataFrame())
    particles.index.name = "source_row_id"
    optics.index.name = "source_row_id"
    report = ImportReport(
        path=str(star_path),
        source_type="relion",
        row_count=len(particles),
        columns=tuple(particles.columns),
    )
    return RelionStarData(
        path=str(star_path),
        optics=optics,
        particles=particles,
        report=report,
    )


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

    star_path = Path(path)
    if not star_path.is_file():
        raise FileNotFoundError(f"RELION STAR file does not exist: {star_path}")
    requested = tuple(dict.fromkeys(columns))
    current_block: str | None = None
    in_loop = False
    headers: list[str] = []
    selected_positions: tuple[tuple[int, str], ...] | None = None
    retained: dict[str, list[str]] = {}
    particle_headers: tuple[str, ...] = ()
    particle_row_count = 0
    particle_loop_seen = False

    def close_loop() -> None:
        nonlocal particle_headers, particle_loop_seen
        if current_block == "data_particles" and headers:
            if particle_loop_seen:
                raise ValueError("RELION STAR contains multiple data_particles loops")
            particle_headers = tuple(headers)
            particle_loop_seen = True

    with star_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.startswith("data_"):
                close_loop()
                current_block = stripped
                in_loop = False
                headers = []
                selected_positions = None
                continue
            if stripped == "loop_":
                close_loop()
                in_loop = True
                headers = []
                selected_positions = None
                continue
            if in_loop and stripped.startswith("_"):
                headers.append(stripped.split()[0])
                continue
            if not in_loop:
                continue
            if not headers:
                raise ValueError(f"STAR row before loop headers at line {line_number}")
            row = stripped.split()
            if len(row) != len(headers):
                raise ValueError(
                    f"STAR row at line {line_number} has {len(row)} fields; "
                    f"expected {len(headers)}"
                )
            if current_block != "data_particles":
                continue
            if selected_positions is None:
                selected_positions = tuple(
                    (index, header) for index, header in enumerate(headers)
                    if header in requested
                )
                retained = {header: [] for _, header in selected_positions}
            for index, header in selected_positions:
                retained[header].append(row[index])
            particle_row_count += 1
    close_loop()

    particles = pd.DataFrame(
        {column: retained[column] for column in requested if column in retained},
        index=pd.RangeIndex(particle_row_count, name="source_row_id"),
    )
    report = ImportReport(
        path=str(star_path),
        source_type="relion",
        row_count=particle_row_count,
        columns=particle_headers,
    )
    return RelionStarData(
        path=str(star_path),
        optics=pd.DataFrame(),
        particles=particles,
        report=report,
    )
