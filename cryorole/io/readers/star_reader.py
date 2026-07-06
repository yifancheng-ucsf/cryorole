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
