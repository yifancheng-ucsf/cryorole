"""Raw import reports."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ImportReport:
    """Summary of a parsed raw input file."""

    path: str
    source_type: str
    row_count: int
    columns: tuple[str, ...]
    status: str = "ok"

