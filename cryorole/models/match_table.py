"""Particle match table and report models."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass
class MatchTable:
    """Pairing result between two domain pose tables."""

    data: pd.DataFrame

    def __post_init__(self) -> None:
        required = {"particle_key", "domain_a_row", "domain_b_row", "match_status"}
        missing = required - set(self.data.columns)
        if missing:
            raise ValueError(f"MatchTable missing required columns: {sorted(missing)}")


@dataclass(frozen=True)
class MatchReport:
    """Summary of matching counts and overlap status."""

    matched_count: int
    unmatched_a: int
    unmatched_b: int
    overlap_ratio: float
    status: str
    match_key: str | None = None
    ref_row_count: int = 0
    mov_row_count: int = 0
    matched_row_count: int = 0
    dropped_ref_only_count: int = 0
    dropped_mov_only_count: int = 0
    matched_rows_reordered: bool = False
    warnings: tuple[str, ...] = ()
