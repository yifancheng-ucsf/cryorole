"""Validation helpers for normalization inputs."""

from __future__ import annotations

import pandas as pd


def require_columns(table: pd.DataFrame, columns: tuple[str, ...], *, context: str) -> None:
    """Raise if a raw table is missing required columns."""

    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{context} missing required columns: {missing}")

