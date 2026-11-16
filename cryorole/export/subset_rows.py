"""Exact integer source-row selection shared by native subset writers."""
from __future__ import annotations

from typing import Sequence

import numpy as np


def validated_subset_rows(row_ids: Sequence[int], *, source_row_count: int, row_id_field: str, source_type: str) -> np.ndarray:
    rows = np.asarray(row_ids)
    if rows.ndim != 1 or (rows.size and (
        rows.dtype.kind not in 'iu' or any(isinstance(row, (bool, np.bool_)) for row in row_ids)
    )):
        raise ValueError(f'{row_id_field} must contain one-dimensional integer source row IDs')
    if np.unique(rows).size != rows.size:
        raise ValueError(f'{row_id_field} contains duplicate source row IDs')
    invalid = rows[(rows < 0) | (rows >= source_row_count)]
    if invalid.size:
        raise ValueError(f'{row_id_field} contains out-of-bounds row IDs for {source_type} source '
                         f'with {source_row_count} rows: {invalid[:5].tolist()}')
    return rows.astype(np.int64, copy=False)
