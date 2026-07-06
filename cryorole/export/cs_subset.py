"""Source-preserving CryoSPARC .cs row subset writer."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class CsSubsetResult:
    """Summary of a CryoSPARC .cs subset export."""

    source_path: str
    output_path: str
    selected_row_count: int
    source_row_count: int
    row_id_field: str
    uid_output_path: str | None = None
    warnings: tuple[str, ...] = ()


def write_cryosparc_cs_subset(
    source_path: str | Path,
    output_path: str | Path,
    row_ids: Sequence[int],
    *,
    row_id_field: str,
    overwrite: bool = False,
    uid_output_path: str | Path | None = None,
) -> CsSubsetResult:
    """Write a native CryoSPARC structured-array subset without changing fields."""

    source = Path(source_path)
    output = Path(output_path)
    if not source.is_file():
        raise ValueError(f"CryoSPARC CS source file does not exist: {source}")
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {output}")
    array = np.load(source, allow_pickle=False)
    if array.dtype.names is None:
        raise ValueError("CryoSPARC CS file must be a structured array")

    _validate_row_ids(row_ids, source_row_count=len(array), row_id_field=row_id_field)
    subset = array[np.asarray(row_ids, dtype=np.int64)]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as handle:
        np.save(handle, subset)

    written_uid_path: str | None = None
    if "uid" in array.dtype.names:
        uid_path = Path(uid_output_path) if uid_output_path is not None else None
        if uid_path is not None:
            if uid_path.exists() and not overwrite:
                raise FileExistsError(f"Output path already exists: {uid_path}")
            uid_path.parent.mkdir(parents=True, exist_ok=True)
            with uid_path.open("w", encoding="utf-8") as handle:
                for uid in subset["uid"]:
                    handle.write(f"{uid}\n")
            written_uid_path = str(uid_path)

    return CsSubsetResult(
        source_path=str(source),
        output_path=str(output),
        selected_row_count=len(row_ids),
        source_row_count=len(array),
        row_id_field=row_id_field,
        uid_output_path=written_uid_path,
    )


def _validate_row_ids(
    row_ids: Sequence[int],
    *,
    source_row_count: int,
    row_id_field: str,
) -> None:
    if len(set(row_ids)) != len(row_ids):
        raise ValueError(f"{row_id_field} contains duplicate source row IDs")
    out_of_bounds = [row_id for row_id in row_ids if row_id < 0 or row_id >= source_row_count]
    if out_of_bounds:
        raise ValueError(
            f"{row_id_field} contains out-of-bounds row IDs for CS source "
            f"with {source_row_count} rows: {out_of_bounds[:5]}"
        )
