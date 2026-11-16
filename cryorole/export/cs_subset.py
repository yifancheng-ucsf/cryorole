"""Source-preserving CryoSPARC .cs row subset writer."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from cryorole.export.subset_rows import validated_subset_rows
from cryorole.io.writers.publication import ArtifactFiles


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
    expected_source_row_count: int | None = None,
) -> CsSubsetResult:
    """Write a native CryoSPARC structured-array subset without changing fields.

    When ``expected_source_row_count`` (the row count the run recorded for this
    source) is given, a mismatch fails loudly instead of exporting.
    """

    source = Path(source_path)
    output = Path(output_path)
    if not source.is_file():
        raise ValueError(f"CryoSPARC CS source file does not exist: {source}")
    resolved_source = source.expanduser().resolve(strict=True)
    array = np.load(resolved_source, allow_pickle=False)
    if array.dtype.names is None or array.ndim != 1:
        raise ValueError("CryoSPARC CS file must be a one-dimensional structured array")

    if expected_source_row_count is not None and len(array) != int(expected_source_row_count):
        raise ValueError(
            f"CryoSPARC CS source {source} has {len(array)} rows but the run recorded "
            f"{int(expected_source_row_count)}; refusing to export. The source file may have been "
            "replaced, or it is not the file the run used."
        )
    rows = validated_subset_rows(row_ids, source_row_count=len(array), row_id_field=row_id_field, source_type='CS')
    subset = array[rows]
    uid_path = Path(uid_output_path) if uid_output_path is not None and 'uid' in array.dtype.names else None
    destinations = {'metadata': output}
    if uid_path is not None:
        destinations['uids'] = uid_path
    with ArtifactFiles(destinations, overwrite=overwrite, inputs=(source, resolved_source)) as publication:
        with publication.paths['metadata'].open('wb') as handle:
            np.save(handle, subset)
        if uid_path is not None:
            with publication.paths['uids'].open("w", encoding="utf-8") as handle:
                for uid in subset["uid"]:
                    handle.write(f"{uid}\n")

    return CsSubsetResult(
        source_path=str(source),
        output_path=str(output),
        selected_row_count=len(rows),
        source_row_count=len(array),
        row_id_field=row_id_field,
        uid_output_path=str(uid_path) if uid_path is not None else None,
    )
