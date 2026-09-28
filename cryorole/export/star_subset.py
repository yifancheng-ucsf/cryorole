"""Source-preserving RELION STAR row subset writer."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from cryorole.io.readers.star_reader import locate_particle_loop


@dataclass(frozen=True)
class StarSubsetResult:
    """Summary of a RELION STAR subset export."""

    source_path: str
    output_path: str
    selected_row_count: int
    source_row_count: int
    row_id_field: str
    warnings: tuple[str, ...] = ()


def write_relion_star_subset(
    source_path: str | Path,
    output_path: str | Path,
    row_ids: Sequence[int],
    *,
    row_id_field: str,
    overwrite: bool = False,
    expected_source_row_count: int | None = None,
    required_headers: Sequence[str] = (),
) -> StarSubsetResult:
    """Write a STAR file containing selected particle-loop rows only.

    The writer preserves the original text outside the particle loop and copies
    selected particle rows verbatim. It does not interpret or modify RELION
    pose, origin, CTF, optics, or image-name values.

    The particle loop is resolved with the same ``data_particles`` block rule
    as the run-time reader, so row IDs recorded by the run index the same rows
    here. When ``expected_source_row_count`` is given (the particle count the
    run recorded for this source) a mismatch fails loudly instead of exporting
    rows from a different table. ``required_headers`` (for example the pose
    columns) must all be present in the resolved loop.
    """

    source = Path(source_path)
    output = Path(output_path)
    if not source.is_file():
        raise ValueError(f"RELION STAR source file does not exist: {source}")
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {output}")

    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    loop = locate_particle_loop(lines)
    if expected_source_row_count is not None and loop.row_count != int(expected_source_row_count):
        raise ValueError(
            f"STAR particle loop in {source} has {loop.row_count} rows but the run recorded "
            f"{int(expected_source_row_count)}; refusing to export. The source file may have been "
            "replaced, or it is not the file the run used."
        )
    missing_headers = [header for header in required_headers if header not in loop.headers]
    if missing_headers:
        raise ValueError(
            f"STAR particle loop in {source} is missing required columns {missing_headers}; "
            "refusing to export."
        )
    _validate_row_ids(row_ids, source_row_count=loop.row_count, row_id_field=row_id_field)
    selected_lines = [_with_newline(lines[loop.data_line_indices[row_id]]) for row_id in row_ids]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "".join(lines[: loop.data_start] + selected_lines + lines[loop.data_end :]),
        encoding="utf-8",
    )
    return StarSubsetResult(
        source_path=str(source),
        output_path=str(output),
        selected_row_count=len(row_ids),
        source_row_count=loop.row_count,
        row_id_field=row_id_field,
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
            f"{row_id_field} contains out-of-bounds row IDs for STAR source "
            f"with {source_row_count} rows: {out_of_bounds[:5]}"
        )


def _with_newline(line: str) -> str:
    return line if line.endswith(("\n", "\r")) else line + "\n"
