"""Source-preserving RELION STAR row subset writer."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


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
) -> StarSubsetResult:
    """Write a STAR file containing selected particle-loop rows only.

    The writer preserves the original text outside the particle loop and copies
    selected particle rows verbatim. It does not interpret or modify RELION
    pose, origin, CTF, optics, or image-name values.
    """

    source = Path(source_path)
    output = Path(output_path)
    if not source.is_file():
        raise ValueError(f"RELION STAR source file does not exist: {source}")
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {output}")

    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    loop = _find_particles_loop(lines)
    _validate_row_ids(row_ids, source_row_count=len(loop.data_lines), row_id_field=row_id_field)
    selected_lines = [loop.data_lines[row_id] for row_id in row_ids]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "".join(lines[: loop.data_start] + selected_lines + lines[loop.data_end :]),
        encoding="utf-8",
    )
    return StarSubsetResult(
        source_path=str(source),
        output_path=str(output),
        selected_row_count=len(row_ids),
        source_row_count=len(loop.data_lines),
        row_id_field=row_id_field,
    )


@dataclass(frozen=True)
class _StarLoop:
    data_start: int
    data_end: int
    headers: tuple[str, ...]
    data_lines: tuple[str, ...]


def _find_particles_loop(lines: list[str]) -> _StarLoop:
    loops: list[_StarLoop] = []
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        header_start = index + 1
        header_index = header_start
        headers: list[str] = []
        while header_index < len(lines):
            stripped = lines[header_index].strip()
            if not stripped or stripped.startswith("#"):
                header_index += 1
                continue
            if stripped.startswith("_"):
                headers.append(stripped.split()[0])
                header_index += 1
                continue
            break
        data_start = header_index
        data_end = data_start
        while data_end < len(lines):
            stripped = lines[data_end].strip()
            if not stripped or stripped.startswith("#"):
                data_end += 1
                continue
            if stripped == "loop_" or stripped.startswith("data_") or stripped.startswith("_"):
                break
            data_end += 1
        data_lines = tuple(
            line for line in lines[data_start:data_end] if line.strip() and not line.strip().startswith("#")
        )
        if headers:
            loops.append(
                _StarLoop(
                    data_start=data_start,
                    data_end=data_end,
                    headers=tuple(headers),
                    data_lines=data_lines,
                )
            )
        index = max(data_end, index + 1)

    for loop in loops:
        if "_rlnImageName" in loop.headers:
            return loop
    for loop in reversed(loops):
        if loop.data_lines and any(header.startswith("_rln") for header in loop.headers):
            return loop
    raise ValueError("STAR particles loop cannot be found")


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
