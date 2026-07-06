"""Compact RELION STAR pre-alignment for cryoROLE run."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shlex
from typing import Any, Mapping, Sequence


AUTO_KEY_CANDIDATES = (
    ("_rlnTomoParticleName",),
    ("_rlnImageName",),
)
LOW_OVERLAP_FRACTION = 0.50
NON_DEFAULT_KEY_MARKERS = (
    "defocus",
    "coordinate",
    "coord",
    "angle",
    "origin",
    "shift",
    "ctf",
)


@dataclass(frozen=True)
class StarRow:
    tokens: tuple[str, ...]
    line: str
    row_number: int
    line_number: int


@dataclass(frozen=True)
class StarLoop:
    block_name: str
    columns: tuple[str, ...]
    row_start_line: int
    row_end_line: int
    rows: tuple[StarRow, ...]

    @property
    def canonical_column_index(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for index, column in enumerate(self.columns):
            out.setdefault(canonical_column_name(column), index)
        return out


@dataclass(frozen=True)
class StarDocument:
    path: Path
    lines: tuple[str, ...]
    loops: tuple[StarLoop, ...]


@dataclass(frozen=True)
class _MatchRows:
    ref_matched: tuple[StarRow, ...]
    mov_matched: tuple[StarRow, ...]
    ref_only: tuple[StarRow, ...]
    mov_only: tuple[StarRow, ...]
    duplicate_ref: tuple[StarRow, ...]
    duplicate_mov: tuple[StarRow, ...]
    match_rows: tuple[dict[str, object], ...]
    duplicate_ref_key_count: int
    duplicate_mov_key_count: int
    warnings: tuple[str, ...]


def align_star_files(
    *,
    ref: str | Path,
    mov: str | Path,
    align_id: str = "default",
    key_columns: Sequence[str] | None = None,
    float_tolerances: Mapping[str, float] | None = None,
    path_mode: str = "exact",
    duplicate_policy: str = "exclude",
    overwrite: bool = False,
) -> dict[str, Any]:
    """Align two RELION-style STAR particle loops and write auditable artifacts."""

    ref_path = Path(ref)
    mov_path = Path(mov)
    _validate_star_input(ref_path, "ref")
    _validate_star_input(mov_path, "mov")
    _validate_align_id(align_id)
    _validate_path_mode(path_mode)
    if duplicate_policy not in {"exclude", "first"}:
        raise ValueError("--duplicate-policy must be exclude or first")
    tolerances = _normalize_float_tolerances(float_tolerances or {})
    output_dir = Path("alignments") / align_id
    _prepare_output_dir(output_dir, overwrite=overwrite)

    ref_doc = parse_star(ref_path)
    mov_doc = parse_star(mov_path)
    key_source = "explicit" if key_columns else "auto"
    resolved_key_columns = (
        list(key_columns)
        if key_columns
        else _resolve_auto_key_columns(ref_doc, mov_doc)
    )
    ref_loop = select_loop(ref_doc, resolved_key_columns)
    mov_loop = select_loop(mov_doc, resolved_key_columns)
    _validate_key_columns(ref_loop, resolved_key_columns, ref_path)
    _validate_key_columns(mov_loop, resolved_key_columns, mov_path)

    key_policy = {
        "key_source": key_source,
        "key_columns": list(resolved_key_columns),
        "auto_candidates": [list(candidate) for candidate in AUTO_KEY_CANDIDATES],
    }
    warnings = list(_explicit_key_warnings(resolved_key_columns, key_source))
    match = _match_rows(
        ref_loop=ref_loop,
        mov_loop=mov_loop,
        key_columns=resolved_key_columns,
        float_tolerances=tolerances,
        path_mode=path_mode,
        duplicate_policy=duplicate_policy,
    )
    warnings.extend(match.warnings)

    matched_count = len(match.ref_matched)
    ref_row_count = len(ref_loop.rows)
    mov_row_count = len(mov_loop.rows)
    if matched_count == 0:
        _cleanup_empty_output_dir(output_dir)
        raise ValueError(
            "Resolved STAR key has zero overlap; specify a different --key, "
            "or adjust --path-mode / --float-tol."
        )
    ref_overlap = matched_count / ref_row_count if ref_row_count else 0.0
    mov_overlap = matched_count / mov_row_count if mov_row_count else 0.0
    low_overlap_warning = 0.0 < min(ref_overlap, mov_overlap) < LOW_OVERLAP_FRACTION
    if low_overlap_warning:
        warnings.append(
            "Low overlap warning: matched rows are below 50% of at least one input."
        )

    output_paths = _output_paths(output_dir)
    _write_lines(output_paths["aligned_ref"], _replace_loop_rows(ref_doc, ref_loop, match.ref_matched))
    _write_lines(output_paths["aligned_mov"], _replace_loop_rows(mov_doc, mov_loop, match.mov_matched))
    _write_lines(output_paths["ref_only"], _replace_loop_rows(ref_doc, ref_loop, match.ref_only))
    _write_lines(output_paths["mov_only"], _replace_loop_rows(mov_doc, mov_loop, match.mov_only))
    _write_lines(output_paths["duplicate_ref"], _replace_loop_rows(ref_doc, ref_loop, match.duplicate_ref))
    _write_lines(output_paths["duplicate_mov"], _replace_loop_rows(mov_doc, mov_loop, match.duplicate_mov))
    _write_match_table(output_paths["match_table"], match.match_rows)

    report: dict[str, Any] = {
        "artifact_type": "align_report",
        "schema_version": "1",
        "status": "ok",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "ref_path": str(ref_path),
        "mov_path": str(mov_path),
        "align_id": align_id,
        "output_dir": str(output_dir),
        "output_paths": {key: str(path) for key, path in output_paths.items()},
        "key_policy": key_policy,
        "key_columns": list(resolved_key_columns),
        "key_source": key_source,
        "path_mode": path_mode,
        "float_tolerances": dict(tolerances),
        "duplicate_policy": duplicate_policy,
        "ref_row_count": ref_row_count,
        "mov_row_count": mov_row_count,
        "matched_count": matched_count,
        "ref_only_count": len(match.ref_only),
        "mov_only_count": len(match.mov_only),
        "duplicate_ref_row_count": len(match.duplicate_ref),
        "duplicate_mov_row_count": len(match.duplicate_mov),
        "duplicate_ref_key_count": match.duplicate_ref_key_count,
        "duplicate_mov_key_count": match.duplicate_mov_key_count,
        "ref_overlap_fraction": ref_overlap,
        "mov_overlap_fraction": mov_overlap,
        "low_overlap_warning": low_overlap_warning,
        "warnings": warnings,
    }
    _write_json(output_paths["align_report"], report)
    return report


def parse_star(path: Path) -> StarDocument:
    lines = tuple(path.read_text(encoding="utf-8", errors="replace").splitlines())
    loops: list[StarLoop] = []
    block = ""
    index = 0
    while index < len(lines):
        stripped = lines[index].strip()
        if stripped.startswith("data_"):
            block = stripped.split()[0]
            index += 1
            continue
        if stripped != "loop_":
            index += 1
            continue
        header_index = index + 1
        while header_index < len(lines) and (
            not lines[header_index].strip()
            or lines[header_index].lstrip().startswith("#")
        ):
            header_index += 1
        columns: list[str] = []
        while header_index < len(lines) and lines[header_index].lstrip().startswith("_"):
            columns.append(lines[header_index].strip().split()[0])
            header_index += 1
        if not columns:
            raise ValueError(f"No STAR loop columns found in {path}")
        row_start = header_index
        rows: list[StarRow] = []
        row_number = 0
        while header_index < len(lines):
            row_text = lines[header_index]
            row_stripped = row_text.strip()
            if row_stripped.startswith("data_") or row_stripped == "loop_":
                break
            if not row_stripped or row_stripped.startswith("#"):
                header_index += 1
                continue
            if row_stripped.startswith("_"):
                break
            tokens = _split_star_row(row_text)
            if len(tokens) != len(columns):
                raise ValueError(
                    f"STAR row token count mismatch at {path}:{header_index + 1}; "
                    f"expected {len(columns)}, got {len(tokens)}"
                )
            rows.append(
                StarRow(
                    tokens=tuple(tokens),
                    line=row_text,
                    row_number=row_number,
                    line_number=header_index,
                )
            )
            row_number += 1
            header_index += 1
        loops.append(
            StarLoop(
                block_name=block,
                columns=tuple(columns),
                row_start_line=row_start,
                row_end_line=header_index,
                rows=tuple(rows),
            )
        )
        index = max(header_index, index + 1)
    return StarDocument(path=path, lines=lines, loops=tuple(loops))


def select_loop(doc: StarDocument, key_columns: Sequence[str]) -> StarLoop:
    preferred = [
        loop
        for loop in doc.loops
        if loop.block_name == "data_particles" and _loop_has_columns(loop, key_columns)
    ]
    if len(preferred) == 1:
        return preferred[0]
    if len(preferred) > 1:
        raise ValueError("Multiple data_particles loops contain the key columns; complex multi-loop selection is not supported")
    candidates = [loop for loop in doc.loops if _loop_has_columns(loop, key_columns)]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        available = sorted({column for loop in doc.loops for column in loop.columns})
        raise ValueError(
            f"Could not find a STAR particle loop containing requested key columns. "
            f"Available columns include: {available[:20]}"
        )
    raise ValueError("Multiple STAR loops contain the key columns; complex multi-loop selection is not supported")


def canonical_column_name(name: str) -> str:
    return name.strip().lstrip("_").lower()


def _validate_star_input(path: Path, label: str) -> None:
    if path.suffix.lower() != ".star":
        raise ValueError(f"cryorole align only supports STAR input for now; {label} is not .star: {path}")
    if not path.is_file():
        raise ValueError(f"{label} STAR file does not exist: {path}")


def _validate_align_id(align_id: str) -> None:
    if not align_id or Path(align_id).name != align_id or align_id in {".", ".."}:
        raise ValueError("--align-id must be a simple directory name")


def _validate_path_mode(path_mode: str) -> None:
    if path_mode in {"exact", "basename"}:
        return
    if re.fullmatch(r"suffix:[1-9][0-9]*", path_mode):
        return
    raise ValueError("--path-mode must be exact, basename, or suffix:N with N > 0")


def _prepare_output_dir(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not path.is_dir():
        raise ValueError(f"Align output path exists and is not a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(f"Align output directory already exists: {path}")


def _cleanup_empty_output_dir(path: Path) -> None:
    try:
        path.rmdir()
        parent = path.parent
        if parent.name == "alignments":
            parent.rmdir()
    except OSError:
        pass


def _split_star_row(line: str) -> list[str]:
    lexer = shlex.shlex(line, posix=False)
    lexer.whitespace_split = True
    lexer.commenters = ""
    return list(lexer)


def _loop_has_columns(loop: StarLoop, key_columns: Sequence[str]) -> bool:
    index = loop.canonical_column_index
    return all(canonical_column_name(column) in index for column in key_columns)


def _validate_key_columns(loop: StarLoop, key_columns: Sequence[str], path: Path) -> None:
    missing = [
        column
        for column in key_columns
        if canonical_column_name(column) not in loop.canonical_column_index
    ]
    if missing:
        raise ValueError(
            f"Missing key columns in {path}: {missing}. "
            f"Available columns include: {list(loop.columns)[:20]}"
        )


def _resolve_auto_key_columns(ref_doc: StarDocument, mov_doc: StarDocument) -> list[str]:
    for candidate in AUTO_KEY_CANDIDATES:
        try:
            select_loop(ref_doc, candidate)
            select_loop(mov_doc, candidate)
        except ValueError:
            continue
        return list(candidate)
    raise ValueError(
        "No safe automatic STAR key was found. Specify an explicit key, for example: "
        "--key _rlnMicrographName _rlnCoordinateX _rlnCoordinateY"
    )


def _normalize_float_tolerances(float_tolerances: Mapping[str, float]) -> dict[str, float]:
    out: dict[str, float] = {}
    for column, tolerance in float_tolerances.items():
        try:
            value = float(tolerance)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid --float-tol value for {column}: {tolerance}") from exc
        if value <= 0:
            raise ValueError(f"--float-tol for {column} must be > 0")
        out[canonical_column_name(column)] = value
    return out


def _explicit_key_warnings(key_columns: Sequence[str], key_source: str) -> tuple[str, ...]:
    if key_source != "explicit":
        return ()
    warnings = []
    for column in key_columns:
        normalized = canonical_column_name(column)
        if any(marker in normalized for marker in NON_DEFAULT_KEY_MARKERS):
            warnings.append(
                f"Explicit key column {column} is a user-selected tracking key, "
                "not a cryoROLE default identity key."
            )
    return tuple(warnings)


def _match_rows(
    *,
    ref_loop: StarLoop,
    mov_loop: StarLoop,
    key_columns: Sequence[str],
    float_tolerances: Mapping[str, float],
    path_mode: str,
    duplicate_policy: str,
) -> _MatchRows:
    ref_keys = [
        _build_key(row, ref_loop, key_columns, float_tolerances, path_mode)
        for row in ref_loop.rows
    ]
    mov_keys = [
        _build_key(row, mov_loop, key_columns, float_tolerances, path_mode)
        for row in mov_loop.rows
    ]
    ref_counts = Counter(ref_keys)
    mov_counts = Counter(mov_keys)
    ref_duplicate_keys = {key for key, count in ref_counts.items() if count > 1}
    mov_duplicate_keys = {key for key, count in mov_counts.items() if count > 1}
    duplicate_keys = ref_duplicate_keys | mov_duplicate_keys
    warnings: list[str] = []
    if duplicate_keys:
        warnings.append(
            f"Duplicate keys detected; duplicate-policy {duplicate_policy} was applied."
        )

    ref_duplicate_rows: list[StarRow] = []
    mov_duplicate_rows: list[StarRow] = []

    if duplicate_policy == "exclude":
        ref_candidate_indices = [
            index for index, key in enumerate(ref_keys) if key not in duplicate_keys
        ]
        mov_candidate_indices = [
            index for index, key in enumerate(mov_keys) if key not in duplicate_keys
        ]
        ref_duplicate_rows = [
            row for row, key in zip(ref_loop.rows, ref_keys) if key in duplicate_keys
        ]
        mov_duplicate_rows = [
            row for row, key in zip(mov_loop.rows, mov_keys) if key in duplicate_keys
        ]
    elif duplicate_policy == "first":
        warnings.append("duplicate-policy first kept the first occurrence and excluded later duplicates.")
        ref_candidate_indices, ref_duplicate_rows = _first_candidate_indices(ref_loop.rows, ref_keys)
        mov_candidate_indices, mov_duplicate_rows = _first_candidate_indices(mov_loop.rows, mov_keys)
    else:
        raise ValueError("--duplicate-policy must be exclude or first")

    mov_by_key: dict[tuple[str, ...], int] = {}
    for index in mov_candidate_indices:
        mov_by_key[mov_keys[index]] = index

    used_mov: set[int] = set()
    ref_matched: list[StarRow] = []
    mov_matched: list[StarRow] = []
    ref_only: list[StarRow] = []
    match_table: list[dict[str, object]] = []
    for ref_index in ref_candidate_indices:
        key = ref_keys[ref_index]
        mov_index = mov_by_key.get(key)
        if mov_index is None or mov_index in used_mov:
            ref_only.append(ref_loop.rows[ref_index])
            continue
        output_row = len(ref_matched)
        ref_matched.append(ref_loop.rows[ref_index])
        mov_matched.append(mov_loop.rows[mov_index])
        used_mov.add(mov_index)
        key_string = _key_to_string(key)
        match_table.append(
            {
                "particle_key": key_string,
                "match_key": key_string,
                "ref_source_row_id": ref_loop.rows[ref_index].row_number,
                "mov_source_row_id": mov_loop.rows[mov_index].row_number,
                "ref_output_row_id": output_row,
                "mov_output_row_id": output_row,
                "match_status": "matched",
            }
        )

    mov_only = [
        mov_loop.rows[index]
        for index in mov_candidate_indices
        if index not in used_mov
    ]

    return _MatchRows(
        ref_matched=tuple(ref_matched),
        mov_matched=tuple(mov_matched),
        ref_only=tuple(ref_only),
        mov_only=tuple(mov_only),
        duplicate_ref=tuple(ref_duplicate_rows),
        duplicate_mov=tuple(mov_duplicate_rows),
        match_rows=tuple(match_table),
        duplicate_ref_key_count=len(ref_duplicate_keys),
        duplicate_mov_key_count=len(mov_duplicate_keys),
        warnings=tuple(warnings),
    )


def _first_candidate_indices(
    rows: Sequence[StarRow],
    keys: Sequence[tuple[str, ...]],
) -> tuple[list[int], list[StarRow]]:
    seen: set[tuple[str, ...]] = set()
    candidates: list[int] = []
    duplicates: list[StarRow] = []
    for index, key in enumerate(keys):
        if key in seen:
            duplicates.append(rows[index])
        else:
            seen.add(key)
            candidates.append(index)
    return candidates, duplicates


def _build_key(
    row: StarRow,
    loop: StarLoop,
    key_columns: Sequence[str],
    float_tolerances: Mapping[str, float],
    path_mode: str,
) -> tuple[str, ...]:
    values: list[str] = []
    column_index = loop.canonical_column_index
    for requested in key_columns:
        canonical = canonical_column_name(requested)
        index = column_index[canonical]
        real_column = loop.columns[index]
        values.append(_normalize_key_value(row.tokens[index], real_column, float_tolerances, path_mode))
    return tuple(values)


def _normalize_key_value(
    raw_value: str,
    column: str,
    float_tolerances: Mapping[str, float],
    path_mode: str,
) -> str:
    canonical = canonical_column_name(column)
    value = _unquote_simple(raw_value)
    if canonical in float_tolerances:
        try:
            numeric = float(value)
        except ValueError as exc:
            raise ValueError(f"Column {column} was given --float-tol but value is not numeric: {raw_value}") from exc
        tolerance = float_tolerances[canonical]
        return f"floatbin:{round(numeric / tolerance):d}"
    if path_mode != "exact" and _looks_like_path_column(column):
        if "@" in value:
            left, right = value.split("@", 1)
            return f"{left}@{_normalize_path(right, path_mode)}"
        return _normalize_path(value, path_mode)
    return value


def _looks_like_path_column(column: str) -> bool:
    canonical = canonical_column_name(column)
    return any(token in canonical for token in ("image", "micrograph", "movie", "path", "filename"))


def _normalize_path(value: str, path_mode: str) -> str:
    normalized = value.replace("\\", "/")
    parts = [part for part in normalized.split("/") if part not in {"", "."}]
    if path_mode == "basename":
        return parts[-1] if parts else normalized
    if path_mode.startswith("suffix:"):
        count = int(path_mode.split(":", 1)[1])
        return "/".join(parts[-count:]) if parts else normalized
    return normalized


def _unquote_simple(value: str) -> str:
    stripped = value.strip()
    if len(stripped) >= 2 and stripped[0] == stripped[-1] and stripped[0] in {"'", '"'}:
        return stripped[1:-1]
    return stripped


def _key_to_string(key: tuple[str, ...]) -> str:
    return " | ".join(key)


def _replace_loop_rows(doc: StarDocument, loop: StarLoop, rows: Sequence[StarRow]) -> list[str]:
    return list(doc.lines[: loop.row_start_line]) + [row.line for row in rows] + list(doc.lines[loop.row_end_line :])


def _output_paths(output_dir: Path) -> dict[str, Path]:
    return {
        "aligned_ref": output_dir / "aligned_ref.star",
        "aligned_mov": output_dir / "aligned_mov.star",
        "match_table": output_dir / "match_table.csv",
        "align_report": output_dir / "align_report.json",
        "ref_only": output_dir / "ref_only.star",
        "mov_only": output_dir / "mov_only.star",
        "duplicate_ref": output_dir / "duplicate_ref.star",
        "duplicate_mov": output_dir / "duplicate_mov.star",
    }


def _write_lines(path: Path, lines: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_match_table(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    fieldnames = [
        "particle_key",
        "match_key",
        "ref_source_row_id",
        "mov_source_row_id",
        "ref_output_row_id",
        "mov_output_row_id",
        "match_status",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
