"""Key building and exact one-to-one matching for ``cryorole align``."""

from __future__ import annotations

import itertools
import re
from collections import Counter
from dataclasses import dataclass, field
from typing import Mapping, Sequence

from cryorole.io.readers.star_reader import unquote_star_value

AUTO_KEY_CANDIDATES = (
    ("_rlnTomoParticleName",),
    ("_rlnImageName",),
)
NON_DEFAULT_KEY_MARKERS = ("defocus", "coordinate", "coord", "angle", "origin", "shift", "ctf")
PATH_LIKE_TOKENS = ("image", "micrograph", "movie", "path", "filename")


def canonical_column_name(name: str) -> str:
    return name.strip().lstrip("_").lower()


def resolve_column(headers: Sequence[str], requested: str, *, path: object, label: str) -> str:
    """Map a requested column (case-insensitive, ``_`` optional) to the real header."""

    index = {}
    for header in headers:
        index.setdefault(canonical_column_name(header), header)
    found = index.get(canonical_column_name(requested))
    if found is None:
        raise ValueError(
            f"Missing key columns in {path} ({label}): ['{requested}']. "
            f"Available columns include: {list(headers)[:20]}"
        )
    return found


@dataclass(frozen=True)
class KeySpec:
    """One composite key: parallel column lists for ref and mov."""

    ref_columns: tuple[str, ...]
    mov_columns: tuple[str, ...]
    float_tolerances: Mapping[str, float] = field(default_factory=dict)  # canonical name → tolerance
    path_mode: str = "exact"

    @property
    def is_cross_column(self) -> bool:
        return any(canonical_column_name(a) != canonical_column_name(b) for a, b in zip(self.ref_columns, self.mov_columns))

    def float_positions(self) -> tuple[int, ...]:
        return tuple(
            index
            for index, (a, b) in enumerate(zip(self.ref_columns, self.mov_columns))
            if canonical_column_name(a) in self.float_tolerances or canonical_column_name(b) in self.float_tolerances
        )

    def describe(self) -> list[str]:
        return [a if canonical_column_name(a) == canonical_column_name(b) else f"{a}={b}"
                for a, b in zip(self.ref_columns, self.mov_columns)]


def validate_path_mode(path_mode: str) -> None:
    if path_mode in {"exact", "basename"} or re.fullmatch(r"suffix:[1-9][0-9]*", path_mode):
        return
    raise ValueError("--path-mode must be exact, basename, or suffix:N with N > 0")


def normalize_float_tolerances(float_tolerances: Mapping[str, float]) -> dict[str, float]:
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


def _looks_like_path_column(column: str) -> bool:
    canonical = canonical_column_name(column)
    return any(token in canonical for token in PATH_LIKE_TOKENS)


def normalize_path(value: str, path_mode: str) -> str:
    normalized = value.replace("\\", "/")
    parts = [part for part in normalized.split("/") if part not in {"", "."}]
    if path_mode == "basename":
        return parts[-1] if parts else normalized
    if path_mode.startswith("suffix:"):
        count = int(path_mode.split(":", 1)[1])
        return "/".join(parts[-count:]) if parts else normalized
    return normalized


def normalize_key_value(raw: str, column: str, spec: KeySpec) -> object:
    """Normalize one key field. Float-tolerance columns become integer bins."""

    value = unquote_star_value(raw)
    canonical = canonical_column_name(column)
    if canonical in spec.float_tolerances:
        try:
            numeric = float(value)
        except ValueError as exc:
            raise ValueError(f"Column {column} was given --float-tol but value is not numeric: {raw}") from exc
        return int(round(numeric / spec.float_tolerances[canonical]))
    if spec.path_mode != "exact" and _looks_like_path_column(column):
        if "@" in value:
            left, right = value.split("@", 1)
            return f"{left}@{normalize_path(right, spec.path_mode)}"
        return normalize_path(value, spec.path_mode)
    return value


def build_keys(columns: Mapping[str, Sequence[str]], names: Sequence[str], spec: KeySpec) -> list[tuple]:
    value_lists = [columns[name] for name in names]
    return [
        tuple(normalize_key_value(values[i], name, spec) for values, name in zip(value_lists, names))
        for i in range(len(value_lists[0]) if value_lists else 0)
    ]


@dataclass
class MatchResult:
    pairs: list[tuple[int, int]]  # (ref_row, mov_row), in ref order
    ref_only: list[int]
    mov_only: list[int]
    duplicate_ref: list[int]
    duplicate_mov: list[int]
    ambiguous_ref: list[int]
    ambiguous_mov: list[int]
    duplicate_ref_key_count: int
    duplicate_mov_key_count: int
    warnings: list[str]
    pair_keys: list[tuple]


def match_keys(
    ref_keys: Sequence[tuple],
    mov_keys: Sequence[tuple],
    *,
    duplicate_policy: str = "exclude",
    float_positions: Sequence[int] = (),
) -> MatchResult:
    """Exact one-to-one matching with duplicate and tolerance-edge ambiguity handling.

    ``float_positions`` are key components that are tolerance bins. A pair is
    *ambiguous* if a neighbouring bin (±1 in any float component) is occupied
    in either file: the particle could belong to a different partner, so it is
    excluded rather than paired.
    """

    if duplicate_policy not in {"exclude", "first"}:
        raise ValueError("--duplicate-policy must be exclude or first")
    ref_counts = Counter(ref_keys)
    mov_counts = Counter(mov_keys)
    ref_dup_keys = {key for key, count in ref_counts.items() if count > 1}
    mov_dup_keys = {key for key, count in mov_counts.items() if count > 1}
    warnings: list[str] = []
    duplicate_ref: list[int] = []
    duplicate_mov: list[int] = []
    if duplicate_policy == "exclude":
        excluded = ref_dup_keys | mov_dup_keys
        if excluded:
            warnings.append(f"Duplicate keys detected; duplicate-policy {duplicate_policy} was applied.")
        ref_candidates = [i for i, key in enumerate(ref_keys) if key not in excluded]
        mov_candidates = [i for i, key in enumerate(mov_keys) if key not in excluded]
        duplicate_ref = [i for i, key in enumerate(ref_keys) if key in excluded]
        duplicate_mov = [i for i, key in enumerate(mov_keys) if key in excluded]
    else:
        if ref_dup_keys or mov_dup_keys:
            warnings.append(f"Duplicate keys detected; duplicate-policy {duplicate_policy} was applied.")
        warnings.append("duplicate-policy first kept the first occurrence and excluded later duplicates.")
        ref_candidates, duplicate_ref = _first_occurrences(ref_keys)
        mov_candidates, duplicate_mov = _first_occurrences(mov_keys)

    mov_by_key = {mov_keys[i]: i for i in mov_candidates}
    occupied_ref = set(ref_keys)
    occupied_mov = set(mov_keys)
    pairs: list[tuple[int, int]] = []
    pair_keys: list[tuple] = []
    ref_only: list[int] = []
    ambiguous_ref: list[int] = []
    ambiguous_mov: list[int] = []
    used_mov: set[int] = set()
    for ref_index in ref_candidates:
        key = ref_keys[ref_index]
        mov_index = mov_by_key.get(key)
        if mov_index is None or mov_index in used_mov:
            ref_only.append(ref_index)
            continue
        if float_positions and _has_occupied_neighbour(key, float_positions, occupied_ref, occupied_mov):
            ambiguous_ref.append(ref_index)
            ambiguous_mov.append(mov_index)
            used_mov.add(mov_index)
            continue
        pairs.append((ref_index, mov_index))
        pair_keys.append(key)
        used_mov.add(mov_index)
    mov_only = [i for i in mov_candidates if i not in used_mov]
    if ambiguous_ref:
        warnings.append(
            f"{len(ambiguous_ref)} matched particle(s) had another candidate within one tolerance step and were "
            "excluded as ambiguous (see ambiguous_ref.star / ambiguous_mov.star)."
        )
    return MatchResult(
        pairs=pairs,
        ref_only=ref_only,
        mov_only=mov_only,
        duplicate_ref=duplicate_ref,
        duplicate_mov=duplicate_mov,
        ambiguous_ref=ambiguous_ref,
        ambiguous_mov=ambiguous_mov,
        duplicate_ref_key_count=len(ref_dup_keys),
        duplicate_mov_key_count=len(mov_dup_keys),
        warnings=warnings,
        pair_keys=pair_keys,
    )


def _first_occurrences(keys: Sequence[tuple]) -> tuple[list[int], list[int]]:
    seen: set[tuple] = set()
    candidates: list[int] = []
    duplicates: list[int] = []
    for index, key in enumerate(keys):
        if key in seen:
            duplicates.append(index)
        else:
            seen.add(key)
            candidates.append(index)
    return candidates, duplicates


def _has_occupied_neighbour(key: tuple, float_positions: Sequence[int], *occupied_sets: set) -> bool:
    for steps in itertools.product((-1, 0, 1), repeat=len(float_positions)):
        if not any(steps):
            continue
        neighbour = list(key)
        for position, step in zip(float_positions, steps):
            neighbour[position] = neighbour[position] + step
        candidate = tuple(neighbour)
        if any(candidate in occupied for occupied in occupied_sets):
            return True
    return False


def key_to_string(key: tuple) -> str:
    return " | ".join(str(part) for part in key)


def explicit_key_warnings(spec: KeySpec) -> list[str]:
    warnings = []
    for column in dict.fromkeys((*spec.ref_columns, *spec.mov_columns)):
        normalized = canonical_column_name(column)
        if any(marker in normalized for marker in NON_DEFAULT_KEY_MARKERS):
            warnings.append(
                f"Explicit key column {column} is a user-selected tracking key, not a cryoROLE default identity "
                "key. Coordinates change after re-extraction or subtraction with recentring; for those cases use "
                "--coordinate-match recentered-exact or a particle-name key."
            )
    return warnings
