"""Suggest an explicit ``cryorole align`` strategy when default STAR matching fails.

The diagnosis is side-effect-free and only *suggests*: it never runs ``align``,
never changes the preflight verdict, and ``align`` repeats every check on all
particles before writing anything.

Overlap estimates sample rows from **one** input (the reference; deterministic
seed) and query them against a full index of the **other** input, so the
estimate is unbiased. Sampling both files independently would find only about
f² of the true matches for a sampling fraction f.
"""

from __future__ import annotations

import shlex
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from cryorole.align.keys import KeySpec, build_keys, normalize_path
from cryorole.align.relion_jobs import (
    extraction_parameters,
    find_job_command,
    subtraction_parameters,
)
from cryorole.io.readers.star_reader import (
    index_star_particles,
    read_star_particle_headers,
    unquote_star_value,
)

DEFAULT_SAMPLE_SIZE = 20_000
DEFAULT_SEED = 0
RECOMMEND_MIN_OVERLAP = 0.5


def stale_subtraction_warnings(paths: Sequence[str | Path]) -> list[dict[str, Any]]:
    """Warn for inputs that sit in a RELION Subtract job recentred with ``--center``."""

    findings = []
    for path in paths:
        star = Path(path)
        command = find_job_command(star.parent)
        params = subtraction_parameters(command) if command is not None else None
        if not params or (params["center_px"] is None and not params["recenter_on_mask"]):
            continue
        findings.append({
            "code": "STALE_SUBTRACTION_COORDINATES",
            "path": str(star),
            "source": params["source"],
            "message": (
                f"{star} is the output of a RELION particle subtraction recentred on "
                f"{params['center_px'] or 'the mask centre'} ({params['source']}). RELION does not update "
                "_rlnCoordinateX/Y after recentring, so the stored coordinates are off by the projected shift. "
                "Do not match on its coordinates; write a corrected copy with "
                f"`cryorole align --fix-subtract-coordinates {shlex.quote(str(star.parent))}`."
            ),
        })
    return findings


def diagnose_star_matching(
    ref: str | Path,
    mov: str | Path,
    *,
    sample_size: int = DEFAULT_SAMPLE_SIZE,
    seed: int = DEFAULT_SEED,
) -> dict[str, Any]:
    """Evaluate the candidate ladder and suggest one explicit ``align`` command."""

    ref_path, mov_path = Path(ref), Path(mov)
    ref_headers = set(read_star_particle_headers(ref_path))
    mov_headers = set(read_star_particle_headers(mov_path))
    ref_count = index_star_particles(ref_path).row_count
    mov_count = index_star_particles(mov_path).row_count
    rng = np.random.default_rng(seed)
    sample = np.sort(rng.choice(ref_count, size=min(sample_size, ref_count), replace=False)) if ref_count else np.empty(0, int)
    candidates: list[dict[str, Any]] = []
    base = ["cryorole", "align", "--ref", str(ref_path), "--mov", str(mov_path)]

    def evaluate(name: str, pairs: Sequence[tuple[str, str]], *, path_mode: str = "exact", tolerances=None,
                 extra_args: Sequence[str] = (), note: str | None = None) -> None:
        if not all(a in ref_headers for a, _ in pairs) or not all(b in mov_headers for _, b in pairs):
            return
        spec = KeySpec(ref_columns=tuple(a for a, _ in pairs), mov_columns=tuple(b for _, b in pairs),
                       float_tolerances=tolerances or {}, path_mode=path_mode)
        ref_index = index_star_particles(ref_path, columns=spec.ref_columns)
        mov_index = index_star_particles(mov_path, columns=spec.mov_columns)
        ref_keys = build_keys(ref_index.columns, spec.ref_columns, spec)
        mov_keys = build_keys(mov_index.columns, spec.mov_columns, spec)
        ref_counts = Counter(ref_keys)
        mov_counts = Counter(mov_keys)
        matched = one_to_many = many_to_one = 0
        for row in sample:
            key = ref_keys[row]
            count = mov_counts.get(key, 0)
            if count == 0:
                continue
            if ref_counts[key] > 1:
                many_to_one += 1
            elif count > 1:
                one_to_many += 1
            else:
                matched += 1
        scale = ref_count / len(sample) if len(sample) else 0.0
        estimate = matched * scale
        record: dict[str, Any] = {
            "name": name,
            "key": spec.describe(),
            "path_mode": path_mode,
            "sampled_ref_rows": int(len(sample)),
            "estimated_matched": int(round(estimate)),
            "ref_overlap_estimate": estimate / ref_count if ref_count else 0.0,
            "mov_overlap_estimate": estimate / mov_count if mov_count else 0.0,
            "one_to_many_estimate": int(round(one_to_many * scale)),
            "many_to_one_estimate": int(round(many_to_one * scale)),
            "command": shlex.join([*base, *extra_args]),
        }
        if path_mode != "exact":
            record["stack_merge"] = _stack_merge_count(ref_index.columns[spec.ref_columns[0]], path_mode) + \
                _stack_merge_count(mov_index.columns[spec.mov_columns[0]], path_mode)
            if record["stack_merge"]:
                record["warning"] = (
                    f"{record['stack_merge']} distinct stack paths collapse to the same {path_mode} name; "
                    "unrelated stacks could be merged."
                )
        if note:
            record["note"] = note
        candidates.append(record)

    for column in ("_rlnTomoParticleName", "_rlnImageName"):
        evaluate(f"exact {column}", [(column, column)], extra_args=["--key", column])
    for mode in ("suffix:2", "basename"):
        evaluate(f"_rlnImageName with --path-mode {mode}", [("_rlnImageName", "_rlnImageName")], path_mode=mode,
                 extra_args=["--key", "_rlnImageName", "--path-mode", mode])
    for a, b in (("_rlnImageName", "_rlnImageOriginalName"), ("_rlnImageOriginalName", "_rlnImageName"),
                 ("_rlnImageOriginalName", "_rlnImageOriginalName")):
        evaluate(f"{a} = {b}", [(a, b)], extra_args=["--key-pair", f"{a}={b}"],
                 note="_rlnImageOriginalName records only the previous processing step.")
    evaluate(
        "micrograph + unchanged coordinates",
        [("_rlnMicrographName", "_rlnMicrographName"), ("_rlnCoordinateX", "_rlnCoordinateX"),
         ("_rlnCoordinateY", "_rlnCoordinateY")],
        tolerances={"rlncoordinatex": 1e-3, "rlncoordinatey": 1e-3},
        extra_args=["--coordinate-match", "unchanged"],
    )
    exact = _recentered_exact_candidate(ref_path, mov_path)
    if exact is not None:
        candidates.append(exact)

    recommended = next(
        (c for c in candidates
         if min(c["ref_overlap_estimate"], c["mov_overlap_estimate"]) >= RECOMMEND_MIN_OVERLAP
         and not c.get("one_to_many_estimate") and not c.get("many_to_one_estimate") and not c.get("stack_merge")),
        None,
    )
    return {
        "artifact_type": "cryorole_align_diagnosis",
        "schema_version": "1",
        "sampling": {"side": "ref", "size": int(len(sample)), "seed": seed, "queried_against": "full mov index"},
        "ref_row_count": ref_count,
        "mov_row_count": mov_count,
        "candidates": candidates,
        "recommended": recommended["name"] if recommended else None,
        "recommended_command": recommended["command"] if recommended else None,
        "note": "Suggestions only. `cryorole align` repeats every check on all particles before writing.",
    }


def _stack_merge_count(values: Sequence[str], path_mode: str) -> int:
    stacks: dict[str, set[str]] = {}
    for value in values:
        text = unquote_star_value(value)
        stack = text.split("@", 1)[1] if "@" in text else text
        stacks.setdefault(normalize_path(stack, path_mode), set()).add(stack)
    return sum(len(paths) - 1 for paths in stacks.values() if len(paths) > 1)


def _recentered_exact_candidate(ref_path: Path, mov_path: Path) -> dict[str, Any] | None:
    """Evaluate recentered-exact when a re-extraction note.txt sits next to either input."""

    from cryorole.align.star_align import _exact_extraction_link, _exact_summary

    for input_path, output_path in ((ref_path, mov_path), (mov_path, ref_path)):
        command = find_job_command(output_path.parent)
        params = extraction_parameters(command) if command is not None else None
        if params is None:
            continue
        try:
            in_index, out_index, link = _exact_extraction_link(
                input_path, output_path, recenter_shift=params["recenter_shift_px"],
                ref_angpix=params["ref_angpix"], micrograph_angpix=None,
            )
        except (ValueError, KeyError) as exc:
            return {"name": "recentered-exact (from note.txt)", "key": ["recentered-exact"], "error": str(exc),
                    "estimated_matched": 0, "ref_overlap_estimate": 0.0, "mov_overlap_estimate": 0.0,
                    "source": params["source"], "command": None}
        summary = _exact_summary(link, in_index.row_count, out_index.row_count)
        args = ["cryorole", "align", "--ref", str(input_path), "--mov", str(output_path),
                "--coordinate-match", "recentered-exact",
                "--recenter-shift", *[f"{v:g}" for v in params["recenter_shift_px"]]]
        if params["ref_angpix"]:
            args += ["--ref-angpix", f"{params['ref_angpix']:g}"]
        verified = summary["verified_count"]
        return {
            "name": "recentered-exact (from note.txt)",
            "key": ["recentered-exact"],
            "source": params["source"],
            "recenter_shift_px": params["recenter_shift_px"],
            "ref_angpix": params["ref_angpix"],
            "estimated_matched": verified,
            "ref_overlap_estimate": verified / _rows_of(in_index, out_index, input_path, ref_path),
            "mov_overlap_estimate": verified / _rows_of(in_index, out_index, input_path, mov_path),
            "verified_fraction": summary["verified_fraction_of_smaller_input"],
            "indistinguishable_duplicates": summary["indistinguishable_input_duplicates"],
            "swapped": input_path == mov_path,
            "command": shlex.join(args),
            "note": "Exact per-row verification of coordinates, origins, angles and one-to-one assignment.",
        }
    return None


def _rows_of(in_index, out_index, input_path: Path, which: Path) -> int:
    index = in_index if which == input_path else out_index
    return max(index.row_count, 1)
