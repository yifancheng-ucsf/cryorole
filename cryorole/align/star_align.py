"""Auditable RELION STAR pre-alignment for ``cryorole run --row-aligned``.

``cryorole align`` *establishes* particle correspondence when the default
identity keys of ``cryorole run`` are not enough. Every strategy is exact:

* ``key``: exact (normalized) key columns, optionally different columns in the
  two files (``--key-pair REF_COL=MOV_COL``);
* ``coordinate-unchanged``: micrograph + identical coordinates;
* ``recentered-exact``: RELION re-extraction input → output, verified per row
  by predicted integer coordinates, predicted residual origins, identical
  angles and one-to-one assignment (see :mod:`cryorole.align.relion_geometry`);
* ``chain``: ``ref ↔ extraction input`` (names) → ``input ↔ output``
  (recentered-exact) → ``output ↔ mov`` (names). Only rows verified on every
  link are paired.

Outputs are verbatim row subsets of the inputs plus ``match_table.csv`` and
``align_report.json``. Source files are never modified.
"""

from __future__ import annotations

import csv
import json
import shlex
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from cryorole.align.keys import (
    AUTO_KEY_CANDIDATES,
    KeySpec,
    build_keys,
    canonical_column_name,
    explicit_key_warnings,
    key_to_string,
    match_keys,
    normalize_float_tolerances,
    normalize_path,
    resolve_column,
    validate_path_mode,
)
from cryorole.align.relion_geometry import predict_reextraction
from cryorole.io.readers.star_reader import (
    StarParticleIndex,
    index_star_particles,
    read_star_particle_headers,
    unquote_star_value,
    write_star_particle_subset,
)
from cryorole.align.star_values import angle_difference_deg, float_columns, optics_per_row
from cryorole.provenance.source_identity import stream_sha256

ALIGN_REPORT_SCHEMA_VERSION = "2"
LOW_OVERLAP_FRACTION = 0.50
EXACT_ACCEPT_FRACTION = 0.99
ORIGIN_TOLERANCE_ANGST = 2e-3
ANGLE_TOLERANCE_DEG = 1e-5
COORDINATE_KEY_STEP_PX = 1e-3
DEFAULT_OUTPUT_PARENT = "cryorole_alignments"

ANGLE_COLUMNS = ("_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi")
ORIGIN_COLUMNS = ("_rlnOriginXAngst", "_rlnOriginYAngst")
COORDINATE_COLUMNS = ("_rlnCoordinateX", "_rlnCoordinateY")
MICROGRAPH_COLUMN = "_rlnMicrographName"
OPTICS_GROUP_COLUMN = "_rlnOpticsGroup"
NAME_LINK_CANDIDATES = (
    ("_rlnImageName", "_rlnImageName"),
    ("_rlnTomoParticleName", "_rlnTomoParticleName"),
    ("_rlnImageName", "_rlnImageOriginalName"),
    ("_rlnImageOriginalName", "_rlnImageName"),
    ("_rlnImageOriginalName", "_rlnImageOriginalName"),
)
COORDINATE_MATCH_MODES = ("unchanged", "recentered-exact")


# --------------------------------------------------------------------------- public entry


def align_star_files(
    *,
    ref: str | Path,
    mov: str | Path,
    align_id: str = "default",
    key_columns: Sequence[str] | None = None,
    key_pairs: Sequence[str] | None = None,
    float_tolerances: Mapping[str, float] | None = None,
    path_mode: str = "exact",
    duplicate_policy: str = "exclude",
    overwrite: bool = False,
    output_dir: str | Path | None = None,
    coordinate_match: str | None = None,
    recenter_shift: Sequence[float] | None = None,
    ref_angpix: float | None = None,
    micrograph_angpix: float | None = None,
    via_extraction: Sequence[str | Path] | None = None,
) -> dict[str, Any]:
    """Align two RELION STAR particle tables and write auditable artifacts."""

    ref_path, mov_path = Path(ref), Path(mov)
    _validate_star_input(ref_path, "ref")
    _validate_star_input(mov_path, "mov")
    _validate_align_id(align_id)
    validate_path_mode(path_mode)
    if duplicate_policy not in {"exclude", "first"}:
        raise ValueError("--duplicate-policy must be exclude or first")
    if coordinate_match is not None and coordinate_match not in COORDINATE_MATCH_MODES:
        raise ValueError(f"--coordinate-match must be one of {COORDINATE_MATCH_MODES}")
    if via_extraction is not None:
        if len(via_extraction) != 2:
            raise ValueError("--via-extraction needs two STAR files: the extraction INPUT and OUTPUT")
        if coordinate_match not in (None, "recentered-exact"):
            raise ValueError("--via-extraction always uses recentered-exact matching for its middle link")
        coordinate_match = "recentered-exact"
    if coordinate_match is not None and (key_columns or key_pairs):
        raise ValueError("--coordinate-match cannot be combined with --key or --key-pair")
    if coordinate_match == "recentered-exact" and recenter_shift is None:
        raise ValueError(
            "--coordinate-match recentered-exact needs --recenter-shift X Y Z (the RELION --recenter_x/y/z "
            "values, in reference pixels). `cryorole preflight` shows them when the job's note.txt is found."
        )
    if recenter_shift is not None and coordinate_match != "recentered-exact":
        raise ValueError("--recenter-shift is only used with --coordinate-match recentered-exact")
    tolerances = normalize_float_tolerances(float_tolerances or {})

    resolved_output = Path(output_dir) if output_dir is not None else ref_path.parent / DEFAULT_OUTPUT_PARENT / align_id
    _prepare_output_dir(resolved_output, overwrite=overwrite)

    warnings: list[str] = []
    if via_extraction is not None:
        plan = _align_chain(
            ref_path=ref_path,
            mov_path=mov_path,
            extraction_input=Path(via_extraction[0]),
            extraction_output=Path(via_extraction[1]),
            recenter_shift=recenter_shift,
            ref_angpix=ref_angpix,
            micrograph_angpix=micrograph_angpix,
            duplicate_policy=duplicate_policy,
        )
    elif coordinate_match == "recentered-exact":
        plan = _align_recentered_exact(
            ref_path=ref_path,
            mov_path=mov_path,
            recenter_shift=recenter_shift,
            ref_angpix=ref_angpix,
            micrograph_angpix=micrograph_angpix,
            path_mode=path_mode,
        )
    else:
        plan = _align_by_keys(
            ref_path=ref_path,
            mov_path=mov_path,
            key_columns=key_columns,
            key_pairs=key_pairs,
            tolerances=tolerances,
            path_mode=path_mode,
            duplicate_policy=duplicate_policy,
            coordinate_unchanged=coordinate_match == "unchanged",
        )
    warnings.extend(plan.warnings)

    matched_count = len(plan.pairs)
    ref_rows, mov_rows = plan.ref_index.row_count, plan.mov_index.row_count
    if matched_count == 0:
        _cleanup_empty_output_dir(resolved_output)
        raise ValueError(plan.zero_overlap_message)
    ref_overlap = matched_count / ref_rows if ref_rows else 0.0
    mov_overlap = matched_count / mov_rows if mov_rows else 0.0
    low_overlap_warning = 0.0 < min(ref_overlap, mov_overlap) < LOW_OVERLAP_FRACTION
    if low_overlap_warning:
        warnings.append("Low overlap warning: matched rows are below 50% of at least one input.")

    paths = _output_paths(resolved_output)
    comment = f"cryoROLE align ({plan.strategy}); report: align_report.json"
    ref_pairs = [pair[0] for pair in plan.pairs]
    mov_pairs = [pair[1] for pair in plan.pairs]
    write_star_particle_subset(plan.ref_index, ref_pairs, paths["aligned_ref"], header_comment=comment)
    write_star_particle_subset(plan.mov_index, mov_pairs, paths["aligned_mov"], header_comment=comment)
    for name, index, rows in (
        ("ref_only", plan.ref_index, plan.ref_only),
        ("mov_only", plan.mov_index, plan.mov_only),
        ("duplicate_ref", plan.ref_index, plan.duplicate_ref),
        ("duplicate_mov", plan.mov_index, plan.duplicate_mov),
        ("ambiguous_ref", plan.ref_index, plan.ambiguous_ref),
        ("ambiguous_mov", plan.mov_index, plan.ambiguous_mov),
    ):
        write_star_particle_subset(index, rows, paths[name], header_comment=comment)
    if plan.unverified_ref is not None:
        paths["unverified_ref"] = resolved_output / "unverified_ref.star"
        paths["unverified_mov"] = resolved_output / "unverified_mov.star"
        write_star_particle_subset(plan.ref_index, plan.unverified_ref, paths["unverified_ref"], header_comment=comment)
        write_star_particle_subset(plan.mov_index, plan.unverified_mov, paths["unverified_mov"], header_comment=comment)
    _write_match_table(paths["match_table"], plan)

    hashes = {
        "ref_source": stream_sha256(ref_path),
        "mov_source": stream_sha256(mov_path),
        "aligned_ref": stream_sha256(paths["aligned_ref"]),
        "aligned_mov": stream_sha256(paths["aligned_mov"]),
        "match_table": stream_sha256(paths["match_table"]),
    }
    next_command = shlex.join([
        "cryorole", "run",
        "--ref", str(paths["aligned_ref"].resolve()),
        "--mov", str(paths["aligned_mov"].resolve()),
        "--row-aligned",
    ])
    report: dict[str, Any] = {
        "artifact_type": "align_report",
        "schema_version": ALIGN_REPORT_SCHEMA_VERSION,
        "status": "ok",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "strategy": plan.strategy,
        "ref_path": str(ref_path),
        "mov_path": str(mov_path),
        "ref_resolved_path": str(ref_path.resolve()),
        "mov_resolved_path": str(mov_path.resolve()),
        "align_id": align_id,
        "output_dir": str(resolved_output),
        "output_paths": {key: str(path) for key, path in paths.items()},
        "key_policy": plan.key_policy,
        "key_columns": plan.key_policy.get("key_columns", []),
        "key_source": plan.key_policy.get("key_source"),
        "path_mode": path_mode,
        "float_tolerances": dict(tolerances),
        "duplicate_policy": duplicate_policy,
        "ref_row_count": ref_rows,
        "mov_row_count": mov_rows,
        "matched_count": matched_count,
        "aligned_row_count": matched_count,
        "ref_only_count": len(plan.ref_only),
        "mov_only_count": len(plan.mov_only),
        "duplicate_ref_row_count": len(plan.duplicate_ref),
        "duplicate_mov_row_count": len(plan.duplicate_mov),
        "duplicate_ref_key_count": plan.duplicate_ref_key_count,
        "duplicate_mov_key_count": plan.duplicate_mov_key_count,
        "ambiguous_ref_row_count": len(plan.ambiguous_ref),
        "ambiguous_mov_row_count": len(plan.ambiguous_mov),
        "ref_overlap_fraction": ref_overlap,
        "mov_overlap_fraction": mov_overlap,
        "low_overlap_warning": low_overlap_warning,
        "sha256": hashes,
        "next_command": next_command,
        "warnings": warnings,
    }
    if plan.details:
        report.update(plan.details)
    report["coverage"] = _coverage(matched_count, ref_rows, mov_rows)
    if plan.ambiguous_groups:
        groups_path = resolved_output / "ambiguous_groups.csv"
        report["suspected_duplicate_groups"] = _write_ambiguous_groups(groups_path, plan)
        report["output_paths"]["ambiguous_groups"] = str(groups_path)
    if report["coverage"]["label"] != "full" and plan.strategy in {"recentered-exact", "chain"}:
        hint = _unique_name_key_hint(ref_path, mov_path, duplicate_policy=duplicate_policy)
        if hint:
            report["unique_key_available"] = hint
            warnings.append(
                f"A unique particle-name key ({hint['key']}) pairs all {hint['matched']} rows of both files. For a "
                f"full paired baseline run: {hint['command']}"
            )
    if report["coverage"]["label"] != "full":
        warnings.insert(0, report["coverage"]["statement"])
    _write_json(paths["align_report"], report)
    return report


def _coverage(paired: int, ref_rows: int, mov_rows: int) -> dict[str, Any]:
    full = paired == ref_rows == mov_rows
    label = "full" if full else "matchable subset"
    statement = (
        f"FULL: all {paired} rows of both files are paired."
        if full
        else (
            f"MATCHABLE SUBSET: {paired} of {ref_rows} ref rows ({paired / max(ref_rows, 1):.1%}) and {paired} of "
            f"{mov_rows} mov rows ({paired / max(mov_rows, 1):.1%}) are paired. A landscape from these files covers "
            "only this subset; it is neither the full data nor a deduplicated dataset."
        )
    )
    return {
        "label": label,
        "paired": paired,
        "ref_rows": ref_rows,
        "mov_rows": mov_rows,
        "ref_fraction": paired / max(ref_rows, 1),
        "mov_fraction": paired / max(mov_rows, 1),
        "statement": statement,
    }


def _write_ambiguous_groups(path: Path, plan: _Plan) -> dict[str, Any]:
    columns = ("_rlnImageName", "_rlnMicrographName", "_rlnCoordinateX", "_rlnCoordinateY", "_rlnRandomSubset",
               "_rlnImageOriginalName")
    sides = {}
    for side, index in (("ref", plan.ref_index), ("mov", plan.mov_index)):
        available = [c for c in columns if c in index.headers]
        sides[side] = (available, index_star_particles(index.path, columns=available).columns if available else {})
    ref_rows_total = mov_rows_total = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["group_id", "side", "source_row_id", *columns])
        for group_id, (ref_rows, mov_rows) in enumerate(plan.ambiguous_groups or []):
            for side, rows in (("ref", ref_rows), ("mov", mov_rows)):
                available, values = sides[side]
                for row in rows:
                    writer.writerow([group_id, side, row, *[values[c][row] if c in values else "" for c in columns]])
            ref_rows_total += len(ref_rows)
            mov_rows_total += len(mov_rows)
    groups = len(plan.ambiguous_groups or [])
    return {
        "groups": groups,
        "rows_excluded_from_geometric_matching": {"ref": ref_rows_total, "mov": mov_rows_total},
        "rows_in_excess_if_each_group_is_one_particle": {"ref": ref_rows_total - groups, "mov": mov_rows_total - groups},
        "evidence": "identical micrograph, predicted coordinates, residual origins (0.002 A) and angles (1e-5 deg) "
                    "within each group, so geometry cannot tell the members apart",
        "interpretation": "suspected duplicates; not confirmed. cryoROLE never removes particles and never guesses a "
                          "pairing inside a group.",
        "next_steps": [
            "Trace the picking and extraction history (joined picking runs, overlapping picks, half-set "
            "assignments) before deciding whether groups are one physical particle.",
            "If a unique particle-name key pairs all rows, align with --key-pair for a full paired baseline.",
            "Any deduplication is a separate, explicit decision; cryoROLE does not perform it.",
        ],
        "file": str(path),
    }


def _unique_name_key_hint(ref_path: Path, mov_path: Path, *, duplicate_policy: str) -> dict[str, Any] | None:
    try:
        mapping, info = _name_link(ref_path, mov_path, duplicate_policy="exclude")
    except ValueError:
        return None
    ref_rows = info["input_rows"][0]
    mov_rows = info["input_rows"][1]
    if not (len(mapping) == ref_rows == mov_rows and len(set(mapping.values())) == len(mapping)):
        return None
    key = info["key"]
    flag = ["--key-pair", key] if "=" in key else ["--key", key]
    return {
        "key": key,
        "matched": len(mapping),
        "command": shlex.join(["cryorole", "align", "--ref", str(ref_path), "--mov", str(mov_path), *flag,
                               "--align-id", "full_baseline"]),
    }


# --------------------------------------------------------------------------- plans


@dataclass
class _Plan:
    strategy: str
    ref_index: StarParticleIndex
    mov_index: StarParticleIndex
    pairs: list[tuple[int, int]]
    pair_keys: list[str]
    ref_only: list[int]
    mov_only: list[int]
    duplicate_ref: list[int]
    duplicate_mov: list[int]
    ambiguous_ref: list[int]
    ambiguous_mov: list[int]
    duplicate_ref_key_count: int
    duplicate_mov_key_count: int
    warnings: list[str]
    key_policy: dict[str, Any]
    zero_overlap_message: str
    unverified_ref: list[int] | None = None
    unverified_mov: list[int] | None = None
    extra_columns: dict[str, list[Any]] | None = None
    details: dict[str, Any] | None = None
    ambiguous_groups: list[tuple[list[int], list[int]]] | None = None


def _align_by_keys(
    *,
    ref_path: Path,
    mov_path: Path,
    key_columns: Sequence[str] | None,
    key_pairs: Sequence[str] | None,
    tolerances: Mapping[str, float],
    path_mode: str,
    duplicate_policy: str,
    coordinate_unchanged: bool,
) -> _Plan:
    ref_headers = read_star_particle_headers(ref_path)
    mov_headers = read_star_particle_headers(mov_path)
    key_source = "explicit"
    if coordinate_unchanged:
        key_source = "coordinate_unchanged"
        pairs = [(MICROGRAPH_COLUMN, MICROGRAPH_COLUMN), *((c, c) for c in COORDINATE_COLUMNS)]
        tolerances = {**tolerances, **{canonical_column_name(c): COORDINATE_KEY_STEP_PX for c in COORDINATE_COLUMNS}}
    elif key_pairs:
        pairs = [_parse_key_pair(item) for item in key_pairs]
        pairs = [*((c, c) for c in (key_columns or ())), *pairs]
    elif key_columns:
        pairs = [(c, c) for c in key_columns]
    else:
        key_source = "auto"
        pairs = _resolve_auto_key(ref_headers, mov_headers)
    ref_columns = tuple(resolve_column(ref_headers, a, path=ref_path, label="ref") for a, _ in pairs)
    mov_columns = tuple(resolve_column(mov_headers, b, path=mov_path, label="mov") for _, b in pairs)
    spec = KeySpec(ref_columns=ref_columns, mov_columns=mov_columns, float_tolerances=tolerances, path_mode=path_mode)
    ref_index = index_star_particles(ref_path, columns=ref_columns)
    mov_index = index_star_particles(mov_path, columns=mov_columns)
    ref_keys = build_keys(ref_index.columns, ref_columns, spec)
    mov_keys = build_keys(mov_index.columns, mov_columns, spec)
    match = match_keys(ref_keys, mov_keys, duplicate_policy=duplicate_policy, float_positions=spec.float_positions())
    warnings = [] if key_source in {"auto", "coordinate_unchanged"} else explicit_key_warnings(spec)
    warnings.extend(match.warnings)
    strategy = "coordinate-unchanged" if coordinate_unchanged else ("key-pair" if spec.is_cross_column else "key")
    return _Plan(
        strategy=strategy,
        ref_index=ref_index,
        mov_index=mov_index,
        pairs=match.pairs,
        pair_keys=[key_to_string(key) for key in match.pair_keys],
        ref_only=match.ref_only,
        mov_only=match.mov_only,
        duplicate_ref=match.duplicate_ref,
        duplicate_mov=match.duplicate_mov,
        ambiguous_ref=match.ambiguous_ref,
        ambiguous_mov=match.ambiguous_mov,
        duplicate_ref_key_count=match.duplicate_ref_key_count,
        duplicate_mov_key_count=match.duplicate_mov_key_count,
        warnings=warnings,
        key_policy={
            "key_source": key_source,
            "key_columns": spec.describe(),
            "ref_key_columns": list(ref_columns),
            "mov_key_columns": list(mov_columns),
            "auto_candidates": [list(candidate) for candidate in AUTO_KEY_CANDIDATES],
        },
        zero_overlap_message=(
            "Resolved STAR key has zero overlap; specify a different --key or --key-pair, "
            "or adjust --path-mode / --float-tol. `cryorole preflight --ref … --mov …` lists candidate keys."
        ),
    )


def _parse_key_pair(item: str) -> tuple[str, str]:
    if "=" not in item:
        raise ValueError(f"--key-pair must be REF_COLUMN=MOV_COLUMN, got {item!r}")
    ref_column, mov_column = (part.strip() for part in item.split("=", 1))
    if not ref_column or not mov_column:
        raise ValueError(f"--key-pair must be REF_COLUMN=MOV_COLUMN, got {item!r}")
    return ref_column, mov_column


def _resolve_auto_key(ref_headers: Sequence[str], mov_headers: Sequence[str]) -> list[tuple[str, str]]:
    ref_canonical = {canonical_column_name(h) for h in ref_headers}
    mov_canonical = {canonical_column_name(h) for h in mov_headers}
    for candidate in AUTO_KEY_CANDIDATES:
        if all(canonical_column_name(c) in ref_canonical and canonical_column_name(c) in mov_canonical for c in candidate):
            return [(c, c) for c in candidate]
    raise ValueError(
        "No safe automatic STAR key was found. Specify an explicit key, for example: "
        "--key _rlnMicrographName _rlnCoordinateX _rlnCoordinateY, or a cross-column key such as "
        "--key-pair _rlnImageName=_rlnImageOriginalName."
    )


# ----- recentered-exact


@dataclass
class _ExactLink:
    pairs: list[tuple[int, int]]
    unverified_in: list[int]
    unverified_out: list[int]
    input_only: list[int]
    output_only: list[int]
    ambiguous_in: list[int]
    ambiguous_out: list[int]
    origin_error: list[float]
    reasons: dict[str, int]
    parameters: dict[str, Any]
    predicted_keys: list[str]
    groups: list[tuple[list[int], list[int]]] | None = None  # indistinguishable (input rows, output rows)




def _pixel_sizes(index: StarParticleIndex, *, micrograph_angpix: float | None) -> tuple[np.ndarray, np.ndarray, dict[str, str]]:
    sources: dict[str, str] = {}
    particle = optics_per_row(index, "_rlnImagePixelSize")
    if particle is None:
        raise ValueError(f"{index.path}: data_optics has no _rlnImagePixelSize (needed for recentring geometry)")
    sources["particle_angpix"] = "_rlnImagePixelSize (input optics)"
    if micrograph_angpix is not None:
        micrograph = np.full(index.row_count, float(micrograph_angpix))
        sources["micrograph_angpix"] = "--micrograph-angpix"
    else:
        micrograph = optics_per_row(index, "_rlnMicrographPixelSize")
        sources["micrograph_angpix"] = "_rlnMicrographPixelSize (input optics)"
        if micrograph is None:
            micrograph = optics_per_row(index, "_rlnMicrographOriginalPixelSize")
            sources["micrograph_angpix"] = (
                "_rlnMicrographOriginalPixelSize (input optics; assumes micrographs were not binned — "
                "verified by the exact origin/coordinate agreement)"
            )
        if micrograph is None:
            raise ValueError(f"{index.path}: no micrograph pixel size in data_optics; give --micrograph-angpix")
    return particle, micrograph, sources




def _exact_extraction_link(
    input_path: Path,
    output_path: Path,
    *,
    recenter_shift: Sequence[float],
    ref_angpix: float | None,
    micrograph_angpix: float | None,
    path_mode: str = "exact",
) -> tuple[StarParticleIndex, StarParticleIndex, _ExactLink]:
    needed = (MICROGRAPH_COLUMN, *COORDINATE_COLUMNS, *ANGLE_COLUMNS, *ORIGIN_COLUMNS, OPTICS_GROUP_COLUMN)
    for path, label in ((input_path, "extraction input"), (output_path, "extraction output")):
        headers = read_star_particle_headers(path)
        missing = [c for c in needed if c != OPTICS_GROUP_COLUMN and c not in headers]
        if missing:
            raise ValueError(f"{path} ({label}) lacks columns needed for recentered-exact matching: {missing}")
    in_index = index_star_particles(input_path, columns=needed)
    out_index = index_star_particles(output_path, columns=needed)
    a_part, a_mic, sources = _pixel_sizes(in_index, micrograph_angpix=micrograph_angpix)
    a_ref = np.full(in_index.row_count, float(ref_angpix)) if ref_angpix else a_part
    sources["ref_angpix"] = "--ref-angpix" if ref_angpix else "particle pixel size (RELION default without --ref_angpix)"
    prediction = predict_reextraction(
        coordinates=float_columns(in_index, COORDINATE_COLUMNS),
        origins_angst=float_columns(in_index, ORIGIN_COLUMNS),
        angles_deg=float_columns(in_index, ANGLE_COLUMNS),
        recenter_ref_px=np.asarray(recenter_shift, dtype=float),
        ref_angpix=a_ref,
        particle_angpix=a_part,
        micrograph_angpix=a_mic,
    )

    def mic_key(value: str) -> str:
        text = unquote_star_value(value)
        return normalize_path(text, path_mode) if path_mode != "exact" else text

    def coord_key(value: float) -> int:
        return int(round(value / COORDINATE_KEY_STEP_PX))

    in_mics = [mic_key(v) for v in in_index.columns[MICROGRAPH_COLUMN]]
    out_mics = [mic_key(v) for v in out_index.columns[MICROGRAPH_COLUMN]]
    out_coords = float_columns(out_index, COORDINATE_COLUMNS)
    in_keys = [(m, coord_key(x), coord_key(y)) for m, (x, y) in zip(in_mics, prediction.coordinates)]
    out_keys = [(m, coord_key(x), coord_key(y)) for m, (x, y) in zip(out_mics, out_coords)]
    in_origin = prediction.origins_angst
    out_origin = float_columns(out_index, ORIGIN_COLUMNS)
    in_angles = float_columns(in_index, ANGLE_COLUMNS)
    out_angles = float_columns(out_index, ANGLE_COLUMNS)
    out_by_key: dict[tuple, list[int]] = {}
    for j, key in enumerate(out_keys):
        out_by_key.setdefault(key, []).append(j)
    reasons = {"origin_mismatch": 0, "angle_mismatch": 0, "several_verified_candidates": 0, "output_claimed_twice": 0}
    signature_counts: dict[tuple, int] = {}
    signatures = []
    for i, key in enumerate(in_keys):
        signature = (key, *np.round(in_origin[i] / ORIGIN_TOLERANCE_ANGST).astype(int), *np.round(in_angles[i] / ANGLE_TOLERANCE_DEG).astype(int))
        signatures.append(signature)
        signature_counts[signature] = signature_counts.get(signature, 0) + 1
    indistinguishable = {i for i, signature in enumerate(signatures) if signature_counts[signature] > 1}
    group_rows: dict[tuple, list[int]] = {}
    for i in sorted(indistinguishable):
        group_rows.setdefault(signatures[i], []).append(i)
    candidate_of: dict[int, int] = {}
    unverified_in: list[int] = []
    unverified_out: list[int] = []
    input_only: list[int] = []
    ambiguous_in: list[int] = []
    for i, key in enumerate(in_keys):
        candidates = out_by_key.get(key)
        if not candidates:
            input_only.append(i)
            continue
        verified = []
        for j in candidates:
            origin_error = float(np.max(np.abs(in_origin[i] - out_origin[j])))
            angle_error = float(np.max(np.abs(angle_difference_deg(in_angles[i], out_angles[j]))))
            if origin_error <= ORIGIN_TOLERANCE_ANGST and angle_error <= ANGLE_TOLERANCE_DEG:
                verified.append(j)
        if i in indistinguishable:
            ambiguous_in.append(i)
        elif len(verified) == 1:
            candidate_of[i] = verified[0]
        elif len(verified) > 1:
            reasons["several_verified_candidates"] += 1
            ambiguous_in.append(i)
        else:
            j = candidates[0]
            origin_error = float(np.max(np.abs(in_origin[i] - out_origin[j])))
            if origin_error > ORIGIN_TOLERANCE_ANGST:
                reasons["origin_mismatch"] += 1
            else:
                reasons["angle_mismatch"] += 1
            unverified_in.append(i)
            if len(candidates) == 1:
                unverified_out.append(j)
    claims: dict[int, list[int]] = {}
    for i, j in candidate_of.items():
        claims.setdefault(j, []).append(i)
    pairs: list[tuple[int, int]] = []
    keys: list[str] = []
    errors: list[float] = []
    for i in sorted(candidate_of):
        j = candidate_of[i]
        if len(claims[j]) > 1:
            reasons["output_claimed_twice"] += 1
            ambiguous_in.append(i)
            continue
        pairs.append((i, j))
        keys.append(key_to_string(in_keys[i]))
        errors.append(float(np.max(np.abs(in_origin[i] - out_origin[j]))))
    used_out = {j for _, j in pairs} | set(unverified_out)
    ambiguous_out = sorted({j for i in ambiguous_in for j in out_by_key.get(in_keys[i], []) if j not in used_out})
    output_only = [j for j in range(out_index.row_count) if j not in used_out and j not in set(ambiguous_out)]
    reasons["indistinguishable_input_duplicates"] = len(indistinguishable)
    link = _ExactLink(
        pairs=pairs,
        unverified_in=unverified_in,
        unverified_out=unverified_out,
        input_only=input_only,
        output_only=output_only,
        ambiguous_in=sorted(ambiguous_in),
        ambiguous_out=ambiguous_out,
        origin_error=errors,
        reasons=reasons,
        parameters={
            "model": "relion_reextraction",
            "recenter_shift_ref_px": [float(v) for v in recenter_shift],
            "ref_angpix": float(ref_angpix) if ref_angpix else None,
            "pixel_size_sources": sources,
            "particle_angpix_values": sorted({float(v) for v in a_part}),
            "micrograph_angpix_values": sorted({float(v) for v in a_mic}),
            "matrix": "RELION Euler_angles2matrix(rot, tilt, psi, A, false)",
            "rounding": "RELION ROUND (half away from zero), micrograph pixels",
            "acceptance": {
                "coordinates": "predicted integer coordinates equal",
                "origin_tolerance_angst": ORIGIN_TOLERANCE_ANGST,
                "angle_tolerance_deg": ANGLE_TOLERANCE_DEG,
                "one_to_one": True,
                "identical_coordinates": "resolved by origin and angle agreement; ambiguous if not unique",
                "min_verified_fraction": EXACT_ACCEPT_FRACTION,
            },
        },
        predicted_keys=keys,
        groups=[
            (rows, [j for j in out_by_key.get(signature[0], []) if j not in used_out])
            for signature, rows in group_rows.items()
        ],
    )
    return in_index, out_index, link




def _exact_summary(link: _ExactLink, in_rows: int, out_rows: int) -> dict[str, Any]:
    verified = len(link.pairs)
    indistinguishable = link.reasons.get("indistinguishable_input_duplicates", 0)
    denominator = max(min(in_rows, out_rows) - indistinguishable, 1)
    return {
        "indistinguishable_input_duplicates": indistinguishable,
        "verified_count": verified,
        "verified_fraction_of_smaller_input": verified / denominator,
        "verified_fraction_note": "denominator excludes indistinguishable duplicates (identical micrograph, coordinates, origins and angles)",
        "input_rows": in_rows,
        "output_rows": out_rows,
        "unverified_count": len(link.unverified_in),
        "unverified_reasons": dict(link.reasons),
        "ambiguous_input_rows": len(link.ambiguous_in),
        "ambiguous_output_rows": len(link.ambiguous_out),
        "unmatched_input_rows": len(link.input_only),
        "unmatched_output_rows": len(link.output_only),
        "max_origin_error_angst": max(link.origin_error) if link.origin_error else None,
    }


def _exact_warnings(summary: Mapping[str, Any]) -> list[str]:
    count = summary.get("indistinguishable_input_duplicates", 0)
    if not count:
        return []
    return [
        f"{count} extraction-input rows are indistinguishable duplicates (identical micrograph, coordinates, origins "
        "and angles as another row), usually duplicate picks of one particle that converged to the same pose. They "
        "cannot be paired unambiguously and were excluded (see ambiguous_*.star). Consider removing duplicates in RELION."
    ]


def _require_exact_acceptance(summary: Mapping[str, Any], *, label: str) -> None:
    if summary["verified_fraction_of_smaller_input"] < EXACT_ACCEPT_FRACTION:
        raise ValueError(
            f"{label}: only {summary['verified_count']} rows ({summary['verified_fraction_of_smaller_input']:.1%} of "
            f"the smaller file) verify as a RELION re-extraction input/output pair for this recentre vector "
            f"(need ≥ {EXACT_ACCEPT_FRACTION:.0%}). These files are not an extraction input/output pair for this "
            f"vector, or the pixel sizes differ from those assumed. Unverified: {summary['unverified_reasons']}; "
            "no aligned files were written."
        )


def _align_recentered_exact(
    *,
    ref_path: Path,
    mov_path: Path,
    recenter_shift: Sequence[float],
    ref_angpix: float | None,
    micrograph_angpix: float | None,
    path_mode: str,
) -> _Plan:
    in_index, out_index, link = _exact_extraction_link(
        ref_path,
        mov_path,
        recenter_shift=recenter_shift,
        ref_angpix=ref_angpix,
        micrograph_angpix=micrograph_angpix,
        path_mode=path_mode,
    )
    summary = _exact_summary(link, in_index.row_count, out_index.row_count)
    _require_exact_acceptance(summary, label="recentered-exact")
    return _Plan(
        strategy="recentered-exact",
        ref_index=in_index,
        mov_index=out_index,
        pairs=link.pairs,
        pair_keys=link.predicted_keys,
        ref_only=link.input_only,
        mov_only=link.output_only,
        duplicate_ref=[],
        duplicate_mov=[],
        ambiguous_ref=link.ambiguous_in,
        ambiguous_mov=link.ambiguous_out,
        duplicate_ref_key_count=0,
        duplicate_mov_key_count=0,
        warnings=_exact_warnings(summary),
        key_policy={"key_source": "recentered_exact", "key_columns": [MICROGRAPH_COLUMN, *COORDINATE_COLUMNS]},
        zero_overlap_message="recentered-exact: no row verified.",
        unverified_ref=link.unverified_in,
        unverified_mov=link.unverified_out,
        extra_columns={"coord_verified": [True] * len(link.pairs), "origin_error_A": link.origin_error},
        details={"coordinate_match": {**link.parameters, **summary}},
        ambiguous_groups=link.groups,
    )


# ----- chain


def _name_link(
    a_path: Path,
    b_path: Path,
    *,
    duplicate_policy: str,
    a_index: StarParticleIndex | None = None,
    b_index: StarParticleIndex | None = None,
) -> tuple[dict[int, int], dict[str, Any]]:
    """Best exact particle-name link between two files (unique keys only)."""

    a_headers = read_star_particle_headers(a_path)
    b_headers = read_star_particle_headers(b_path)
    best: tuple[int, dict[int, int], dict[str, Any]] | None = None
    tried = []
    for a_col, b_col in NAME_LINK_CANDIDATES:
        if a_col not in a_headers or b_col not in b_headers:
            continue
        a_values = index_star_particles(a_path, columns=(a_col,)).columns[a_col]
        b_values = index_star_particles(b_path, columns=(b_col,)).columns[b_col]
        spec = KeySpec(ref_columns=(a_col,), mov_columns=(b_col,))
        match = match_keys(
            [(unquote_star_value(v),) for v in a_values],
            [(unquote_star_value(v),) for v in b_values],
            duplicate_policy=duplicate_policy,
        )
        mapping = dict(match.pairs)
        info = {
            "key": spec.describe()[0],
            "input_rows": [len(a_values), len(b_values)],
            "matched": len(mapping),
            "duplicate_rows": [len(match.duplicate_ref), len(match.duplicate_mov)],
            "unmatched_rows": [len(match.ref_only), len(match.mov_only)],
        }
        tried.append(info)
        if best is None or len(mapping) > best[0]:
            best = (len(mapping), mapping, info)
    if best is None or best[0] == 0:
        raise ValueError(
            f"No exact particle-name link between {a_path} and {b_path} (tried {[c for c in NAME_LINK_CANDIDATES]}); "
            "the chain needs name links on both sides of the extraction pair."
        )
    return best[1], {**best[2], "candidates_tried": tried}


def _align_chain(
    *,
    ref_path: Path,
    mov_path: Path,
    extraction_input: Path,
    extraction_output: Path,
    recenter_shift: Sequence[float],
    ref_angpix: float | None,
    micrograph_angpix: float | None,
    duplicate_policy: str,
) -> _Plan:
    for path, label in ((extraction_input, "extraction input"), (extraction_output, "extraction output")):
        _validate_star_input(path, label)
    ref_to_in, link1 = _name_link(ref_path, extraction_input, duplicate_policy=duplicate_policy)
    in_index, out_index, exact = _exact_extraction_link(
        extraction_input,
        extraction_output,
        recenter_shift=recenter_shift,
        ref_angpix=ref_angpix,
        micrograph_angpix=micrograph_angpix,
    )
    exact_summary = _exact_summary(exact, in_index.row_count, out_index.row_count)
    _require_exact_acceptance(exact_summary, label="--via-extraction middle link")
    in_to_out = dict(exact.pairs)
    out_to_mov, link3 = _name_link(extraction_output, mov_path, duplicate_policy=duplicate_policy)
    ref_index = index_star_particles(ref_path)
    mov_index = index_star_particles(mov_path)
    pairs: list[tuple[int, int]] = []
    keys: list[str] = []
    in_rows: list[int] = []
    out_rows: list[int] = []
    used_mov: set[int] = set()
    broken = {"ref_to_input": 0, "input_to_output": 0, "output_to_mov": 0}
    for ref_row in range(ref_index.row_count):
        in_row = ref_to_in.get(ref_row)
        if in_row is None:
            broken["ref_to_input"] += 1
            continue
        out_row = in_to_out.get(in_row)
        if out_row is None:
            broken["input_to_output"] += 1
            continue
        mov_row = out_to_mov.get(out_row)
        if mov_row is None:
            broken["output_to_mov"] += 1
            continue
        if mov_row in used_mov:
            continue
        used_mov.add(mov_row)
        pairs.append((ref_row, mov_row))
        in_rows.append(in_row)
        out_rows.append(out_row)
        keys.append(f"chain:{in_row}->{out_row}")
    paired_ref = {p[0] for p in pairs}
    return _Plan(
        strategy="chain",
        ref_index=ref_index,
        mov_index=mov_index,
        pairs=pairs,
        pair_keys=keys,
        ref_only=[i for i in range(ref_index.row_count) if i not in paired_ref],
        mov_only=[i for i in range(mov_index.row_count) if i not in used_mov],
        duplicate_ref=[],
        duplicate_mov=[],
        ambiguous_ref=[],
        ambiguous_mov=[],
        duplicate_ref_key_count=0,
        duplicate_mov_key_count=0,
        warnings=_exact_warnings(exact_summary),
        key_policy={"key_source": "chain", "key_columns": [link1["key"], "recentered-exact", link3["key"]]},
        zero_overlap_message="--via-extraction: no particle is linked on all three steps.",
        extra_columns={"extraction_input_row_id": in_rows, "extraction_output_row_id": out_rows},
        ambiguous_groups=_chain_groups(exact.groups or [], ref_to_in, out_to_mov),
        details={
            "chain": {
                "extraction_input": str(extraction_input),
                "extraction_output": str(extraction_output),
                "links": [
                    {"step": "ref -> extraction input", "method": "exact particle name", **link1},
                    {"step": "extraction input -> output", "method": "recentered-exact", **exact_summary},
                    {"step": "extraction output -> mov", "method": "exact particle name", **link3},
                ],
                "rows_lost_per_link": broken,
                "final_pair": {
                    "ref_rows": ref_index.row_count,
                    "mov_rows": mov_index.row_count,
                    "paired": len(pairs),
                    "ref_unpaired": ref_index.row_count - len(pairs),
                    "mov_unpaired": mov_index.row_count - len(pairs),
                },
                "sha256": {
                    "extraction_input": stream_sha256(extraction_input),
                    "extraction_output": stream_sha256(extraction_output),
                },
            },
            "coordinate_match": exact.parameters,
        },
    )


def _chain_groups(groups, ref_to_in: dict[int, int], out_to_mov: dict[int, int]):
    in_to_ref = {v: k for k, v in ref_to_in.items()}
    mapped = []
    for in_rows, out_rows in groups:
        mapped.append((
            [in_to_ref[i] for i in in_rows if i in in_to_ref],
            [out_to_mov[j] for j in out_rows if j in out_to_mov],
        ))
    return mapped


# --------------------------------------------------------------------------- outputs and helpers


def _validate_star_input(path: Path, label: str) -> None:
    if path.suffix.lower() != ".star":
        raise ValueError(f"cryorole align only supports STAR input for now; {label} is not .star: {path}")
    if not path.is_file():
        raise ValueError(f"{label} STAR file does not exist: {path}")


def _validate_align_id(align_id: str) -> None:
    if not align_id or Path(align_id).name != align_id or align_id in {".", ".."}:
        raise ValueError("--align-id must be a simple directory name")


def _prepare_output_dir(path: Path, *, overwrite: bool) -> None:
    if path.exists() and not path.is_dir():
        raise ValueError(f"Align output path exists and is not a directory: {path}")
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Align output directory already exists: {path}. Use --overwrite to replace it, "
            "or choose another --align-id / --output-dir."
        )
    path.mkdir(parents=True, exist_ok=True)


def _cleanup_empty_output_dir(path: Path) -> None:
    try:
        path.rmdir()
        if path.parent.name == DEFAULT_OUTPUT_PARENT:
            path.parent.rmdir()
    except OSError:
        pass


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
        "ambiguous_ref": output_dir / "ambiguous_ref.star",
        "ambiguous_mov": output_dir / "ambiguous_mov.star",
    }


MATCH_TABLE_COLUMNS = (
    "particle_key",
    "match_key",
    "ref_source_row_id",
    "mov_source_row_id",
    "ref_output_row_id",
    "mov_output_row_id",
    "match_status",
)


def _write_match_table(path: Path, plan: _Plan) -> None:
    extra = plan.extra_columns or {}
    fieldnames = [*MATCH_TABLE_COLUMNS, *extra.keys()]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for output_row, ((ref_row, mov_row), key) in enumerate(zip(plan.pairs, plan.pair_keys)):
            row = {
                "particle_key": key,
                "match_key": key,
                "ref_source_row_id": ref_row,
                "mov_source_row_id": mov_row,
                "ref_output_row_id": output_row,
                "mov_output_row_id": output_row,
                "match_status": "matched",
            }
            for name, values in extra.items():
                row[name] = values[output_row]
            writer.writerow(row)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, default=str)
        handle.write("\n")
