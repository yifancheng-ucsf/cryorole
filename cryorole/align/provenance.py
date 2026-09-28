"""Attach ``cryorole align`` lineage to an explicit ``run --row-aligned``.

``--row-aligned`` is the user's assertion and never depends on this module:
nothing here can block or change a run. When an ``align_report.json`` sits next
to the run inputs, its recorded lineage is attached only if

1. both run inputs have the SHA-256 the report recorded for ``aligned_ref`` /
   ``aligned_mov``;
2. ``match_table.csv`` has the recorded SHA-256;
3. the match table has exactly one row per aligned row, output row IDs
   0…N−1 in order, and unique source row IDs on each side.

The current state of the *original* files is reported separately as
``verified`` / ``unavailable`` / ``mismatch`` and never affects attachment.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Mapping

from cryorole.provenance.source_identity import stream_sha256

REPORT_NAME = "align_report.json"


def discover_alignment_provenance(
    ref: str | Path,
    mov: str | Path,
    *,
    ref_sha256: str | None = None,
    mov_sha256: str | None = None,
    ref_row_count: int | None = None,
    mov_row_count: int | None = None,
) -> dict[str, Any] | None:
    """Return an ``alignment_provenance`` record, or ``None`` when no report is next to the inputs."""

    ref_path, mov_path = Path(ref), Path(mov)
    candidates = []
    for directory in dict.fromkeys((ref_path.resolve().parent, mov_path.resolve().parent)):
        report_path = directory / REPORT_NAME
        if report_path.is_file():
            candidates.append(report_path)
    if not candidates:
        return None
    report_path = candidates[0]
    record: dict[str, Any] = {"align_report": str(report_path), "attached": False}
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        return {**record, "reason": f"align_report.json is not readable JSON ({exc})"}
    if not isinstance(report, Mapping) or report.get("artifact_type") != "align_report":
        return {**record, "reason": "align_report.json is not a cryoROLE align report"}
    hashes = report.get("sha256") or {}
    required = ("aligned_ref", "aligned_mov", "match_table")
    if not all(hashes.get(key) for key in required):
        return {**record, "reason": "align report has no recorded hashes (created by an older cryoROLE); lineage not verifiable"}

    observed_ref = ref_sha256 or stream_sha256(ref_path)
    observed_mov = mov_sha256 or stream_sha256(mov_path)
    if observed_ref != hashes["aligned_ref"] or observed_mov != hashes["aligned_mov"]:
        return {**record, "reason": "run inputs are not the aligned files recorded in align_report.json (SHA-256 differs)"}
    table_path = Path((report.get("output_paths") or {}).get("match_table") or report_path.with_name("match_table.csv"))
    if not table_path.is_file():
        table_path = report_path.with_name("match_table.csv")
    if not table_path.is_file():
        return {**record, "reason": "match_table.csv not found"}
    if stream_sha256(table_path) != hashes["match_table"]:
        return {**record, "reason": "match_table.csv SHA-256 differs from align_report.json"}
    problem = _check_match_table(table_path, ref_row_count=ref_row_count, mov_row_count=mov_row_count)
    if problem:
        return {**record, "reason": problem}

    lineage = {
        "strategy": report.get("strategy"),
        "key_policy": report.get("key_policy"),
        "original_ref_path": report.get("ref_resolved_path") or report.get("ref_path"),
        "original_mov_path": report.get("mov_resolved_path") or report.get("mov_path"),
        "original_ref_sha256": hashes.get("ref_source"),
        "original_mov_sha256": hashes.get("mov_source"),
        "match_table": str(table_path),
        "match_table_sha256": hashes["match_table"],
        "aligned_row_count": report.get("aligned_row_count", report.get("matched_count")),
        "chain": report.get("chain"),
        "coordinate_match": report.get("coordinate_match"),
        "coverage": report.get("coverage"),
        "suspected_duplicate_groups": report.get("suspected_duplicate_groups"),
    }
    original_files = {
        "ref": _original_status(lineage["original_ref_path"], lineage["original_ref_sha256"]),
        "mov": _original_status(lineage["original_mov_path"], lineage["original_mov_sha256"]),
    }
    return {**record, "attached": True, "reason": None, "lineage": lineage, "original_files": original_files}


def _check_match_table(path: Path, *, ref_row_count: int | None, mov_row_count: int | None) -> str | None:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    n = len(rows)
    for count, label in ((ref_row_count, "ref"), (mov_row_count, "mov")):
        if count is not None and count != n:
            return f"match_table.csv has {n} rows but the aligned {label} file has {count}"
    try:
        ref_out = [int(r["ref_output_row_id"]) for r in rows]
        mov_out = [int(r["mov_output_row_id"]) for r in rows]
        ref_src = [int(r["ref_source_row_id"]) for r in rows]
        mov_src = [int(r["mov_source_row_id"]) for r in rows]
    except (KeyError, ValueError):
        return "match_table.csv lacks row-id columns"
    expected = list(range(n))
    if ref_out != expected or mov_out != expected:
        return "match_table.csv output row IDs are not 0..N-1 in order"
    if len(set(ref_src)) != n or len(set(mov_src)) != n:
        return "match_table.csv source row IDs are not unique"
    return None


def _original_status(path: str | None, expected_sha256: str | None) -> dict[str, Any]:
    if not path:
        return {"path": None, "status": "unavailable"}
    candidate = Path(path)
    if not candidate.is_file():
        return {"path": path, "status": "unavailable"}
    if expected_sha256 and stream_sha256(candidate) == expected_sha256:
        return {"path": path, "status": "verified"}
    return {"path": path, "status": "mismatch"}
