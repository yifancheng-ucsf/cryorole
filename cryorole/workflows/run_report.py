"""Concise human-readable navigation report for a completed run."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from cryorole.models.density_report import DensityReport


@dataclass(frozen=True)
class RunReportRequest:
    output_path: str | Path
    summary: Mapping[str, Any]
    density_report: DensityReport | None
    overwrite: bool = False


def write_run_report(request: RunReportRequest) -> Path:
    """Write a compact Markdown guide derived from typed run results."""

    path = Path(request.output_path)
    if path.exists() and not request.overwrite:
        raise FileExistsError(f"Output path already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_render_run_report(request), encoding="utf-8")
    return path


def _render_run_report(request: RunReportRequest) -> str:
    summary = request.summary
    density = request.density_report
    lines = [
        "# cryoROLE Run Report",
        "",
        f"Run ID: `{summary.get('run_id')}`",
        f"Matched particles: **{summary.get('matched_count')}**",
        "",
        "## Inputs and matching",
        "",
        f"- Reference: `{summary.get('input_paths', {}).get('ref')}` ({summary.get('source_types', {}).get('ref')})",
        f"- Moving: `{summary.get('input_paths', {}).get('mov')}` ({summary.get('source_types', {}).get('mov')})",
        f"- Policy: `{summary.get('match_key')}`; reordered={summary.get('matched_rows_reordered')}",
        f"- Dropped: ref={summary.get('dropped_ref_only_count')}, mov={summary.get('dropped_mov_only_count')}",
        *_alignment_lines(summary.get("alignment_provenance")),
        "",
        "## Analysis",
        "",
        "- Relative orientation: `RO = R_ref^-1 R_mov`",
        f"- Euler convention: `{summary.get('euler_convention')}`",
        f"- SLD: `{summary.get('resolved_sld_metric')}`, k={summary.get('k_neighbors')}",
        "",
        *_input_sanity_lines(summary.get("input_sanity")),
        "## Core files",
        "",
        "- `data/raw_landscape.npz` — machine-readable landscape",
        "- `data/raw_landscape.csv` — user-facing table",
        "- `data/match_table.csv` — source-row provenance",
        "- `reports/density_report.json` — SLD diagnostics",
        "- `run_summary.json` and `run_manifest.json` — policy and provenance",
        "",
        "## Quick-look",
        "",
    ]
    if summary.get("raw_visualization_performed"):
        lines.extend(
            [
                "- `visualizations/quicklook/all_euler_3view_projection.png` — all particles in Euler coordinates",
                "- `visualizations/quicklook/all_rotvec_3view_projection.png` — all particles in rotation-vector coordinates",
                "- `visualizations/quicklook/sld_ge_1_euler_3view_projection.png` — particles with `sld_raw >= 1`",
                "- `visualizations/quicklook/sld_ge_1_rotvec_3view_projection.png` — the same density range in rotation-vector coordinates",
                "- `visualizations/quicklook/top_40pct_euler_3view_projection.png` — highest-density 40% with cutoff ties retained",
                "- `visualizations/quicklook/top_40pct_rotvec_3view_projection.png` — the same high-density range in rotation-vector coordinates",
                "- `visualizations/quicklook/sld_log_distribution.png` — full-landscape `log10(sld_raw)` distribution",
                "",
                "Colors saturate above SLD 100; stored SLD values are not clipped.",
            ]
        )
    else:
        lines.append("Quick-look was disabled with `--no-visualize`.")
    lines.extend(["", "## Diagnostics", ""])
    warnings = list(summary.get("match_warnings") or ())
    if density is not None:
        lines.append(
            f"SLD: P99={density.p99_sld_raw:g}, max={density.max_sld_raw:g}, "
            f">100={density.n_high_sld_points}/{density.n_points}, "
            f"distance-floored={density.n_floored_points}/{density.n_points}."
        )
        if density.n_inf_sld_unfloored:
            lines.extend(
                [
                    "",
                    f"{density.n_inf_sld_unfloored}/{density.n_points} particles have all of their "
                    f"k={density.effective_k_neighbors} nearest neighbours at exactly the same relative "
                    "orientation, so their unfloored SLD (`sld_unfloored`) is +inf. "
                    + (
                        "`sld_raw` stays finite for them because of the distance floor. "
                        if not density.n_inf_sld_raw
                        else f"{density.n_inf_sld_raw} of them also have +inf `sld_raw`. "
                    )
                    + "This usually means duplicated particles, symmetry expansion, or a rigid "
                    "subpopulation larger than k.",
                ]
            )
        diagnostic = density.ro_coordinate_diagnostics
        if diagnostic and diagnostic["severity"] == "info":
            lines.extend(["", str(diagnostic["message"])])
        warnings.extend(density.warnings)
    if warnings:
        lines.append("")
        lines.extend(f"- {warning}" for warning in warnings)
    else:
        lines.append("No matching or SLD warnings were reported.")
    lines.extend(
        [
            "",
            "## Next",
            "",
            "```bash",
            "cryorole status --run-dir .",
            "cryorole visualize --run-dir .",
            "cryorole canonicalize --run-dir .",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def _input_sanity_lines(sanity: Mapping[str, Any] | None) -> list[str]:
    if not sanity:
        return []
    lines = ["## Input sanity", ""]
    summary = sanity.get("ro_angle_summary") or {}
    if summary.get("message"):
        lines.append(str(summary["message"]))
    findings = list(sanity.get("findings") or ())
    if findings:
        lines.append("")
        for finding in findings:
            label = "STRONG WARNING" if finding.get("level") == "strong_warning" else "WARNING"
            lines.append(f"- **{label}** [{finding.get('code')}] {finding.get('message')}")
    else:
        lines.extend(["", "No input-sanity warnings (thresholds are heuristics recorded in `run_summary.json`)."])
    lines.append("")
    return lines


def _alignment_lines(provenance: Mapping[str, Any] | None) -> list[str]:
    if not provenance:
        return []
    if not provenance.get("attached"):
        return [f"- Alignment lineage: not attached ({provenance.get('reason')})"]
    lineage = provenance.get("lineage") or {}
    states = provenance.get("original_files") or {}
    status = ", ".join(f"{k} {v.get('status')}" for k, v in states.items())
    return [f"- Alignment lineage: `{lineage.get('strategy')}` from `{provenance.get('align_report')}` (originals: {status})"]
