"""Input-sanity heuristics for a two-body RO run.

cryoROLE assumes the user has refined two near-rigid bodies and does not
decide which pose field "is" a domain. It can still notice patterns that
usually mean an input mistake, such as the same refinement selected twice.

The RO-angle summary is always reported as information. Warnings are raised
only for the tiered patterns below. All thresholds are explicit heuristics,
recorded with every result, and never change the analysis:

* same file for ``--ref`` and ``--mov`` (same resolved path or SHA-256):
  blocked in ``preflight``, strong warning in ``run``;
* >= ``identical_min_fraction`` of particles with RO angle below
  ``identical_angle_rad``: strong warning (identical poses);
* median RO angle below ``near_identical_median_deg`` *and* 99th percentile
  below ``near_identical_p99_deg``: warning (nearly identical orientations).

Real two-body data are also close to identity (the J75/J80 example has a
median of 8.8 degrees), so the "nearly identical" rule is deliberately narrow
and its wording deliberately non-accusatory: a genuinely rigid complex can
produce a narrow distribution without any user mistake.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

INPUT_SANITY_SCHEMA_VERSION = "1"
_PERCENTILES = (1, 10, 50, 90, 99)


@dataclass(frozen=True)
class InputSanityPolicy:
    """Explicit heuristic thresholds for the RO-angle sanity checks."""

    identical_angle_rad: float = 1e-6
    identical_min_fraction: float = 0.99
    near_identical_median_deg: float = 1.0
    near_identical_p99_deg: float = 2.0

    def __post_init__(self) -> None:
        if not self.identical_angle_rad > 0:
            raise ValueError("identical_angle_rad must be positive")
        if not 0 < self.identical_min_fraction <= 1:
            raise ValueError("identical_min_fraction must be in (0, 1]")
        if not (self.near_identical_median_deg > 0 and self.near_identical_p99_deg > 0):
            raise ValueError("near-identical thresholds must be positive")


def summarize_ro_angles(angle_rad: np.ndarray) -> dict[str, Any]:
    """Return the always-reported RO-angle summary (degrees)."""

    angles = np.asarray(angle_rad, dtype=float).reshape(-1)
    n = int(angles.size)
    if n == 0:
        return {"n": 0, "percentiles_deg": {}, "median_deg": None, "p90_deg": None, "p99_deg": None,
                "message": "RO angle: no matched particles."}
    degrees = np.degrees(angles)
    values = np.percentile(degrees, _PERCENTILES)
    percentiles = {f"p{p}": float(v) for p, v in zip(_PERCENTILES, values)}
    summary = {
        "n": n,
        "percentiles_deg": percentiles,
        "median_deg": percentiles["p50"],
        "p90_deg": percentiles["p90"],
        "p99_deg": percentiles["p99"],
        "fraction_below_1_deg": float(np.mean(degrees < 1.0)),
        "fraction_below_5_deg": float(np.mean(degrees < 5.0)),
    }
    summary["message"] = (
        f"RO angle: median {_deg(summary['median_deg'])}, 90% {_deg(summary['p90_deg'])}, "
        f"99% {_deg(summary['p99_deg'])} (n = {n:,})."
    )
    return summary


def same_input_file(
    ref_identity: Mapping[str, Any] | None,
    mov_identity: Mapping[str, Any] | None,
    *,
    ref_path: str | Path | None = None,
    mov_path: str | Path | None = None,
) -> dict[str, Any] | None:
    """Return a finding when ref and mov are the same file, else ``None``.

    Uses the resolved path first, then the SHA-256 recorded in the source
    identities (so a copy of the same file under another name is caught).
    """

    ref_identity = ref_identity or {}
    mov_identity = mov_identity or {}
    ref_resolved = ref_identity.get("resolved_path") or _resolve(ref_path)
    mov_resolved = mov_identity.get("resolved_path") or _resolve(mov_path)
    reason = None
    if ref_resolved and mov_resolved and str(ref_resolved) == str(mov_resolved):
        reason = "same_resolved_path"
    else:
        ref_hash = ref_identity.get("sha256")
        mov_hash = mov_identity.get("sha256")
        if ref_hash and mov_hash and ref_hash == mov_hash:
            reason = "same_sha256"
    if reason is None:
        return None
    return {
        "code": "SAME_INPUT_FILE",
        "reason": reason,
        "ref_resolved_path": str(ref_resolved) if ref_resolved else None,
        "mov_resolved_path": str(mov_resolved) if mov_resolved else None,
        "message": (
            "The reference and moving inputs are the same file"
            + (" (identical content)" if reason == "same_sha256" else "")
            + ". Select the two different refinements."
        ),
    }


def assess_ro_angles(
    angle_rad: np.ndarray | None,
    *,
    policy: InputSanityPolicy | None = None,
    same_file: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Assess RO angles (and an optional same-file finding) into one report."""

    policy = policy or InputSanityPolicy()
    if angle_rad is None:
        angles = np.empty(0, dtype=float)
        summary = {"n": None, "message": "RO angle summary unavailable for this run backend."}
    else:
        angles = np.asarray(angle_rad, dtype=float).reshape(-1)
        summary = summarize_ro_angles(angles)
    findings: list[dict[str, Any]] = []
    if same_file:
        findings.append({"level": "strong_warning", **dict(same_file)})
    if angles.size:
        identical_fraction = float(np.mean(angles < policy.identical_angle_rad))
        summary["fraction_identical"] = identical_fraction
        if identical_fraction >= policy.identical_min_fraction:
            findings.append(
                {
                    "level": "strong_warning",
                    "code": "IDENTICAL_POSES",
                    "fraction": identical_fraction,
                    "message": (
                        f"Reference and moving poses are identical for {_pct(identical_fraction)} of "
                        "particles. The two inputs are probably the same refinement (e.g. the same "
                        "job selected twice)."
                    ),
                }
            )
        elif (
            summary["median_deg"] < policy.near_identical_median_deg
            and summary["p99_deg"] < policy.near_identical_p99_deg
        ):
            findings.append(
                {
                    "level": "warning",
                    "code": "NEARLY_IDENTICAL_ORIENTATIONS",
                    "message": (
                        "The two refinements give nearly identical orientations "
                        f"(median {_deg(summary['median_deg'])}, 99% {_deg(summary['p99_deg'])}). "
                        "Either there is no measurable inter-domain motion, or both refinements "
                        "aligned the same region. Check the masks."
                    ),
                }
            )
    level = "none"
    if any(f["level"] == "strong_warning" for f in findings):
        level = "strong_warning"
    elif findings:
        level = "warning"
    return {
        "artifact_type": "cryorole_input_sanity",
        "schema_version": INPUT_SANITY_SCHEMA_VERSION,
        "level": level,
        "ro_angle_summary": summary,
        "findings": findings,
        "policy": {"heuristic": True, **asdict(policy)},
    }


def sanity_warning_lines(report: Mapping[str, Any] | None) -> tuple[str, ...]:
    """Return one ``[CODE] message`` line per finding."""

    if not report:
        return ()
    return tuple(
        f"[{finding['code']}] {finding['message']}" for finding in report.get("findings", ())
    )


def _resolve(path: str | Path | None) -> str | None:
    if path is None:
        return None
    try:
        return str(Path(path).expanduser().resolve())
    except OSError:
        return None


def _deg(value: float | None) -> str:
    if value is None:
        return "n/a"
    return f"{value:.1f}°" if value >= 0.1 else f"{value:.2g}°"


def _pct(fraction: float) -> str:
    return f"{fraction:.1%}"
