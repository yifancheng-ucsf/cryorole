"""Match particles by resolved identity keys."""

from __future__ import annotations

import pandas as pd

from cryorole.match.identity import IdentityResolutionError, ResolvedIdentity
from cryorole.match.overlap_report import overlap_status
from cryorole.models.match_table import MatchReport, MatchTable
from cryorole.models.policies import MatchPolicy


def _ensure_unique(identity: ResolvedIdentity, label: str) -> None:
    if identity.keys.duplicated(keep=False).any():
        raise IdentityResolutionError(f"{label} identity contains duplicate keys")


def match_resolved_identities(
    domain_a: ResolvedIdentity,
    domain_b: ResolvedIdentity,
    *,
    policy: MatchPolicy | None = None,
) -> tuple[MatchTable, MatchReport]:
    """Join two domains by explicit resolved identity keys."""

    policy = policy or MatchPolicy()
    if policy.join_type != "inner":
        raise ValueError("Phase 1 supports only inner identity matching")
    if policy.duplicate_handling == "fail":
        _ensure_unique(domain_a, "domain_a")
        _ensure_unique(domain_b, "domain_b")

    a = pd.DataFrame({"particle_key": domain_a.keys, "domain_a_row": domain_a.keys.index})
    b = pd.DataFrame({"particle_key": domain_b.keys, "domain_b_row": domain_b.keys.index})
    matched = a.merge(b, on="particle_key", how="inner")
    matched["match_status"] = "matched"

    matched_count = len(matched)
    unmatched_a = len(a) - matched_count
    unmatched_b = len(b) - matched_count
    denominator = min(len(a), len(b)) if min(len(a), len(b)) else 1
    ratio = matched_count / denominator
    status = overlap_status(
        ratio,
        threshold=policy.overlap_threshold,
        low_overlap_behavior=policy.low_overlap_behavior,
    )
    matched_rows_reordered = bool(
        matched_count
        and not (matched["domain_a_row"].to_numpy() == matched["domain_b_row"].to_numpy()).all()
    )
    identity_columns = domain_a.report.identity_columns or domain_b.report.identity_columns
    match_key = "+".join(identity_columns) if identity_columns else None
    warnings: list[str] = []
    if matched_rows_reordered:
        warnings.append(
            "inputs_matched_by_key_and_reordered_before_ro_computation"
        )
    if unmatched_a or unmatched_b:
        warnings.append(
            "inputs_matched_by_key_with_unmatched_particles_dropped"
        )

    report = MatchReport(
        matched_count=matched_count,
        unmatched_a=unmatched_a,
        unmatched_b=unmatched_b,
        overlap_ratio=ratio,
        status=status,
        match_key=match_key,
        ref_row_count=len(a),
        mov_row_count=len(b),
        matched_row_count=matched_count,
        dropped_ref_only_count=unmatched_a,
        dropped_mov_only_count=unmatched_b,
        matched_rows_reordered=matched_rows_reordered,
        warnings=tuple(warnings),
    )
    if status == "failed_low_overlap":
        raise IdentityResolutionError("Particle overlap is below MatchPolicy threshold")

    return MatchTable(matched), report
