"""Single-source public matching safety policy."""

from __future__ import annotations

from cryorole.models.policies import DEFAULT_MIN_MATCH_OVERLAP, MatchPolicy


def public_match_policy(*, allow_low_overlap: bool = False) -> MatchPolicy:
    """Return the public run matching policy from one centralized default."""

    return MatchPolicy(
        overlap_threshold=DEFAULT_MIN_MATCH_OVERLAP,
        low_overlap_behavior="warn" if allow_low_overlap else "fail",
        low_overlap_allowed=allow_low_overlap,
    )
