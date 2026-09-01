"""Particle identity resolution and matching."""

from cryorole.match.identity import (
    IdentityResolutionError,
    ResolvedIdentity,
    resolve_identity,
)
from cryorole.match.matcher import match_resolved_identities
from cryorole.match.preflight import DEFAULT_MIN_MATCH_OVERLAP, public_match_policy

__all__ = [
    "IdentityResolutionError",
    "ResolvedIdentity",
    "match_resolved_identities",
    "DEFAULT_MIN_MATCH_OVERLAP",
    "public_match_policy",
    "resolve_identity",
]
