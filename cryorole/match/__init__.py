"""Particle identity resolution and matching."""

from cryorole.match.identity import (
    IdentityResolutionError,
    ResolvedIdentity,
    resolve_identity,
)
from cryorole.match.matcher import match_resolved_identities

__all__ = [
    "IdentityResolutionError",
    "ResolvedIdentity",
    "match_resolved_identities",
    "resolve_identity",
]

