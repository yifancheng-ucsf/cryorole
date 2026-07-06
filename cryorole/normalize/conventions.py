"""Centralized source convention normalization.

Phase 2 intentionally supports RELION conventions only. New source convention
support must be added here, not in readers, matching, core computation, or
frontends.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial.transform import Rotation

from cryorole.models.policies import ConventionPolicy


IMPLEMENTED_RELION_CONVERSION_RULE = (
    "active_matrix = scipy Rotation.from_euler('ZYZ', angles).as_matrix().T"
)


@dataclass(frozen=True)
class ConventionResolution:
    """Resolved convention details for audit and future manifest writing."""

    policy: ConventionPolicy
    source_sequence_used: str
    conversion_rule_used: str


def relion_passive_euler_to_active_matrix(
    rot: float,
    tilt: float,
    psi: float,
) -> np.ndarray:
    """Convert RELION passive intrinsic ZYZ Euler angles to an active matrix."""

    passive_matrix = Rotation.from_euler("ZYZ", [rot, tilt, psi], degrees=True).as_matrix()
    return passive_matrix.T


class ConventionResolver:
    """Resolve source pose conventions into internal active matrices.

    Current scope is RELION-only. This class rejects any non-RELION convention
    so that future source support is added explicitly and testably.
    """

    def __init__(self, policy: ConventionPolicy) -> None:
        self.policy = policy
        self._validate_policy()

    def _validate_policy(self) -> None:
        if self.policy.source_software != "relion":
            raise ValueError(f"Unsupported source software: {self.policy.source_software}")
        if self.policy.source_euler_sequence != "ZYZ":
            raise ValueError("RELION Euler parsing must use uppercase 'ZYZ'")
        if self.policy.source_semantics != "passive":
            raise ValueError("RELION source semantics must be passive")
        if self.policy.internal_semantics != "active":
            raise ValueError("Internal rotation semantics must be active")
        if self.policy.conversion_rule != IMPLEMENTED_RELION_CONVERSION_RULE:
            raise ValueError(
                "ConventionPolicy conversion_rule does not match the implemented "
                "RELION passive-to-active bridge"
            )

    def resolve(self) -> ConventionResolution:
        """Return the auditable convention resolution."""

        return ConventionResolution(
            policy=self.policy,
            source_sequence_used="ZYZ",
            conversion_rule_used=self.policy.conversion_rule,
        )

    def euler_to_active_matrix(self, rot: float, tilt: float, psi: float) -> np.ndarray:
        """Convert source Euler angles to an internal active matrix."""

        return relion_passive_euler_to_active_matrix(rot, tilt, psi)
