"""Centralized source convention normalization.

Supported sources are RELION (passive intrinsic ZYZ Euler angles) and
CryoSPARC (``alignments3D/pose`` axis-angle vectors). New source convention
support must be added here, not in readers, matching, core computation, or
frontends.

CryoSPARC bridge: the CryoSPARC pose is the rotation applied during image
back-projection, the inverse (transpose) of RELION's reference-to-image
rotation. The bridge follows pyem ``csparc2star.py``
(``rot2euler(expmap(pose))``; pyem's ``expmap(v)`` equals SciPy
``Rotation.from_rotvec(v).as_matrix().T``), so a ``.cs`` input and its
pyem-converted STAR give the same internal active matrix. Rotation vectors of
any norm (including > pi) are valid; ``from_rotvec`` handles them exactly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial.transform import Rotation

from cryorole.models.policies import ConventionPolicy


IMPLEMENTED_RELION_CONVERSION_RULE = (
    "active_matrix = scipy Rotation.from_euler('ZYZ', angles).as_matrix().T"
)
IMPLEMENTED_CRYOSPARC_CONVERSION_RULE = (
    "active_matrix = scipy Rotation.from_rotvec(pose).as_matrix().T"
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


def relion_passive_euler_to_active_matrices(angles: np.ndarray) -> np.ndarray:
    """Vectorized RELION passive intrinsic-ZYZ to internal active matrices."""

    values = np.asarray(angles, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("RELION Euler angles must have shape (n, 3)")
    if not np.isfinite(values).all():
        raise ValueError("RELION Euler angles contain non-finite values")
    passive = Rotation.from_euler("ZYZ", values, degrees=True).as_matrix()
    return np.swapaxes(passive, 1, 2)


def cryosparc_rotvec_to_active_matrices(poses: np.ndarray) -> np.ndarray:
    """Vectorized CryoSPARC ``alignments3D/pose`` to internal active matrices."""

    values = np.asarray(poses, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("CryoSPARC pose rotation vectors must have shape (n, 3)")
    if not np.isfinite(values).all():
        raise ValueError("CryoSPARC pose rotation vectors contain non-finite values")
    return np.swapaxes(Rotation.from_rotvec(values).as_matrix(), 1, 2)


class ConventionResolver:
    """Resolve source pose conventions into internal active matrices.

    Supports RELION and CryoSPARC. Any other source, or a policy whose
    conversion rule does not match the implemented bridge, is rejected so that
    new source support is added explicitly and testably.
    """

    def __init__(self, policy: ConventionPolicy) -> None:
        self.policy = policy
        self._validate_policy()

    def _validate_policy(self) -> None:
        if self.policy.source_software == "cryosparc":
            self._validate_cryosparc_policy()
            return
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

    def _validate_cryosparc_policy(self) -> None:
        if self.policy.source_euler_sequence != "rotvec":
            raise ValueError("CryoSPARC poses must be parsed as rotation vectors ('rotvec')")
        if self.policy.source_semantics != "backprojection":
            raise ValueError("CryoSPARC source semantics must be 'backprojection'")
        if self.policy.internal_semantics != "active":
            raise ValueError("Internal rotation semantics must be active")
        if self.policy.conversion_rule != IMPLEMENTED_CRYOSPARC_CONVERSION_RULE:
            raise ValueError(
                "ConventionPolicy conversion_rule does not match the implemented "
                "CryoSPARC pyem-consistent bridge"
            )

    def resolve(self) -> ConventionResolution:
        """Return the auditable convention resolution."""

        return ConventionResolution(
            policy=self.policy,
            source_sequence_used=(
                "rotvec" if self.policy.source_software == "cryosparc" else "ZYZ"
            ),
            conversion_rule_used=self.policy.conversion_rule,
        )

    def _require(self, software: str) -> None:
        if self.policy.source_software != software:
            raise ValueError(
                f"Convention resolver for {self.policy.source_software} cannot convert {software} poses"
            )

    def euler_to_active_matrix(self, rot: float, tilt: float, psi: float) -> np.ndarray:
        """Convert source Euler angles to an internal active matrix."""

        self._require("relion")
        return relion_passive_euler_to_active_matrix(rot, tilt, psi)

    def euler_to_active_matrices(self, angles: np.ndarray) -> np.ndarray:
        """Convert a batch while preserving the same centralized convention bridge."""

        self._require("relion")
        return relion_passive_euler_to_active_matrices(angles)

    def rotvec_to_active_matrices(self, poses: np.ndarray) -> np.ndarray:
        """Convert CryoSPARC pose rotation vectors to internal active matrices."""

        self._require("cryosparc")
        return cryosparc_rotvec_to_active_matrices(poses)
