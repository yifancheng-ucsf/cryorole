"""Field mappings from source metadata to normalized pose fields."""

from __future__ import annotations

RELION_EULER_FIELDS = ("_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi")
RELION_SHIFT_FIELDS = ("_rlnOriginXAngst", "_rlnOriginYAngst")

CRYOSPARC_UID_FIELD = "uid"
CRYOSPARC_POSE_FIELD = "alignments3D/pose"
CRYOSPARC_SHIFT_FIELD = "alignments3D/shift"


def relion_required_pose_fields() -> tuple[str, ...]:
    """Return RELION fields required for Phase 1 pose normalization."""

    return RELION_EULER_FIELDS
