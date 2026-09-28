"""Array-native production core for ``cryorole run``."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from cryorole.core.euler_conventions import resolve_euler_convention
from cryorole.io.readers import read_relion_star_particle_columns
from cryorole.match import ResolvedIdentity, match_resolved_identities, resolve_identity
from cryorole.models.match_table import MatchReport, MatchTable
from cryorole.models.policies import ConventionPolicy, IdentityPolicy, MatchPolicy
from cryorole.models.pose_arrays import MatchedPoseArrays, ROArrays
from cryorole.normalize.conventions import ConventionResolver
from cryorole.normalize.field_mapper import CRYOSPARC_POSE_FIELD, RELION_EULER_FIELDS


@dataclass(frozen=True)
class ArrayNativePreflightResult:
    """Matched and normalized production inputs plus audit records."""

    matched_poses: MatchedPoseArrays
    identity_ref: ResolvedIdentity
    identity_mov: ResolvedIdentity
    match_table: MatchTable
    match_report: MatchReport
    ref_row_count: int
    mov_row_count: int
    source_ref: str
    source_mov: str


@dataclass(frozen=True)
class _RawSource:
    source_type: str
    row_count: int
    identity_table: pd.DataFrame
    pose_source: Any


def inspect_source_identity(
    path: str | Path,
    *,
    identity_policy: IdentityPolicy,
) -> tuple[int, str, ResolvedIdentity]:
    """Read minimal fields and return an identity report even for duplicates."""

    source = _read_minimal_source(path, identity_policy)
    report_policy = replace(identity_policy, failure_mode="report")
    identity = resolve_identity(
        source.identity_table,
        source_type=source.source_type,
        policy=report_policy,
    )
    return source.row_count, source.source_type, identity


def preflight_and_normalize_matched_arrays(
    ref_path: str | Path,
    mov_path: str | Path,
    *,
    ref_domain: str,
    mov_domain: str,
    identity_policy_ref: IdentityPolicy,
    identity_policy_mov: IdentityPolicy,
    match_policy: MatchPolicy,
    convention_policy_ref: ConventionPolicy | None = None,
    convention_policy_mov: ConventionPolicy | None = None,
) -> ArrayNativePreflightResult:
    """Read identity/pose fields, match first, then normalize matched rows only."""

    ref = _read_minimal_source(ref_path, identity_policy_ref)
    mov = _read_minimal_source(mov_path, identity_policy_mov)
    _validate_row_aligned_counts(
        ref.row_count,
        mov.row_count,
        identity_policy_ref,
        identity_policy_mov,
    )
    identity_ref = resolve_identity(
        ref.identity_table,
        source_type=ref.source_type,
        policy=identity_policy_ref,
    )
    identity_mov = resolve_identity(
        mov.identity_table,
        source_type=mov.source_type,
        policy=identity_policy_mov,
    )
    match_table, match_report = match_resolved_identities(
        identity_ref,
        identity_mov,
        policy=match_policy,
    )
    ref_rows = match_table.data["domain_a_row"].to_numpy(dtype=np.int64)
    mov_rows = match_table.data["domain_b_row"].to_numpy(dtype=np.int64)
    keys = match_table.data["particle_key"].astype(str).to_numpy()
    ref_matrices = _normalize_rows(
        ref,
        ref_rows,
        convention_policy=convention_policy_ref,
    )
    mov_matrices = _normalize_rows(
        mov,
        mov_rows,
        convention_policy=convention_policy_mov,
    )
    matched = MatchedPoseArrays(
        particle_key=keys,
        ref_rotation_matrix_active=ref_matrices,
        mov_rotation_matrix_active=mov_matrices,
        ref_source_row_id=ref_rows,
        mov_source_row_id=mov_rows,
    )
    return ArrayNativePreflightResult(
        matched_poses=matched,
        identity_ref=identity_ref,
        identity_mov=identity_mov,
        match_table=match_table,
        match_report=match_report,
        ref_row_count=ref.row_count,
        mov_row_count=mov.row_count,
        source_ref=ref.source_type,
        source_mov=mov.source_type,
    )


def compute_relative_orientation_arrays(
    matched: MatchedPoseArrays,
    *,
    euler_sequence: str | None = None,
) -> ROArrays:
    """Compute ``R_ref.T @ R_mov`` and all derived representations in batches."""

    matrices = np.matmul(
        np.swapaxes(matched.ref_rotation_matrix_active, 1, 2),
        matched.mov_rotation_matrix_active,
    )
    rotations = Rotation.from_matrix(matrices)
    resolved = resolve_euler_convention(scipy_euler_sequence=euler_sequence)
    return ROArrays(
        particle_key=matched.particle_key,
        rotation_matrix=matrices,
        quaternion_xyzw=rotations.as_quat(),
        rotation_vector=rotations.as_rotvec(),
        euler_zyx=rotations.as_euler(resolved.scipy_euler_sequence, degrees=True),
        angle_rad=rotations.magnitude(),
        ref_source_row_id=matched.ref_source_row_id,
        mov_source_row_id=matched.mov_source_row_id,
    )


def _read_minimal_source(path: str | Path, identity_policy: IdentityPolicy) -> _RawSource:
    source_path = Path(path)
    suffix = source_path.suffix.lower()
    if suffix == ".star":
        identity_columns: tuple[str, ...]
        if identity_policy.identity_mode == "relion_image_name":
            identity_columns = ("_rlnTomoParticleName", "_rlnImageName", "rlnImageName")
        elif identity_policy.identity_mode == "relion_user_columns":
            identity_columns = tuple(identity_policy.identity_columns)
        else:
            identity_columns = ()
        raw = read_relion_star_particle_columns(
            source_path,
            columns=tuple(dict.fromkeys((*RELION_EULER_FIELDS, *identity_columns))),
        )
        _validate_pose_source_schema("relion", raw.particles)
        return _RawSource("relion", len(raw.particles), raw.particles, raw.particles)
    if suffix == ".cs":
        array = np.load(source_path, allow_pickle=False, mmap_mode="r")
        if not array.dtype.names:
            raise ValueError("CryoSPARC .cs input must be a structured array")
        required = {CRYOSPARC_POSE_FIELD}
        if identity_policy.identity_mode != "row_aligned":
            required.add("uid")
        missing = sorted(required - set(array.dtype.names))
        if missing:
            raise ValueError(f"CryoSPARC input missing required fields: {missing}")
        _validate_pose_source_schema("cryosparc", array)
        identity_columns: dict[str, Any] = {}
        if "uid" in array.dtype.names:
            identity_columns["uid"] = np.asarray(array["uid"])
        identity_table = pd.DataFrame(identity_columns, index=np.arange(len(array)))
        return _RawSource("cryosparc", len(array), identity_table, array)
    raise ValueError(f"Unsupported input suffix for run: {suffix}")


def _validate_pose_source_schema(source_type: str, pose_source: Any) -> None:
    """Validate all native pose values without performing convention conversion."""

    if source_type == "relion":
        missing = [field for field in RELION_EULER_FIELDS if field not in pose_source.columns]
        if missing:
            raise ValueError(f"RELION particles missing required pose fields: {missing}")
        values = pose_source.loc[:, list(RELION_EULER_FIELDS)].to_numpy(dtype=float)
        if values.ndim != 2 or values.shape[1] != 3:
            raise ValueError("RELION Rot/Tilt/Psi pose values must have shape (n, 3)")
    else:
        values = np.asarray(pose_source[CRYOSPARC_POSE_FIELD])
        if values.ndim != 2 or values.shape[1] != 3:
            raise ValueError(
                f"CryoSPARC {CRYOSPARC_POSE_FIELD} must have shape (n, 3), got {values.shape}"
            )
        values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"{source_type} pose values contain NaN or Inf")


def _normalize_rows(
    source: _RawSource,
    row_ids: np.ndarray,
    *,
    convention_policy: ConventionPolicy | None,
) -> np.ndarray:
    if source.source_type == "relion":
        table = source.pose_source
        missing = [field for field in RELION_EULER_FIELDS if field not in table.columns]
        if missing:
            raise ValueError(f"RELION particles missing required pose fields: {missing}")
        angles = table.loc[row_ids, list(RELION_EULER_FIELDS)].to_numpy(dtype=float)
        resolver = ConventionResolver(convention_policy or ConventionPolicy.relion_default())
        return resolver.euler_to_active_matrices(angles)
    poses = np.asarray(source.pose_source[CRYOSPARC_POSE_FIELD][row_ids], dtype=float)
    if poses.shape != (len(row_ids), 3):
        raise ValueError(
            f"CryoSPARC {CRYOSPARC_POSE_FIELD} must have shape (n, 3), got {poses.shape}"
        )
    if not np.isfinite(poses).all():
        raise ValueError(f"CryoSPARC {CRYOSPARC_POSE_FIELD} contains non-finite values")
    resolver = ConventionResolver(convention_policy or ConventionPolicy.cryosparc_default())
    return resolver.rotvec_to_active_matrices(poses)


def _validate_row_aligned_counts(
    ref_count: int,
    mov_count: int,
    ref_policy: IdentityPolicy,
    mov_policy: IdentityPolicy,
) -> None:
    ref_aligned = ref_policy.identity_mode == "row_aligned"
    mov_aligned = mov_policy.identity_mode == "row_aligned"
    if ref_aligned != mov_aligned:
        raise ValueError("row_aligned identity mode must be requested for both domains")
    if ref_aligned and ref_count != mov_count:
        raise ValueError(
            "row_aligned identity mode requires equal reference and moving row counts; "
            f"got {ref_count} and {mov_count}"
        )
