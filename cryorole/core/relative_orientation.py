"""Relative-orientation computation from normalized matched inputs."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cryorole.core.representations import derive_rotation_representations
from cryorole.models.match_table import MatchTable
from cryorole.models.pose_table import PoseTable
from cryorole.models.ro_result import ROResult


def _require_normalized_matched_inputs(
    reference: PoseTable,
    moving: PoseTable,
    matches: MatchTable,
) -> None:
    if not isinstance(reference, PoseTable):
        raise TypeError("reference must be a normalized PoseTable")
    if not isinstance(moving, PoseTable):
        raise TypeError("moving must be a normalized PoseTable")
    if not isinstance(matches, MatchTable):
        raise TypeError("matches must be a MatchTable produced by explicit matching")
    if matches.data.empty:
        raise ValueError("matches must contain at least one matched particle")
    if not (matches.data["match_status"] == "matched").all():
        raise ValueError("RO computation accepts matched rows only")


def compute_relative_orientations(
    reference: PoseTable,
    moving: PoseTable,
    matches: MatchTable,
    *,
    euler_sequence: str | None = None,
    degrees: bool = True,
) -> ROResult:
    """Compute RO as ``R_ref^-1 R_mov`` from normalized active matrices."""

    _require_normalized_matched_inputs(reference, moving, matches)
    reference_by_source_row = reference.data.set_index("source_row_id", drop=False)
    moving_by_source_row = moving.data.set_index("source_row_id", drop=False)
    rows = []
    for _, match in matches.data.iterrows():
        ref_row = reference_by_source_row.loc[int(match["domain_a_row"])]
        mov_row = moving_by_source_row.loc[int(match["domain_b_row"])]
        if ref_row["particle_key"] != match["particle_key"]:
            raise ValueError("reference PoseTable particle_key does not match MatchTable")
        if mov_row["particle_key"] != match["particle_key"]:
            raise ValueError("moving PoseTable particle_key does not match MatchTable")

        r_ref = np.asarray(ref_row["rotation_matrix_active"], dtype=float)
        r_mov = np.asarray(mov_row["rotation_matrix_active"], dtype=float)
        r_ro = np.linalg.inv(r_ref) @ r_mov
        derived = derive_rotation_representations(
            r_ro,
            euler_sequence=euler_sequence,
            degrees=degrees,
        )
        rows.append(
            {
                "particle_key": match["particle_key"],
                "R_ro": r_ro,
                "reference_source_row_id": ref_row["source_row_id"],
                "moving_source_row_id": mov_row["source_row_id"],
                **derived,
            }
        )
    return ROResult(pd.DataFrame(rows))
