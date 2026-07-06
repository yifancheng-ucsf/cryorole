"""Normalize raw source tables into PoseTable objects."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation

from cryorole.models.pose_table import PoseTable
from cryorole.normalize.conventions import ConventionResolver
from cryorole.normalize.field_mapper import (
    CRYOSPARC_POSE_FIELD,
    CRYOSPARC_SHIFT_FIELD,
    CRYOSPARC_UID_FIELD,
    RELION_EULER_FIELDS,
    RELION_SHIFT_FIELDS,
    relion_required_pose_fields,
)
from cryorole.normalize.validators import require_columns


class PoseNormalizer:
    """Convert parsed raw metadata to the shared internal pose schema."""

    def normalize_relion(
        self,
        particles: pd.DataFrame,
        *,
        domain_name: str,
        convention_resolver: ConventionResolver,
    ) -> PoseTable:
        """Normalize RELION particles using the centralized convention resolver."""

        require_columns(particles, relion_required_pose_fields(), context="RELION particles")
        rows = []
        for source_row_id, row in particles.iterrows():
            source_row_id_int = int(source_row_id)
            rot, tilt, psi = (float(row[field]) for field in RELION_EULER_FIELDS)
            matrix = convention_resolver.euler_to_active_matrix(rot, tilt, psi)
            if all(field in particles.columns for field in RELION_SHIFT_FIELDS):
                shift_xy = np.array([float(row[RELION_SHIFT_FIELDS[0]]), float(row[RELION_SHIFT_FIELDS[1]])])
            else:
                shift_xy = None
            rows.append(
                {
                    "particle_key": None,
                    "domain_name": domain_name,
                    "rotation_matrix_active": matrix,
                    "shift_xy": shift_xy,
                    "source_type": "relion",
                    "source_row_id": source_row_id_int,
                    "uid": None,
                    "raw_metadata": {
                        **row.to_dict(),
                        "_cryorole_source_row_id": source_row_id_int,
                    },
                }
            )
        return PoseTable(pd.DataFrame(rows))

    def normalize_cryosparc(
        self,
        particles: pd.DataFrame,
        *,
        domain_name: str,
        pose_field: str = CRYOSPARC_POSE_FIELD,
        shift_field: str = CRYOSPARC_SHIFT_FIELD,
    ) -> PoseTable:
        """Normalize native CryoSPARC .cs rotvec poses into active matrices."""

        require_columns(
            particles,
            (CRYOSPARC_UID_FIELD, pose_field),
            context="CryoSPARC particles",
        )
        rows = []
        for source_row_id, row in particles.iterrows():
            source_row_id_int = int(source_row_id)
            rotvec = np.asarray(row[pose_field], dtype=float)
            if rotvec.shape != (3,):
                raise ValueError(
                    f"CryoSPARC {pose_field} at source row {source_row_id_int} "
                    f"must be a length-3 rotation vector, got shape {rotvec.shape}"
                )
            shift_xy = (
                np.asarray(row[shift_field], dtype=float)
                if shift_field in particles.columns
                else None
            )
            if shift_xy is not None and shift_xy.shape != (2,):
                raise ValueError(
                    f"CryoSPARC {shift_field} at source row {source_row_id_int} "
                    f"must be a length-2 shift vector, got shape {shift_xy.shape}"
                )
            rows.append(
                {
                    "particle_key": None,
                    "domain_name": domain_name,
                    "rotation_matrix_active": Rotation.from_rotvec(rotvec).as_matrix(),
                    "shift_xy": shift_xy,
                    "source_type": "cryosparc",
                    "source_row_id": source_row_id_int,
                    "uid": row[CRYOSPARC_UID_FIELD],
                    "raw_metadata": {
                        **row.to_dict(),
                        "_cryorole_source_row_id": source_row_id_int,
                    },
                }
            )
        return PoseTable(pd.DataFrame(rows))
