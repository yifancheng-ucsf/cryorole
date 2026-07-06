"""Explicit particle identity resolution."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from cryorole.models.identity_report import IdentityReport
from cryorole.models.policies import IdentityPolicy


class IdentityResolutionError(ValueError):
    """Raised when particle identity cannot be resolved safely."""


@dataclass(frozen=True)
class ResolvedIdentity:
    """Resolved keys and their audit report."""

    keys: pd.Series
    report: IdentityReport


RELION_FORBIDDEN_IDENTITY_COLUMNS = {
    "_rlnAngleRot",
    "_rlnAngleTilt",
    "_rlnAnglePsi",
    "_rlnCenteredCoordinateXAngst",
    "_rlnCenteredCoordinateYAngst",
    "_rlnCenteredCoordinateZAngst",
    "_rlnCoordinateX",
    "_rlnCoordinateY",
    "_rlnCoordinateZ",
    "_rlnDefocusU",
    "_rlnDefocusV",
    "_rlnDefocusAngle",
    "_rlnOriginX",
    "_rlnOriginY",
    "_rlnOriginZ",
    "_rlnOriginXAngst",
    "_rlnOriginYAngst",
    "_rlnOriginZAngst",
}


def _normalize_value(value: Any, rule: str | None) -> str:
    if pd.isna(value):
        return ""
    normalized = str(value)
    if not rule:
        return normalized
    for step in rule.split("|"):
        if step == "strip":
            normalized = normalized.strip()
        elif step == "casefold":
            normalized = normalized.casefold()
        elif step.startswith("round:"):
            digits = int(step.split(":", 1)[1])
            normalized = str(round(float(normalized), digits))
        else:
            raise IdentityResolutionError(f"Unsupported identity normalization rule: {step}")
    return normalized


def _build_keys(table: pd.DataFrame, policy: IdentityPolicy) -> pd.Series:
    missing = [column for column in policy.identity_columns if column not in table.columns]
    if missing:
        raise IdentityResolutionError(f"Identity columns missing from input table: {missing}")

    key_parts = []
    for column in policy.identity_columns:
        rule = policy.column_normalization_rules.get(column)
        key_parts.append(table[column].map(lambda value: _normalize_value(value, rule)))

    if len(key_parts) == 1:
        return key_parts[0].rename("particle_key")

    keys = key_parts[0]
    for part in key_parts[1:]:
        keys = keys + "|" + part
    return keys.rename("particle_key")


def _build_keys_from_mapping_file(table: pd.DataFrame, policy: IdentityPolicy) -> pd.Series:
    if not policy.mapping_file:
        raise IdentityResolutionError("explicit_mapping_file identity requires mapping_file")

    mapping_path = Path(policy.mapping_file)
    if not mapping_path.is_file():
        raise IdentityResolutionError(f"Mapping file does not exist: {mapping_path}")

    mapping = pd.read_csv(mapping_path, sep=None, engine="python")
    required = {"source_row_id", "particle_key"}
    missing = required - set(mapping.columns)
    if missing:
        raise IdentityResolutionError(
            f"Mapping file missing required columns: {sorted(missing)}"
        )
    if mapping["source_row_id"].duplicated(keep=False).any():
        raise IdentityResolutionError("Mapping file contains duplicate source_row_id values")

    mapped = mapping.set_index("source_row_id")["particle_key"]
    missing_rows = [row_id for row_id in table.index if row_id not in mapped.index]
    if missing_rows:
        raise IdentityResolutionError(
            f"Mapping file does not cover source rows: {missing_rows[:5]}"
        )
    return table.index.to_series().map(mapped).map(str).rename("particle_key")


def _build_row_aligned_keys(table: pd.DataFrame) -> pd.Series:
    return table.index.to_series().map(lambda row_id: f"row:{row_id}").rename("particle_key")


def _build_relion_image_name_keys(table: pd.DataFrame) -> tuple[pd.Series, tuple[str, ...]]:
    for column in ("_rlnTomoParticleName", "_rlnImageName", "rlnImageName"):
        if column in table.columns:
            policy = IdentityPolicy.relion_user_columns((column,))
            return _build_keys(table, policy), (column,)
    raise IdentityResolutionError(
        "RELION default identity requires _rlnTomoParticleName, _rlnImageName, "
        "or rlnImageName. "
        "Use cryorole align or manually pre-align metadata before running RO computation."
    )


def _collision_examples(table: pd.DataFrame, keys: pd.Series) -> tuple[dict[str, Any], ...]:
    duplicate_mask = keys.duplicated(keep=False)
    if not duplicate_mask.any():
        return ()
    examples = []
    duplicate_rows = table.loc[duplicate_mask].copy()
    duplicate_rows["particle_key"] = keys.loc[duplicate_mask]
    for _, row in duplicate_rows.head(5).iterrows():
        examples.append(row.to_dict())
    return tuple(examples)


def resolve_identity(
    table: pd.DataFrame,
    *,
    source_type: str,
    policy: IdentityPolicy | None = None,
) -> ResolvedIdentity:
    """Resolve particle identity according to an explicit policy."""

    if source_type == "cryosparc" and policy is None:
        policy = IdentityPolicy.cryosparc_uid()

    if policy is None:
        raise IdentityResolutionError(
            "IdentityPolicy is required; RELION row-order matching is not allowed"
        )

    if source_type == "relion":
        if policy.identity_mode not in {
            "relion_image_name",
            "relion_user_columns",
            "explicit_mapping_file",
            "row_aligned",
        }:
            raise IdentityResolutionError(
                "RELION identity requires relion_image_name, relion_user_columns, row_aligned, "
                "or an explicit mapping policy"
            )
        if policy.identity_mode == "relion_user_columns":
            forbidden = RELION_FORBIDDEN_IDENTITY_COLUMNS.intersection(policy.identity_columns)
            if forbidden:
                raise IdentityResolutionError(
                    f"RELION refinement-result columns cannot be identity keys: {sorted(forbidden)}"
                )
    elif source_type == "cryosparc":
        if policy.identity_mode not in {"cryosparc_uid", "row_aligned"}:
            raise IdentityResolutionError(
                "CryoSPARC Phase 1 identity mode must be cryosparc_uid or row_aligned"
            )
    else:
        raise IdentityResolutionError(f"Unsupported source type for identity: {source_type}")

    resolved_identity_columns = policy.identity_columns
    if policy.identity_mode == "explicit_mapping_file":
        keys = _build_keys_from_mapping_file(table, policy)
    elif policy.identity_mode == "row_aligned":
        keys = _build_row_aligned_keys(table)
    elif policy.identity_mode == "relion_image_name":
        keys, resolved_identity_columns = _build_relion_image_name_keys(table)
    else:
        keys = _build_keys(table, policy)
    duplicate_count = int(keys.duplicated(keep=False).sum())
    unique_rate = float(keys.nunique(dropna=False) / len(keys)) if len(keys) else 1.0
    status = "ok"
    if duplicate_count and policy.duplicate_policy == "fail":
        status = "failed_duplicates"

    report = IdentityReport(
        identity_mode=policy.identity_mode,
        identity_columns=(
            ("source_row_id",) if policy.identity_mode == "row_aligned" else resolved_identity_columns
        ),
        column_normalization_rules=policy.column_normalization_rules,
        unique_rate=unique_rate,
        duplicate_count=duplicate_count,
        collision_examples=_collision_examples(table, keys),
        status=status,
    )

    if status.startswith("failed") and policy.failure_mode == "fail":
        raise IdentityResolutionError(f"Identity resolution failed: {status}")

    return ResolvedIdentity(keys=keys, report=report)
