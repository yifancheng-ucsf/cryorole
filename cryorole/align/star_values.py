"""Small typed accessors for indexed STAR particle columns."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from cryorole.io.readers.star_reader import StarParticleIndex, unquote_star_value

OPTICS_GROUP_COLUMN = "_rlnOpticsGroup"


def optics_per_row(index: StarParticleIndex, name: str) -> np.ndarray | None:
    """Per-particle value of optics column ``name`` (via ``_rlnOpticsGroup``), or ``None`` if absent."""

    values = index.optics.get(name)
    if not values:
        return None
    numeric = np.asarray([float(unquote_star_value(v)) for v in values], dtype=float)
    groups = index.optics.get(OPTICS_GROUP_COLUMN)
    if OPTICS_GROUP_COLUMN in index.columns and groups:
        decoded = [unquote_star_value(g) for g in groups]
        if len(set(decoded)) != len(decoded):
            raise ValueError(f"{index.path}: duplicate optics group IDs in data_optics")
        lookup = {g: numeric[i] for i, g in enumerate(decoded)}
        try:
            return np.asarray([lookup[unquote_star_value(g)] for g in index.columns[OPTICS_GROUP_COLUMN]], dtype=float)
        except KeyError as exc:
            raise ValueError(f"{index.path}: particle optics group {exc} is not in data_optics") from exc
    if len(numeric) == 1:
        return np.full(index.row_count, numeric[0])
    raise ValueError(f"{index.path}: several optics groups but no {OPTICS_GROUP_COLUMN} column")


def float_columns(index: StarParticleIndex, names: Sequence[str]) -> np.ndarray:
    """(n, len(names)) float array of retained particle columns."""

    return np.column_stack([
        np.asarray([float(unquote_star_value(v)) for v in index.columns[name]], dtype=float) for name in names
    ])


def angle_difference_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Signed difference wrapped to [-180, 180)."""

    return (np.asarray(a) - np.asarray(b) + 180.0) % 360.0 - 180.0
