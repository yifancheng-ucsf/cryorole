"""RELION re-extraction and subtraction recentring geometry.

This module reproduces, in NumPy, how RELION moves particle coordinates and
origins when it recentres on a 3D point. It is deliberately independent of
cryoROLE's relative-orientation convention: the matrix is RELION's own
``Euler_angles2matrix(rot, tilt, psi, A, false)`` computed from the STAR angles.

Verified against real data (356,280 particles, 100 % exact) and the RELION
source (``src/preprocessing.cpp``, ``src/particle_subtractor.cpp``):

Re-extraction (``relion_preprocess --reextract_data_star … --recenter``)::

    center_part = r · (a_ref / a_part)            # --recenter_x/y/z, reference pixels → particle pixels
    off_part    = o_in / a_part − Π(A · center_part)
    off_mic     = off_part · (a_part / a_mic)     # rounding happens in MICROGRAPH pixels
    c_out       = c_in − ROUND(off_mic)
    o_out       = a_mic · (off_mic − ROUND(off_mic))

Particle subtraction (``relion_particle_subtract --center_x/y/z``)::

    center_part = r · (a_model / a_part)          # --center_x/y/z, model pixels → particle pixels
    off_part    = o_in / a_part − Π(A · center_part)
    o_sub       = a_part · (off_part − ROUND(off_part))   # rounding in PARTICLE pixels
    c_sub       = c_in                            # RELION does not update the coordinates
    shift_mic   = ROUND(off_part) · (a_part / a_mic)      # the box moved by this; c_corrected = c_in − shift_mic

``ROUND`` is RELION's rounding (half away from zero), not NumPy's
round-half-to-even.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def relion_round(values: np.ndarray) -> np.ndarray:
    """RELION ``ROUND``: round half away from zero."""

    x = np.asarray(values, dtype=float)
    return np.where(x > 0, np.floor(x + 0.5), np.ceil(x - 0.5))


def relion_euler_matrices(rot_deg, tilt_deg, psi_deg) -> np.ndarray:
    """RELION ``Euler_angles2matrix(rot, tilt, psi, A, false)`` for arrays of angles (degrees)."""

    alpha = np.radians(np.asarray(rot_deg, dtype=float))
    beta = np.radians(np.asarray(tilt_deg, dtype=float))
    gamma = np.radians(np.asarray(psi_deg, dtype=float))
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta), np.sin(beta)
    cg, sg = np.cos(gamma), np.sin(gamma)
    cc, cs, sc, ss = cb * ca, cb * sa, sb * ca, sb * sa
    matrices = np.empty(alpha.shape + (3, 3), dtype=float)
    matrices[..., 0, 0] = cg * cc - sg * sa
    matrices[..., 0, 1] = cg * cs + sg * ca
    matrices[..., 0, 2] = -cg * sb
    matrices[..., 1, 0] = -sg * cc - cg * sa
    matrices[..., 1, 1] = -sg * cs + cg * ca
    matrices[..., 1, 2] = sg * sb
    matrices[..., 2, 0] = sc
    matrices[..., 2, 1] = ss
    matrices[..., 2, 2] = cb
    return matrices


def _projected_center(matrices: np.ndarray, center_part_px: np.ndarray) -> np.ndarray:
    """Π(A · center) in particle pixels; ``center_part_px`` is (3,) or (n, 3)."""

    center = np.asarray(center_part_px, dtype=float)
    if center.ndim == 1:
        return np.einsum("nij,j->ni", matrices[:, :2, :], center)
    return np.einsum("nij,nj->ni", matrices[:, :2, :], center)


@dataclass(frozen=True)
class ExtractionPrediction:
    coordinates: np.ndarray  # (n, 2) micrograph pixels
    origins_angst: np.ndarray  # (n, 2) Å


def predict_reextraction(
    *,
    coordinates: np.ndarray,
    origins_angst: np.ndarray,
    angles_deg: np.ndarray,
    recenter_ref_px: np.ndarray,
    ref_angpix: np.ndarray | float,
    particle_angpix: np.ndarray | float,
    micrograph_angpix: np.ndarray | float,
) -> ExtractionPrediction:
    """Predict re-extraction output coordinates and origins (see module docstring)."""

    matrices = relion_euler_matrices(angles_deg[:, 0], angles_deg[:, 1], angles_deg[:, 2])
    a_part = np.asarray(particle_angpix, dtype=float).reshape(-1, 1)
    a_ref = np.asarray(ref_angpix, dtype=float).reshape(-1, 1)
    a_mic = np.asarray(micrograph_angpix, dtype=float).reshape(-1, 1)
    center_part = np.asarray(recenter_ref_px, dtype=float)[None, :] * (a_ref / a_part)
    offset_part = np.asarray(origins_angst, dtype=float) / a_part - _projected_center(matrices, center_part)
    offset_mic = offset_part * (a_part / a_mic)
    rounded = relion_round(offset_mic)
    return ExtractionPrediction(
        coordinates=np.asarray(coordinates, dtype=float) - rounded,
        origins_angst=a_mic * (offset_mic - rounded),
    )


@dataclass(frozen=True)
class SubtractionPrediction:
    origins_angst: np.ndarray  # (n, 2) Å, what RELION writes
    box_shift_mic: np.ndarray  # (n, 2) micrograph pixels the box moved (not written by RELION)
    corrected_coordinates: np.ndarray  # (n, 2) c_in − box_shift_mic


def predict_subtraction(
    *,
    coordinates: np.ndarray,
    origins_angst: np.ndarray,
    angles_deg: np.ndarray,
    center_model_px: np.ndarray,
    model_angpix: np.ndarray | float,
    particle_angpix: np.ndarray | float,
    micrograph_angpix: np.ndarray | float,
) -> SubtractionPrediction:
    """Predict subtraction origins and the coordinate shift RELION omits."""

    matrices = relion_euler_matrices(angles_deg[:, 0], angles_deg[:, 1], angles_deg[:, 2])
    a_part = np.asarray(particle_angpix, dtype=float).reshape(-1, 1)
    a_model = np.asarray(model_angpix, dtype=float).reshape(-1, 1)
    a_mic = np.asarray(micrograph_angpix, dtype=float).reshape(-1, 1)
    center_part = np.asarray(center_model_px, dtype=float)[None, :] * (a_model / a_part)
    offset_part = np.asarray(origins_angst, dtype=float) / a_part - _projected_center(matrices, center_part)
    rounded = relion_round(offset_part)
    shift_mic = rounded * (a_part / a_mic)
    return SubtractionPrediction(
        origins_angst=a_part * (offset_part - rounded),
        box_shift_mic=shift_mic,
        corrected_coordinates=np.asarray(coordinates, dtype=float) - shift_mic,
    )
