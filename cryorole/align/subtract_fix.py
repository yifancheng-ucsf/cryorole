"""Write a coordinate-corrected copy of a RELION subtracted particle STAR.

``relion_particle_subtract --center_x/y/z`` moves each box to the projection
of a 3D point but writes only the sub-pixel remainder to
``_rlnOriginX/YAngst``; ``_rlnCoordinateX/Y`` keep their old values (confirmed
in ``src/particle_subtractor.cpp`` and on real data). The stored particle
position is therefore off by the whole-pixel shift, which affects any
coordinate-based step (re-extraction, polishing, distance-based duplicate
removal, coordinate matching).

This module recomputes that shift exactly from the subtraction *input* (its
angles and original origins), the centre and the pixel sizes, verifies the
prediction against the subtracted STAR (residual origins and angles must
agree), and only then writes a copy in which ``_rlnCoordinateX/Y`` are
corrected. The original files are never modified.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from cryorole.align.relion_geometry import predict_subtraction
from cryorole.align.relion_jobs import find_job_command, resolve_project_path, subtraction_parameters
from cryorole.io.readers.star_reader import (
    StarParticleIndex,
    index_star_particles,
    read_star_particle_headers,
    unquote_star_value,
    write_star_particle_subset,
)
from cryorole.align.star_values import angle_difference_deg, float_columns, optics_per_row
from cryorole.provenance.source_identity import stream_sha256

ORIGIN_TOLERANCE_ANGST = 2e-3
ANGLE_TOLERANCE_DEG = 1e-5
ACCEPT_FRACTION = 0.99
REPORT_NAME = "coordinate_correction_report.json"

_COORDS = ("_rlnCoordinateX", "_rlnCoordinateY")
_ORIGINS = ("_rlnOriginXAngst", "_rlnOriginYAngst")
_ANGLES = ("_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi")
_NEEDED = ("_rlnImageName", "_rlnImageOriginalName", "_rlnMicrographName", *_COORDS, *_ORIGINS, *_ANGLES, "_rlnOpticsGroup")


def fix_subtract_coordinates(
    target: str | Path,
    *,
    subtract_input: str | Path | None = None,
    center: Sequence[float] | None = None,
    model_angpix: float | None = None,
    micrograph_angpix: float | None = None,
    apply_to: str | Path | None = None,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    log=print,
) -> dict[str, Any]:
    """Correct the stale coordinates of a recentred RELION subtraction (see module docstring)."""

    target_path = Path(target)
    if target_path.is_dir():
        job_dir = target_path
        subtracted = job_dir / "particles_subtracted.star"
    elif target_path.is_file():
        subtracted = target_path
        job_dir = target_path.parent
    else:
        raise FileNotFoundError(f"--fix-subtract-coordinates: {target_path} is neither a Subtract job folder nor a STAR file")
    if not subtracted.is_file():
        raise FileNotFoundError(f"Subtracted STAR not found: {subtracted}")

    command = find_job_command(job_dir)
    recorded = subtraction_parameters(command) if command is not None else None
    sources: dict[str, str] = {"subtracted_star": str(subtracted)}

    if center is not None:
        center_px = [float(v) for v in center]
        sources["center"] = "--center"
    elif recorded and recorded["center_px"] is not None:
        center_px = [float(v) for v in recorded["center_px"]]
        sources["center"] = recorded["source"]
    elif recorded and recorded["recenter_on_mask"]:
        raise ValueError(
            "This subtraction recentred on the mask centre (--recenter_on_mask); cryoROLE cannot recover that "
            "point yet. Give the centre explicitly with --center X Y Z (reference-model pixels)."
        )
    else:
        raise ValueError(
            f"No recentring centre found for {job_dir} (no note.txt with --center_x/y/z). If the subtraction was "
            "recentred, give --center X Y Z (reference-model pixels); without recentring the coordinates are "
            "already consistent and need no correction."
        )

    optimiser_path: Path | None = None
    if subtract_input is not None:
        input_path = Path(subtract_input)
        sources["subtract_input"] = "--subtract-input"
    elif recorded and recorded.get("optimiser"):
        optimiser_path = resolve_project_path(recorded["optimiser"], job_dir)
        input_path = optimiser_path.with_name(optimiser_path.name.replace("_optimiser.star", "_data.star"))
        sources["subtract_input"] = f"{recorded['source']} (--i {recorded['optimiser']} → *_data.star)"
    else:
        raise ValueError("Cannot locate the subtraction input STAR; give --subtract-input Refine3D/jobNNN/run_data.star")
    if not input_path.is_file():
        raise FileNotFoundError(
            f"Subtraction input STAR not found: {input_path}. Give --subtract-input with the refinement's run_data.star."
        )

    for path in (subtracted, input_path):
        headers = read_star_particle_headers(path)
        missing = [c for c in (*_COORDS, *_ORIGINS, *_ANGLES, "_rlnMicrographName") if c not in headers]
        if missing:
            raise ValueError(f"{path} lacks columns needed for the correction: {missing}")
    sub_index = index_star_particles(subtracted, columns=_NEEDED)
    in_index = index_star_particles(input_path, columns=_NEEDED)

    a_part = optics_per_row(in_index, "_rlnImagePixelSize")
    if a_part is None:
        raise ValueError(f"{input_path}: data_optics has no _rlnImagePixelSize")
    sources["particle_angpix"] = "_rlnImagePixelSize (subtraction input optics)"
    if micrograph_angpix is not None:
        a_mic = np.full(in_index.row_count, float(micrograph_angpix))
        sources["micrograph_angpix"] = "--micrograph-angpix"
    else:
        a_mic = optics_per_row(in_index, "_rlnMicrographPixelSize")
        sources["micrograph_angpix"] = "_rlnMicrographPixelSize (input optics)"
        if a_mic is None:
            a_mic = optics_per_row(in_index, "_rlnMicrographOriginalPixelSize")
            sources["micrograph_angpix"] = (
                "_rlnMicrographOriginalPixelSize (input optics; assumes unbinned micrographs, give "
                "--micrograph-angpix otherwise)"
            )
        if a_mic is None:
            raise ValueError("No micrograph pixel size in data_optics; give --micrograph-angpix")
    if model_angpix is not None:
        a_model = np.full(in_index.row_count, float(model_angpix))
        sources["model_angpix"] = "--model-angpix"
    else:
        found = _model_pixel_size(optimiser_path)
        if found is not None:
            a_model = np.full(in_index.row_count, found[0])
            sources["model_angpix"] = found[1]
        else:
            a_model = a_part
            sources["model_angpix"] = (
                "assumed equal to the particle pixel size (no run_model.star found); verified below by the "
                "residual-origin agreement"
            )

    in_row_of_sub, mapping_rule = _map_subtracted_to_input(sub_index, in_index)
    sources["row_mapping"] = mapping_rule
    in_rows = np.asarray(in_row_of_sub, dtype=np.int64)
    prediction = predict_subtraction(
        coordinates=float_columns(in_index, _COORDS)[in_rows],
        origins_angst=float_columns(in_index, _ORIGINS)[in_rows],
        angles_deg=float_columns(in_index, _ANGLES)[in_rows],
        center_model_px=np.asarray(center_px, dtype=float),
        model_angpix=a_model[in_rows],
        particle_angpix=a_part[in_rows],
        micrograph_angpix=a_mic[in_rows],
    )
    sub_origins = float_columns(sub_index, _ORIGINS)
    origin_error = np.max(np.abs(prediction.origins_angst - sub_origins), axis=1)
    angle_error = np.max(np.abs(angle_difference_deg(float_columns(in_index, _ANGLES)[in_rows], float_columns(sub_index, _ANGLES))), axis=1)
    origin_ok = origin_error <= ORIGIN_TOLERANCE_ANGST
    angle_ok = angle_error <= ANGLE_TOLERANCE_DEG
    verified = origin_ok & angle_ok
    sub_coords = float_columns(sub_index, _COORDS)
    stale = np.all(np.isclose(sub_coords, float_columns(in_index, _COORDS)[in_rows], atol=1e-6), axis=1)
    fraction = float(verified.mean()) if verified.size else 0.0
    shift = np.linalg.norm(prediction.box_shift_mic, axis=1)
    verification = {
        "rows": int(sub_index.row_count),
        "verified_rows": int(verified.sum()),
        "verified_fraction": fraction,
        "origin_mismatch_rows": int((~origin_ok).sum()),
        "angle_mismatch_rows": int((~angle_ok).sum()),
        "max_origin_error_angst": float(origin_error[verified].max()) if verified.any() else None,
        "stale_coordinate_rows": int(stale.sum()),
        "shift_px_percentiles": {
            "p50": float(np.percentile(shift, 50)), "p99": float(np.percentile(shift, 99)), "max": float(shift.max())
        } if shift.size else {},
    }
    parameters = {
        "center_model_px": center_px,
        "model_angpix_values": sorted({float(v) for v in a_model}),
        "particle_angpix_values": sorted({float(v) for v in a_part}),
        "micrograph_angpix_values": sorted({float(v) for v in a_mic}),
        "matrix": "RELION Euler_angles2matrix(rot, tilt, psi, A, false) of the subtraction input",
        "rounding": "RELION ROUND (half away from zero), particle pixels",
        "correction": "c_corrected = c_input - ROUND(offset_particle_px) * a_part / a_mic",
    }
    log("cryoROLE subtract-coordinate correction")
    for key, value in sources.items():
        log(f"  {key}: {value}")
    log(f"  centre (model px): {center_px}")
    log(
        f"  verification: {verification['verified_rows']}/{verification['rows']} rows reproduce the subtracted "
        f"origins and angles ({fraction:.2%})"
    )
    if fraction < ACCEPT_FRACTION:
        raise ValueError(
            f"Only {fraction:.1%} of rows reproduce the subtracted origins/angles with these parameters "
            f"(need ≥ {ACCEPT_FRACTION:.0%}); origin mismatches {verification['origin_mismatch_rows']}, angle "
            f"mismatches {verification['angle_mismatch_rows']}. Check --center, --model-angpix, --subtract-input. "
            "Nothing was written."
        )

    out_dir = Path(output_dir) if output_dir is not None else _default_output_dir(job_dir)
    if out_dir.exists() and any(out_dir.iterdir()) and not overwrite:
        raise FileExistsError(f"Output directory already exists: {out_dir}. Use --overwrite or choose --output-dir.")
    out_dir.mkdir(parents=True, exist_ok=True)
    corrected = prediction.corrected_coordinates
    formatted = [[f"{value:.6f}" for value in column] for column in corrected.T]
    verified_rows = np.flatnonzero(verified)
    unverified_rows = np.flatnonzero(~verified)
    comment = f"cryoROLE subtract-coordinate correction; report: {REPORT_NAME}"
    outputs: dict[str, str] = {}
    warnings: list[str] = []
    if unverified_rows.size:
        warnings.append(f"{unverified_rows.size} rows did not verify and were left out (unverified_subtracted.star).")
        path = write_star_particle_subset(sub_index, unverified_rows.tolist(), out_dir / "unverified_subtracted.star",
                                          header_comment=comment)
        outputs["unverified_subtracted"] = str(path)
    if apply_to is None:
        output = out_dir / f"{subtracted.stem}_coords_corrected.star"
        write_star_particle_subset(
            sub_index,
            verified_rows.tolist(),
            output,
            header_comment=comment,
            replace_columns={
                "_rlnCoordinateX": [formatted[0][i] for i in verified_rows],
                "_rlnCoordinateY": [formatted[1][i] for i in verified_rows],
            },
        )
        outputs["corrected_star"] = str(output)
        applied = None
    else:
        apply_path = Path(apply_to)
        apply_index = index_star_particles(apply_path, columns=("_rlnImageName", *_COORDS))
        corrected_by_name = {
            unquote_star_value(sub_index.columns["_rlnImageName"][i]): i for i in verified_rows
        }
        rows, xs, ys, unmatched = [], [], [], []
        for row, name in enumerate(apply_index.columns["_rlnImageName"]):
            source_row = corrected_by_name.get(unquote_star_value(name))
            if source_row is None:
                unmatched.append(row)
                continue
            rows.append(row)
            xs.append(formatted[0][source_row])
            ys.append(formatted[1][source_row])
        if not rows:
            raise ValueError(f"--apply-to: no row of {apply_path} carries a subtracted _rlnImageName")
        output = out_dir / f"{apply_path.stem}_coords_corrected.star"
        write_star_particle_subset(apply_index, rows, output, header_comment=comment,
                                   replace_columns={"_rlnCoordinateX": xs, "_rlnCoordinateY": ys})
        outputs["corrected_star"] = str(output)
        if unmatched:
            warnings.append(f"--apply-to: {len(unmatched)} rows have no verified subtracted image name and were left out.")
            path = write_star_particle_subset(apply_index, unmatched, out_dir / "apply_to_unmatched.star",
                                              header_comment=comment)
            outputs["apply_to_unmatched"] = str(path)
        applied = {"path": str(apply_path), "rows": apply_index.row_count, "corrected_rows": len(rows),
                   "unmatched_rows": len(unmatched), "join": "_rlnImageName = subtracted _rlnImageName",
                   "sha256": stream_sha256(apply_path)}

    report = {
        "artifact_type": "cryorole_subtract_coordinate_correction",
        "schema_version": "1",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "subtract_job_dir": str(job_dir),
        "sources": sources,
        "parameters": parameters,
        "verification": verification,
        "apply_to": applied,
        "outputs": outputs,
        "sha256": {
            "subtracted_star": stream_sha256(subtracted),
            "subtract_input": stream_sha256(input_path),
            "corrected_star": stream_sha256(outputs["corrected_star"]),
        },
        "warnings": warnings,
    }
    report_path = out_dir / REPORT_NAME
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report["report_path"] = str(report_path)
    log(f"  wrote {outputs['corrected_star']}")
    return report


def _default_output_dir(job_dir: Path) -> Path:
    job_dir = job_dir.resolve()
    if job_dir.parent.name == "Subtract":
        return job_dir.parent.parent / "cryorole_alignments" / f"fix_subtract_{job_dir.name}"
    return job_dir / "cryorole_alignments" / "fix_subtract"


def _model_pixel_size(optimiser_path: Path | None) -> tuple[float, str] | None:
    if optimiser_path is None:
        return None
    model = optimiser_path.with_name(optimiser_path.name.replace("_optimiser.star", "_model.star"))
    if not model.is_file():
        return None
    with model.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "_rlnPixelSize":
                return float(parts[1]), f"_rlnPixelSize in {model}"
    return None


def _map_subtracted_to_input(sub: StarParticleIndex, inp: StarParticleIndex) -> tuple[list[int], str]:
    names = inp.columns.get("_rlnImageName")
    originals = sub.columns.get("_rlnImageOriginalName")
    if names and originals:
        lookup: dict[str, int] = {}
        duplicated = set()
        for row, name in enumerate(names):
            key = unquote_star_value(name)
            if key in lookup:
                duplicated.add(key)
            lookup[key] = row
        mapped = [lookup.get(unquote_star_value(value)) for value in originals]
        if all(row is not None for row in mapped) and not duplicated and len(set(mapped)) == len(mapped):
            return [int(row) for row in mapped], "subtracted _rlnImageOriginalName = input _rlnImageName (unique)"
    if sub.row_count == inp.row_count and sub.columns.get("_rlnMicrographName") == inp.columns.get("_rlnMicrographName"):
        return list(range(sub.row_count)), "row order (equal counts and identical micrograph sequence)"
    raise ValueError(
        "Cannot map subtracted particles to the subtraction input: _rlnImageOriginalName does not identify input "
        "rows uniquely and the row order differs."
    )





