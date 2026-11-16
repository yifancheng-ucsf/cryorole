"""Quoted STAR cells have the same meaning while exports retain source text."""
import numpy as np
import pytest

from cryorole import api
from cryorole.export.star_subset import write_relion_star_subset
from cryorole.io.readers.star_reader import (index_star_particles, read_relion_star,
                                            read_relion_star_particle_columns, write_star_particle_subset)
from cryorole.io.writers.landscape_store import read_landscape_npz_arrays


POSE = ("_rlnAngleRot", "_rlnAngleTilt", "_rlnAnglePsi")


def _write_star(path, *, moving=False, quote="", quote_angles=False, reverse=False, stack="stack.mrcs"):
    lines = ["data_optics", "loop_", "_rlnOpticsGroup", "_rlnVoltage", "'1' \"300\"", "",
             "data_particles", "loop_", "_rlnImageName", *POSE, "_rlnParticleName"]
    rows = []
    for i in (range(6, 0, -1) if reverse else range(1, 7)):
        image = f"{i}@{stack}"
        if quote:
            image = f"{quote}{image}{quote}"
        angles = [str(i*(3 if moving else 2)), "22" if moving else "20", "1" if moving else "0"]
        if quote_angles:
            angles = [f'"{angles[0]}"', f"'{angles[1]}'", f'"{angles[2]}"']
        rows.append(f"{image} {' '.join(angles)} '  label {i}  '")
    path.write_text("\n".join(lines + rows) + "\n", encoding="utf-8")
    return rows


def test_tables_unquote_values_while_byte_index_keeps_original_tokens(tmp_path):
    path = tmp_path / "quoted.star"
    _write_star(path, quote='"', quote_angles=True, stack="stack with spaces.mrcs")
    full = read_relion_star(path)
    compact = read_relion_star_particle_columns(path, columns=("_rlnImageName", *POSE))
    assert full.particles.loc[0, "_rlnImageName"] == "1@stack with spaces.mrcs"
    assert full.particles.loc[0, "_rlnParticleName"] == "  label 1  "
    assert full.optics.loc[0, "_rlnOpticsGroup"] == "1"
    assert full.optics.loc[0, "_rlnVoltage"] == "300"
    for column in ("_rlnImageName", *POSE):
        assert full.particles[column].tolist() == compact.particles[column].tolist()
    assert compact.particles.loc[0, "_rlnAngleRot"] == "2"
    index = index_star_particles(path, columns=("_rlnImageName", "_rlnAngleRot"))
    assert index.columns["_rlnImageName"][0] == '"1@stack with spaces.mrcs"'
    assert index.columns["_rlnAngleRot"][0] == '"2"'


@pytest.mark.parametrize("backend", ["array_native", "dataframe_compat"])
@pytest.mark.parametrize("quote_angles", [False, True])
def test_real_run_matches_quoted_particle_names_and_preserves_numeric_results(tmp_path, backend, quote_angles):
    outputs = []
    for label, quoted in [("plain", False), ("quoted", True)]:
        ref, mov = tmp_path / f"{label}_ref.star", tmp_path / f"{label}_mov.star"
        _write_star(ref, quote='"' if quoted else "", quote_angles=quote_angles and quoted)
        _write_star(mov, moving=True, quote="'" if quoted else "", quote_angles=quote_angles and quoted, reverse=True)
        run = tmp_path / label
        api.run(ref, mov, output_dir=run, run_backend=backend, k_neighbors=2, no_visualize=True)
        outputs.append(read_landscape_npz_arrays(run / "data/raw_landscape.npz"))
    plain, quoted = outputs
    assert plain.particle_key.tolist() == quoted.particle_key.tolist() == [f"{i}@stack.mrcs" for i in range(1, 7)]
    np.testing.assert_array_equal(plain.ref_source_row_id, quoted.ref_source_row_id)
    np.testing.assert_array_equal(plain.mov_source_row_id, quoted.mov_source_row_id)
    assert quoted.mov_source_row_id.tolist() == [5, 4, 3, 2, 1, 0]
    np.testing.assert_allclose(plain.coordinates_analysis, quoted.coordinates_analysis, atol=1e-14)
    np.testing.assert_allclose(plain.sld_raw, quoted.sld_raw, atol=1e-14)


@pytest.mark.parametrize("writer", ["byte_index", "metadata_export"])
def test_subset_writers_preserve_quoted_rows_and_original_source(tmp_path, writer):
    source, output = tmp_path / "source.star", tmp_path / "subset.star"
    rows = _write_star(source, quote='"', quote_angles=True, stack="stack with spaces.mrcs")
    original = source.read_bytes()
    if writer == "byte_index":
        write_star_particle_subset(index_star_particles(source), [4, 0], output)
    else:
        write_relion_star_subset(source, output, [4, 0], row_id_field="ref_source_row_id")
    assert rows[4] in output.read_text() and rows[0] in output.read_text()
    loaded = read_relion_star(output)
    assert loaded.particles._rlnImageName.tolist() == ["5@stack with spaces.mrcs", "1@stack with spaces.mrcs"]
    assert source.read_bytes() == original
