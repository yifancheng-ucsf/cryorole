"""Inline comments must not be mistaken for metadata columns."""
import numpy as np
import pytest

from cryorole import api
from cryorole.export.star_subset import write_relion_star_subset
from cryorole.io.readers.star_reader import (index_star_particles, read_relion_star,
                                            split_star_fields, write_star_particle_subset)
from cryorole.io.writers.landscape_store import read_landscape_npz_arrays


def test_comment_delimiter_is_recognized_outside_quotes_only():
    assert split_star_fields('1 "#label with space" path#part # note')==['1','"#label with space"','path#part']


def test_source_reader_and_subset_keep_comment_rows_and_original_text(tmp_path):
    source=tmp_path/'source.star'
    output=tmp_path/'selected.star'
    row='1@stack.mrcs 10 20 30 # original annotation\n'
    source.write_text('data_particles\n\nloop_\n_rlnImageName #1\n_rlnAngleRot #2\n_rlnAngleTilt #3\n_rlnAnglePsi #4\n'
                      +row+'2@stack.mrcs 40 50 60\n')
    data=read_relion_star(source)
    assert data.particles['_rlnImageName'].tolist()==['1@stack.mrcs','2@stack.mrcs']
    write_relion_star_subset(source,output,[0],row_id_field='source_row')
    assert row in output.read_text()
    assert len(read_relion_star(output).particles)==1


def test_indexed_field_replacement_preserves_inline_comment_and_spacing(tmp_path):
    source=tmp_path/'source.star'
    output=tmp_path/'changed.star'
    source.write_text('data_particles\nloop_\n_rlnImageName\n_rlnAngleRot\n_rlnAngleTilt\n_rlnAnglePsi\n'
                      '1@stack.mrcs  10 20 30 # original annotation\n')
    index=index_star_particles(source)
    write_star_particle_subset(index,[0],output,replace_columns={'_rlnAngleRot':['90']})
    assert '1@stack.mrcs  90 20 30 # original annotation\n' in output.read_text()
    assert read_relion_star(output).particles._rlnAngleRot.tolist()==['90']


@pytest.mark.parametrize('backend',['array_native','dataframe_compat'])
def test_real_run_preserves_results_and_source_rows_with_inline_comments(tmp_path,backend):
    landscapes=[]
    for annotated in [False,True]:
        sources=[]
        for moving in [False,True]:
            path=tmp_path/f'{annotated}-{moving}.star'
            rows=[f'{i}@stack.mrcs {i*(3 if moving else 2)} 20 1'
                  +(' # source annotation'if annotated else '')for i in range(1,7)]
            path.write_text('data_particles\nloop_\n_rlnImageName\n_rlnAngleRot\n_rlnAngleTilt\n_rlnAnglePsi\n'
                            +'\n'.join(rows)+'\n')
            sources.append(path)
        output=tmp_path/f'run-{annotated}'
        api.run(*sources,output_dir=output,run_backend=backend,k_neighbors=2,no_visualize=True)
        landscapes.append(read_landscape_npz_arrays(output/'data/raw_landscape.npz'))
    plain,annotated=landscapes
    assert plain.particle_key.tolist()==annotated.particle_key.tolist()
    np.testing.assert_array_equal(plain.ref_source_row_id,annotated.ref_source_row_id)
    np.testing.assert_array_equal(plain.mov_source_row_id,annotated.mov_source_row_id)
    np.testing.assert_array_equal(plain.coordinates_analysis,annotated.coordinates_analysis)
    np.testing.assert_array_equal(plain.sld_raw,annotated.sld_raw)
