"""Alignment numeric policies and transfer identities must be explicit and valid."""
import numpy as np
import pytest

from cryorole.align.keys import KeySpec, normalize_float_tolerances, normalize_key_value
from cryorole.align.star_align import align_star_files
from cryorole.align.subtract_fix import fix_subtract_coordinates
from cryorole.io.readers.star_reader import read_relion_star
from test_alignment_publication import subtraction_pair, correct


@pytest.mark.parametrize('tolerance', [float('nan'), float('inf'), float('-inf'), 0, -1])
def test_float_tolerance_must_be_finite_and_positive(tolerance):
    with pytest.raises(ValueError, match='finite.*> 0'):
        normalize_float_tolerances({'_coord': tolerance})


@pytest.mark.parametrize('raw,tolerance', [('nan', .1), ('inf', .1), ('-inf', .1), ('1e308', 1e-300)])
def test_invalid_float_key_value_reports_a_value_error(raw, tolerance):
    spec = KeySpec(ref_columns=('_coord',), mov_columns=('_coord',), float_tolerances={'coord': tolerance})
    with pytest.raises(ValueError, match='finite'):
        normalize_key_value(raw, '_coord', spec)


@pytest.mark.parametrize('side', ['ref', 'mov'])
def test_cross_column_tolerance_on_either_side_applies_to_the_pair(tmp_path, side):
    ref, mov = tmp_path / 'ref.star', tmp_path / 'mov.star'
    ref.write_text('data_particles\n\nloop_\n_refCoord #1\n10.0\n')
    mov.write_text('data_particles\n\nloop_\n_movCoord #1\n10.04\n')
    report = align_star_files(ref=ref, mov=mov, key_pairs=['_refCoord=_movCoord'],
        float_tolerances={f'_{side}Coord': .1}, output_dir=tmp_path/'out')
    assert report['matched_count'] == 1
    assert report['key_policy']['resolved_float_tolerances'] == {'refcoord': .1, 'movcoord': .1}


def test_cross_column_conflicting_scales_cannot_act_as_unit_conversion(tmp_path):
    ref, mov = tmp_path/'ref.star', tmp_path/'mov.star'
    ref.write_text('data_particles\n\nloop_\n_refCoord #1\n10.0\n')
    mov.write_text('data_particles\n\nloop_\n_movCoord #1\n20.0\n')
    with pytest.raises(ValueError, match='same tolerance'):
        align_star_files(ref=ref, mov=mov, key_pairs=['_refCoord=_movCoord'],
            float_tolerances={'_refCoord':10., '_movCoord':20.}, output_dir=tmp_path/'out')
    assert not (tmp_path/'out').exists()


def test_cross_column_tolerance_normalization_is_transitive_and_does_not_mutate_input():
    supplied = {'A': .1}
    spec = KeySpec(ref_columns=('A','B'), mov_columns=('B','C'), float_tolerances=supplied)
    assert spec.float_tolerances == {'a':.1,'b':.1,'c':.1}
    assert supplied == {'A':.1}
    with pytest.raises(ValueError, match='same tolerance'):
        KeySpec(ref_columns=('A','B'),mov_columns=('B','C'),float_tolerances={'A':.1,'C':.2})


@pytest.mark.parametrize('style', ['numeric', 'groups'])
def test_quoted_optics_are_decoded_as_values_and_source_bytes_are_preserved(tmp_path, style):
    source, subtracted = subtraction_pair(tmp_path/'inputs')
    for path in [source,subtracted]:
        replacement='1 "1" \'1\'' if style=='numeric' else '"1" 1 1'
        text=path.read_text().replace('\n1 1 1\n','\n'+replacement+'\n')
        if style=='groups':
            text=text.replace(' 0 0 0 0 0 1 ',' 0 0 0 0 0 \'1\' ')
        path.write_text(text)
    original={path:path.read_bytes()for path in [source,subtracted]}
    report=correct(source,subtracted,tmp_path/'out',log=lambda *_:None)
    table=read_relion_star(report['outputs']['corrected_star']).particles
    np.testing.assert_array_equal(table['_rlnCoordinateX'].astype(float),np.arange(103,109))
    assert all(path.read_bytes()==data for path,data in original.items())


def test_duplicate_optics_group_ids_are_rejected(tmp_path):
    source,subtracted=subtraction_pair(tmp_path/'inputs')
    source.write_text(source.read_text().replace('\n1 1 1\n','\n1 1 1\n1 1 1\n'))
    with pytest.raises(ValueError,match='duplicate optics group'):
        correct(source,subtracted,tmp_path/'out',log=lambda *_:None)
    assert not (tmp_path/'out').exists()


def test_apply_to_rejects_duplicate_verified_source_names(tmp_path):
    source,subtracted=subtraction_pair(tmp_path/'inputs')
    text=subtracted.read_text()
    for number in range(2,7):
        text=text.replace(f'{number}@sub.mrcs','1@sub.mrcs')
    subtracted.write_text(text)
    target=tmp_path/'downstream.star'
    target.write_text('data_particles\n\nloop_\n_rlnImageName #1\n_rlnCoordinateX #2\n_rlnCoordinateY #3\n1@sub.mrcs 105 205\n')
    original={path:path.read_bytes()for path in [source,subtracted,target]}
    with pytest.raises(ValueError,match='verified.*not unique'):
        fix_subtract_coordinates(subtracted,subtract_input=source,center=[3,0,0],model_angpix=1,
            micrograph_angpix=1,apply_to=target,output_dir=tmp_path/'out',log=lambda *_:None)
    assert not (tmp_path/'out').exists()
    assert all(path.read_bytes()==data for path,data in original.items())
