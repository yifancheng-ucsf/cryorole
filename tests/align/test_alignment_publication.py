"""Input/output collisions, source changes and interruption cannot publish false pairs."""
import json
import os
from pathlib import Path

import numpy as np
import pytest

from cryorole.align import star_align
from cryorole.io.readers.star_reader import index_star_particles, write_star_particle_subset
from cryorole.workflows.run_service import RunExecutionContext, RunRequest, execute_run
from cryorole.io.readers import star_reader
from cryorole.align import subtract_fix


def write_star(path, names=('a', 'b', 'c'), offset=0):
    path.write_text('\n'.join(['data_particles', '', 'loop_', '_rlnImageName #1', '_rlnAngleRot #2',
        '_rlnAngleTilt #3', '_rlnAnglePsi #4',
        *[f'{name} {row + offset} {row + 2} {row + 3}' for row, name in enumerate(names)]]) + '\n')


def files(directory):
    return {str(path.relative_to(directory)): path.read_bytes() for path in directory.rglob('*') if path.is_file()}


@pytest.mark.parametrize('alias', ['direct', 'symlink', 'hardlink'])
def test_star_subset_refuses_overwriting_input_alias(tmp_path, alias):
    source, target = tmp_path / 'source.star', tmp_path / 'target.star'
    write_star(source)
    original = source.read_bytes()
    index = index_star_particles(source)
    if alias == 'direct':
        target = source
    elif alias == 'symlink':
        target.symlink_to(source)
    else:
        os.link(source, target)  # Path.hardlink_to needs Python 3.10
    with pytest.raises(ValueError, match='input source'):
        write_star_particle_subset(index, [0, 1], target)
    assert source.read_bytes() == original


def test_run_refuses_replacing_a_directory_that_contains_its_inputs(tmp_path):
    out = tmp_path / 'run'
    out.mkdir()
    paths = [out / 'ref.cs', out / 'mov.cs']
    for seed, path in enumerate(paths):
        data = np.zeros(6, dtype=[('uid', '<u8'), ('alignments3D/pose', '<f8', (3,))])
        data['uid'] = np.arange(6)
        data['alignments3D/pose'] = np.random.default_rng(seed).normal(0, .2, (6, 3))
        with path.open('wb') as stream:
            np.save(stream, data)
    original = files(out)
    with pytest.raises(ValueError, match='contains an input source'):
        execute_run(RunRequest(ref=str(paths[0]), mov=str(paths[1]), output_dir=str(out), overwrite=True,
                              no_visualize=True, no_raw_csv=True, k_neighbors=2, quiet=True),
                    context=RunExecutionContext(stream_progress=False))
    assert files(out) == original


def test_alignment_refuses_output_ancestor_of_source(tmp_path):
    ref, mov = tmp_path / 'aligned_ref.star', tmp_path / 'mov.star'
    write_star(ref)
    write_star(mov, offset=1)
    original = files(tmp_path)
    with pytest.raises(ValueError, match='contains an input source'):
        star_align.align_star_files(ref=ref, mov=mov, output_dir=tmp_path, overwrite=True)
    assert files(tmp_path) == original


def test_alignment_input_change_after_plan_publishes_nothing(tmp_path, monkeypatch):
    ref, mov, out = tmp_path / 'ref.star', tmp_path / 'mov.star', tmp_path / 'aligned'
    write_star(ref)
    write_star(mov, ('c', 'b', 'a'), offset=1)
    planner = star_align._align_by_keys

    def changed(**kwargs):
        result = planner(**kwargs)
        write_star(ref, ('c', 'b', 'a'))
        return result

    monkeypatch.setattr(star_align, '_align_by_keys', changed)
    with pytest.raises(ValueError, match='source.*changed'):
        star_align.align_star_files(ref=ref, mov=mov, output_dir=out)
    assert not out.exists()


def test_interrupted_alignment_overwrite_keeps_all_old_files(tmp_path, monkeypatch):
    ref, mov, out = tmp_path / 'ref.star', tmp_path / 'mov.star', tmp_path / 'aligned'
    write_star(ref)
    write_star(mov, offset=1)
    star_align.align_star_files(ref=ref, mov=mov, output_dir=out)
    old = files(out)
    write_star(ref, offset=5)
    writer = star_align.write_star_particle_subset

    def interrupted(*args, **kwargs):
        writer(*args, **kwargs)
        raise KeyboardInterrupt

    monkeypatch.setattr(star_align, 'write_star_particle_subset', interrupted)
    with pytest.raises(KeyboardInterrupt):
        star_align.align_star_files(ref=ref, mov=mov, output_dir=out, overwrite=True)
    assert files(out) == old


def test_successful_alignment_keeps_unrelated_files_and_public_paths(tmp_path):
    ref, mov, out = tmp_path / 'ref.star', tmp_path / 'mov.star', tmp_path / 'aligned'
    write_star(ref)
    write_star(mov, ('c', 'b', 'a'), offset=1)
    out.mkdir()
    (out / 'notes.txt').write_text('retain this note')
    original = {path: path.read_bytes() for path in [ref, mov]}
    report = star_align.align_star_files(ref=ref, mov=mov, output_dir=out, overwrite=True)
    assert report['output_dir'] == str(out)
    assert Path(report['output_paths']['align_report']).is_file()
    assert (out / 'notes.txt').read_text() == 'retain this note'
    assert all(path.read_bytes() == data for path, data in original.items())
    saved = json.loads((out / 'align_report.json').read_text())
    assert saved == report
    assert '.stage-' not in saved['next_command']


@pytest.mark.parametrize('rows', [[-1], [3], [1.0], [True], [[0]], ['0']])
def test_invalid_subset_rows_preserve_existing_output(tmp_path, rows):
    source, target = tmp_path / 'source.star', tmp_path / 'target.star'
    write_star(source)
    target.write_bytes(b'existing artifact')
    original = source.read_bytes()
    with pytest.raises(ValueError, match='integer source rows'):
        write_star_particle_subset(index_star_particles(source), rows, target)
    assert source.read_bytes() == original
    assert target.read_bytes() == b'existing artifact'


def test_interrupted_field_replacement_preserves_existing_output(tmp_path, monkeypatch):
    source, target = tmp_path / 'source.star', tmp_path / 'target.star'
    write_star(source)
    target.write_bytes(b'existing artifact')

    def fail(*args, **kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(star_reader, '_replace_fields_in_place', fail)
    with pytest.raises(KeyboardInterrupt):
        write_star_particle_subset(index_star_particles(source), [0, 1], target,
                                    replace_columns={'_rlnAngleRot': ['4', '5']})
    assert target.read_bytes() == b'existing artifact'
    assert not list(tmp_path.glob('.target.star.*.tmp'))


def test_preexisting_temporary_name_is_never_removed(tmp_path, monkeypatch):
    source, target = tmp_path / 'source.star', tmp_path / 'target.star'
    write_star(source)
    collision = tmp_path / '.target.star.collision.tmp'
    collision.write_text('another writer owns this file')
    monkeypatch.setattr(star_reader, 'uuid4', lambda: type('Identifier', (), {'hex': 'collision'})())
    with pytest.raises(FileExistsError):
        write_star_particle_subset(index_star_particles(source), [0], target)
    assert collision.read_text() == 'another writer owns this file'
    assert not target.exists()


def subtraction_pair(root):
    """Identity angles, zero origins and centre (3,0,0): corrected x is x+3."""
    root.mkdir(parents=True, exist_ok=True)
    source, subtracted = root / 'input.star', root / 'subtracted.star'
    common = ['data_optics', '', 'loop_', '_rlnOpticsGroup #1', '_rlnImagePixelSize #2',
              '_rlnMicrographOriginalPixelSize #3', '1 1 1', '', 'data_particles', '', 'loop_']
    headers = ['_rlnImageName', '_rlnMicrographName', '_rlnCoordinateX', '_rlnCoordinateY',
               '_rlnAngleRot', '_rlnAngleTilt', '_rlnAnglePsi', '_rlnOriginXAngst',
               '_rlnOriginYAngst', '_rlnOpticsGroup', '_rlnImageOriginalName']
    for path, stack in [(source, 'raw'), (subtracted, 'sub')]:
        body = [f'{i+1}@{stack}.mrcs mic.mrc {100+i} {200+i} 0 0 0 0 0 1 {i+1}@raw.mrcs' for i in range(6)]
        path.write_text('\n'.join([*common, *[f'{name} #{i+1}' for i, name in enumerate(headers)], *body]) + '\n')
    return source, subtracted


def correct(source, subtracted, output, **kwargs):
    return subtract_fix.fix_subtract_coordinates(subtracted, subtract_input=source,
        center=[3, 0, 0], model_angpix=1, micrograph_angpix=1, output_dir=output, **kwargs)


def test_subtract_success_reports_committed_paths_and_preserves_inputs(tmp_path):
    source, subtracted = subtraction_pair(tmp_path / 'inputs')
    original = files(source.parent)
    messages = []
    output = tmp_path / 'out'
    report = correct(source, subtracted, output, log=messages.append)
    result = star_reader.read_relion_star(report['outputs']['corrected_star']).particles
    np.testing.assert_array_equal(result['_rlnCoordinateX'].astype(float), np.arange(100, 106) + 3)
    assert files(source.parent) == original
    assert Path(report['report_path']).is_file()
    assert '.stage-' not in messages[-1]
    assert messages[-1].strip() == 'wrote ' + report['outputs']['corrected_star']


def test_subtract_interrupt_preserves_previous_group(tmp_path, monkeypatch):
    source, subtracted = subtraction_pair(tmp_path / 'inputs')
    output = tmp_path / 'out'
    correct(source, subtracted, output, log=lambda *_: None)
    original = files(output)
    write = subtract_fix.write_star_particle_subset

    def interrupt(*args, **kwargs):
        write(*args, **kwargs)
        raise KeyboardInterrupt

    monkeypatch.setattr(subtract_fix, 'write_star_particle_subset', interrupt)
    with pytest.raises(KeyboardInterrupt):
        correct(source, subtracted, output, overwrite=True, log=lambda *_: None)
    assert files(output) == original


def test_subtract_input_change_during_write_publishes_nothing(tmp_path, monkeypatch):
    source, subtracted = subtraction_pair(tmp_path / 'inputs')
    output = tmp_path / 'out'
    write = subtract_fix.write_star_particle_subset

    def change(*args, **kwargs):
        result = write(*args, **kwargs)
        source.write_text(source.read_text() + '# another process changed the input\n')
        return result

    monkeypatch.setattr(subtract_fix, 'write_star_particle_subset', change)
    with pytest.raises(ValueError, match='source.*changed'):
        correct(source, subtracted, output, log=lambda *_: None)
    assert not output.exists()


def test_subtract_discovered_input_is_protected_from_output_replacement(tmp_path):
    source, subtracted = subtraction_pair(tmp_path / 'original')
    project = tmp_path / 'project'
    refine, subtract = project / 'Refine3D/job001', project / 'Subtract/job002'
    refine.mkdir(parents=True)
    subtract.mkdir(parents=True)
    discovered = refine / 'run_data.star'
    discovered.write_bytes(source.read_bytes())
    target = subtract / 'particles_subtracted.star'
    target.write_bytes(subtracted.read_bytes())
    (subtract / 'note.txt').write_text('relion_particle_subtract --i Refine3D/job001/run_optimiser.star '
        '--center_x 3 --center_y 0 --center_z 0\n')
    original = files(project)
    with pytest.raises(ValueError, match='contains an input source'):
        subtract_fix.fix_subtract_coordinates(subtract, output_dir=refine, overwrite=True, log=lambda *_: None)
    assert files(project) == original
