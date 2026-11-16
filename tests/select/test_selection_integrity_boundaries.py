"""Public selection entries preserve exact identities and complete output groups."""
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cryorole import api
from cryorole.export import export_selection
from cryorole.io.readers.star_reader import read_relion_star
from cryorole.io.writers.landscape_store import read_landscape_npz_arrays
from cryorole.models.landscape import Landscape
from cryorole.models.policies import SelectionExportPolicy, SelectionPolicy
from cryorole.select import select_particles
from cryorole.select import artifacts
from cryorole.select.artifacts import (
    SelectedRowProvenance, SelectionArtifactRequest, write_selection_artifact,
)


LABELS = [2**53, 2**53 + 1, 2**64 - 1, 2**64 - 2, 10, 11]


def _write_star(path, *, moving=False):
    header = ['data_particles', 'loop_', '_rlnImageName', '_rlnAngleRot',
              '_rlnAngleTilt', '_rlnAnglePsi', '_rlnParticleId']
    rows = [f'{i+1}@particles.mrcs {i*(3 if moving else 2)} '
            f'{22 if moving else 20} {i if moving else 0} {label}'
            for i, label in enumerate(LABELS)]
    path.write_text('\n'.join(header + (rows[::-1] if moving else rows)) + '\n')


@pytest.fixture
def real_run(tmp_path):
    def make(backend='array_native'):
        ref, mov, run = tmp_path / 'ref.star', tmp_path / 'mov.star', tmp_path / 'run'
        _write_star(ref)
        _write_star(mov, moving=True)
        api.run(ref, mov, output_dir=run, run_backend=backend,
                k_neighbors=2, no_visualize=True, no_raw_csv=True, quiet=True)
        return run, ref, mov
    return make


@pytest.mark.parametrize('backend', ['array_native', 'dataframe_compat'])
@pytest.mark.parametrize('domain', ['ref', 'mov'])
@pytest.mark.parametrize('label', [2**53 + 1, 2**64 - 1, 2**64 - 2])
def test_real_star_selection_and_export_keep_adjacent_large_ids(real_run, backend, domain, label):
    run, ref, mov = real_run(backend)
    protected = [ref, mov, run / 'data/raw_landscape.npz']
    hashes = {p: sha256(p.read_bytes()).hexdigest() for p in protected}
    result = api.select(run, selection_id='exact', mode='metadata', metadata_domain=domain,
                        metadata_column='rlnParticleId', metadata_value=str(label))
    assert result.selected_counts == (1,)
    api.export(run, selection_id='exact')
    for selected_domain in ['ref', 'mov']:
        path = run / 'exports/exact' / selected_domain / f'selected_{selected_domain}.star'
        table = read_relion_star(path).particles
        assert table['_rlnParticleId'].tolist() == [str(label)]
        assert table['_rlnImageName'].tolist() == [f'{LABELS.index(label)+1}@particles.mrcs']
    assert hashes == {p: sha256(p.read_bytes()).hexdigest() for p in protected}


@pytest.mark.parametrize('backend', ['array_native', 'dataframe_compat'])
@pytest.mark.parametrize('domain', ['ref', 'mov'])
def test_real_star_split_keeps_every_large_integer_category(real_run, backend, domain):
    run, _, _ = real_run(backend)
    result = api.select(run, selection_id='split', mode='metadata', metadata_domain=domain,
                        metadata_column='rlnParticleId', split_by_value=True)
    assert len(result.selection_dirs) == 6
    assert result.selected_counts == (1,) * 6
    assert {p.name for p in result.selection_dirs} == {f'split_{label}' for label in LABELS}


def _landscape():
    return Landscape(data=pd.DataFrame({
        'particle_key': ['p0', 'p1', 'p2'],
        'coordinates_analysis': [np.zeros(3)] * 3,
        'coordinates_display': [np.zeros(3)] * 3,
        'sld_raw': [1., 3., 5.],
        'sld_unfloored': [1., 3., 5.],
        'sld_display': [1., 3., 5.],
        'sld_was_floored': [False] * 3,
        'sld_local_k_mean': [1.] * 3,
        'sld_effective_local_k_mean': [1.] * 3,
        'sld_distance_floor': [1e-4] * 3,
    }))


@pytest.mark.parametrize('bounds,expected', [
    ({'sld_min': 2.}, ('p1', 'p2')),
    ({'sld_max': 3.}, ('p0', 'p1')),
    ({'sld_min': 2., 'sld_max': 4.}, ('p1',)),
])
def test_library_density_bounds_get_a_default_name(bounds, expected):
    selection = select_particles(_landscape(), policy=SelectionPolicy(
        selection_mode='threshold_by_density', **bounds))
    assert selection.selected_particle_keys == expected
    assert selection.selection_id.startswith('selection:threshold_by_density:sld_raw:')
    assert 'None' not in selection.selection_id


def test_library_threshold_default_name_stays_compatible():
    selection = select_particles(_landscape(), policy=SelectionPolicy(
        selection_mode='threshold_by_density', threshold=2.))
    assert selection.selection_id == 'selection:threshold_by_density:sld_raw:>=2'


@pytest.mark.parametrize('values,query,expected', [
    (['01', '1.0', '1e0'], '1', ('p0', 'p1', 'p2')),
    (['9007199254740992.0', '9007199254740993.0', '2'],
     '9007199254740993', ('p1',)),
    (['1.000000000000000000001', '1', '2'], '1', ('p1',)),
])
def test_library_numeric_metadata_canonicalization_is_exact(values, query, expected):
    landscape = _landscape()
    landscape.data['ref_source_row_id'] = [0, 1, 2]
    selected = select_particles(landscape, policy=SelectionPolicy(
        selection_mode='metadata_value', selection_id='numeric', metadata_domain='ref',
        metadata_column='label', metadata_values=(query,)),
        source_metadata=pd.DataFrame({'label': values}))
    assert selected.selected_particle_keys == expected


@pytest.mark.parametrize('row_type', [int, str, np.uint64])
@pytest.mark.parametrize('metadata_kind', ['dataframe', 'mapping'])
def test_library_custom_source_row_labels_keep_integer_precision(row_type, metadata_kind):
    landscape = _landscape()
    rows = [2**53, 2**53 + 1, 2**64 - 1]
    landscape.data['ref_source_row_id'] = [row_type(v) for v in rows]
    values = ['a', 'b', 'c']
    metadata = (pd.DataFrame({'label': values}, index=rows)
                if metadata_kind == 'dataframe' else dict(zip(rows, values)))
    selected = select_particles(landscape, policy=SelectionPolicy(
        selection_mode='metadata_value', selection_id='custom', metadata_domain='ref',
        metadata_column='label', metadata_values=('b',)), source_metadata=metadata)
    assert selected.selected_particle_keys == ('p1',)
    assert selected.metadata_missing_count == 0


def _tree(path):
    return {str(p.relative_to(path)): p.read_bytes() for p in path.rglob('*') if p.is_file()}


def _assert_no_staging(parent):
    assert not [p for p in parent.iterdir() if p.name.startswith('.')]


@pytest.mark.parametrize('existing', [False, True])
@pytest.mark.parametrize('derived', [False, True])
@pytest.mark.parametrize('fail_at', ['rows_csv', 'summary'])
def test_api_selection_failure_preserves_complete_previous_group(
    real_run, monkeypatch, existing, derived, fail_at,
):
    run, _, _ = real_run()
    output = run / 'selections/chosen'
    if existing:
        api.select(run, selection_id='chosen', center=(0, 0, 0), radius=180,
                   write_selected_landscape=derived)
        (output / 'user-notes.txt').write_text('keep this unrelated file')
    before = _tree(output)

    def fail(*args, **kwargs):
        raise OSError('injected selection output failure')

    monkeypatch.setattr(artifacts, '_write_selected_rows' if fail_at == 'rows_csv'
                        else 'write_json_artifact', fail)
    with pytest.raises(api.CryoroleError, match='injected selection output failure'):
        api.select(run, selection_id='chosen', mode='random', fraction=.5, seed=23,
                   write_selected_landscape=derived, overwrite=existing)
    assert _tree(output) == before
    assert output.exists() == existing
    _assert_no_staging(output.parent)


@pytest.mark.parametrize('existing', [False, True])
def test_api_keyboard_interrupt_keeps_previous_selection(real_run, monkeypatch, existing):
    run, _, _ = real_run()
    output = run / 'selections/chosen'
    if existing:
        api.select(run, selection_id='chosen', center=(0, 0, 0), radius=180,
                   write_selected_landscape=True)
    before = _tree(output)

    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt('injected selection interruption')

    monkeypatch.setattr(artifacts, 'write_json_artifact', interrupt)
    with pytest.raises(KeyboardInterrupt, match='injected selection interruption'):
        api.select(run, selection_id='chosen', mode='random', fraction=.5, seed=23,
                   write_selected_landscape=True, overwrite=existing)
    assert _tree(output) == before
    assert output.exists() == existing
    _assert_no_staging(output.parent)


def test_api_overwrite_replaces_selected_landscape_and_public_paths(real_run):
    run, _, _ = real_run()
    output = api.select(run, selection_id='chosen', center=(0, 0, 0), radius=180,
                        write_selected_landscape=True).output_dir
    (output / 'user-notes.txt').write_text('keep me')
    api.select(run, selection_id='chosen', mode='random', fraction=.5, seed=23,
               write_selected_landscape=True, overwrite=True)
    payload = json.loads((output / 'selection.json').read_text())
    derived = read_landscape_npz_arrays(output / 'selected_landscape/landscape.npz')
    assert derived.particle_key.tolist() == payload['selected_particle_keys']
    assert derived.n_points == 3
    assert (output / 'user-notes.txt').read_text() == 'keep me'
    for path in output.rglob('*.json'):
        assert '.stage-' not in path.read_text()


def test_api_overwrite_without_derived_output_removes_stale_landscape(real_run):
    run, _, _ = real_run()
    output = api.select(run, selection_id='chosen', center=(0, 0, 0), radius=180,
                        write_selected_landscape=True).output_dir
    (output / 'user-notes.txt').write_text('keep me')
    api.select(run, selection_id='chosen', mode='random', fraction=.5, seed=23,
               overwrite=True)
    assert not (output / 'selected_landscape').exists()
    assert (output / 'user-notes.txt').read_text() == 'keep me'


def _artifact_request(output, *, overwrite=False):
    selection = select_particles(_landscape(), policy=SelectionPolicy(
        top_fraction=1., selection_id='writer'))
    return SelectionArtifactRequest(selection=selection, output_dir=output,
        overwrite=overwrite, selected_rows=SelectedRowProvenance(selection.selected_particle_keys))


@pytest.mark.parametrize('existing', [False, True])
def test_standard_writer_failure_is_one_group(tmp_path, monkeypatch, existing):
    output = tmp_path / 'selection'
    request = _artifact_request(output)
    if existing:
        write_selection_artifact(request)
    before = _tree(output)

    def fail(*args, **kwargs):
        raise OSError('injected final summary failure')

    monkeypatch.setattr(artifacts, 'write_json_artifact', fail)
    changed = replace(request.selection, selected_particle_keys=('p1',), selected_count=1)
    with pytest.raises(OSError, match='injected final summary failure'):
        write_selection_artifact(replace(request, overwrite=existing, selection=changed,
            selected_rows=SelectedRowProvenance(changed.selected_particle_keys)))
    assert _tree(output) == before
    assert output.exists() == existing
    _assert_no_staging(output.parent)


def test_standard_writer_returns_published_paths(tmp_path):
    output = tmp_path / 'selection'
    result = write_selection_artifact(_artifact_request(output))
    assert result.output_dir == output
    for path in [result.selection_json, result.selected_particle_keys_csv,
                 result.selection_csv, result.selected_landscape_rows_csv,
                 result.selection_summary_json]:
        assert path.is_file() and path.parent == output
    for path in result.export_report.output_paths.values():
        assert Path(path).is_file() and '.stage-' not in path
    assert '.stage-' not in json.dumps(dict(result.export_report.export_policy))


@pytest.mark.parametrize('existing', [False, True])
def test_low_level_export_late_failure_keeps_all_previous_files(tmp_path, monkeypatch, existing):
    from cryorole.export import selection_export as module
    output = tmp_path / 'export'
    request = _artifact_request(output)
    policy = SelectionExportPolicy(output_dir=output)
    if existing:
        export_selection(request.selection, policy=policy)
    before = _tree(output)

    def fail(*args, **kwargs):
        raise OSError('injected export report failure')

    monkeypatch.setattr(module, '_write_export_report_json', fail)
    changed = replace(request.selection, selected_particle_keys=('p1',), selected_count=1)
    with pytest.raises(OSError, match='injected export report failure'):
        export_selection(changed, policy=replace(policy, overwrite=existing))
    assert _tree(output) == before
    _assert_no_staging(output.parent)


@pytest.mark.parametrize('left,right', [
    ('selected_keys_filename', 'selection_json_filename'),
    ('selected_keys_filename', 'export_report_filename'),
    ('selection_json_filename', 'export_report_filename'),
])
def test_low_level_export_rejects_colliding_filenames_before_writing(tmp_path, left, right):
    request = _artifact_request(tmp_path / 'selection')
    policy = SelectionExportPolicy(output_dir=tmp_path / 'export')
    policy = replace(policy, **{left: 'same.txt', right: 'same.txt'})
    with pytest.raises(ValueError, match='distinct'):
        export_selection(request.selection, policy=policy)
    assert not (tmp_path / 'export').exists()
