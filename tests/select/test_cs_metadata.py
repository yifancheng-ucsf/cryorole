"""Native CS categorical selection, precision, and source provenance."""

from dataclasses import replace
import hashlib
import json

import numpy as np
import pytest

from cryorole.cli.main import main
from cryorole.io.readers.cs_reader import read_cryosparc_cs_column
from cryorole.select.metadata import CsMetadataColumn
from cryorole.select.service import SelectRequest, create_selection


@pytest.mark.parametrize("values,requested,expected", [
    (np.array([2**64 - 1, 2**64 - 2], dtype='u8'), [str(2**64 - 1)], [True, False]),
    (np.array([True, False]), ['TRUE'], [True, False]),
    (np.array([True, False]), ['0'], [False, True]),
    (np.array(['01', '1', '', ' x ']), ['01', ' x '], [True, False, False, True]),
    (np.array([b'01', b'1', b'']), ['01'], [True, False, False]),
])
def test_typed_exact_matches(values, requested, expected):
    column = CsMetadataColumn(values)
    mask, resolved = column.select(requested)
    np.testing.assert_array_equal(mask, expected)
    assert resolved


@pytest.mark.parametrize("values,requested", [
    (np.array([1], dtype='u8'), ['-1']),
    (np.array([1], dtype='u8'), [str(2**64)]),
    (np.array([1], dtype='u8'), ['1.0']),
    (np.array([True]), ['yes']),
    (np.array(['a']), ['']),
])
def test_invalid_values_fail(values, requested):
    with pytest.raises(ValueError, match='metadata'):
        CsMetadataColumn(values).select(requested)


def test_column_reader_is_mapped_and_does_not_materialize_vector_fields(tmp_path):
    path = tmp_path / 'wide.cs'
    data = np.lib.format.open_memmap(path, mode='w+', shape=(10000,),
        dtype=[('class', 'i4'), ('unused_vector', 'f8', (256,))])
    data['class'] = np.arange(len(data)) % 3
    data.flush()
    column = read_cryosparc_cs_column(path, 'class')
    assert isinstance(column, np.memmap)
    assert column.nbytes == 40000
    assert column.strides[0] > column.dtype.itemsize
    np.testing.assert_array_equal(column[:5], [0, 1, 2, 0, 1])
    with pytest.raises(ValueError, match='scalar'):
        read_cryosparc_cs_column(path, 'unused_vector')


@pytest.fixture
def cs_run(tmp_path):
    data = np.zeros(80, dtype=[('uid', 'u8'), ('alignments3D/pose', 'f8', (3,)),
        ('alignments3D/class', 'i4'), ('flag', '?'), ('label', 'S8'), ('score', 'f4'),
        ('big_id', 'u8')])
    data['uid'] = np.arange(80)
    data['alignments3D/class'] = np.arange(80) % 3
    data['big_id'] = np.uint64(2**64 - 100) + np.arange(80, dtype='u8')
    data['label'] = np.array([b'01', b'1', b'', b'a/b', b'a_b'] * 16)
    data['flag'] = np.arange(80) % 2 == 0
    ref = tmp_path / 'ref.cs'
    mov = tmp_path / 'mov.cs'
    with ref.open('wb') as handle:
        np.save(handle, data)
    data['alignments3D/pose'] = np.random.default_rng(7).normal(0, 0.15, (80, 3))
    # Partial matching, different domain values, and reversed source row order.
    data['alignments3D/class'] = (data['alignments3D/class'] + 1) % 3
    with mov.open('wb') as handle:
        np.save(handle, data[10:][::-1])
    run = tmp_path / 'run'
    assert main(['run', '--ref', str(ref), '--mov', str(mov), '--output-dir', str(run),
        '--no-visualize', '--k-neighbors', '15']) == 0
    return run, ref, mov


def request(run, **kwargs):
    return SelectRequest(run_dir=str(run), selection_id='chosen', selection_mode='metadata',
        metadata_domain='ref', metadata_column='alignments3D/class', **kwargs)


@pytest.mark.parametrize('domain', ['ref', 'mov'])
def test_cs_selection_visualization_and_export(cs_run, domain):
    run, ref, mov = cs_run
    protected = [ref, mov, run / 'data/raw_landscape.npz']
    hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
    result = create_selection(replace(request(run, metadata_value='0,2'), metadata_domain=domain))
    summary = json.loads((result.output_dir / 'selection_summary.json').read_text())
    details = summary['selection_policy']['metadata_source_details']
    assert details['format'] == 'cryosparc_cs'
    assert details['verification_status'] == 'verified_sha256'
    assert summary['metadata_invalid_count'] == 0
    assert main(['visualize', '--run-dir', str(run), '--selection-id', 'chosen', '--all',
        '--representation', 'rotvec', '--max-points', '10']) == 0
    assert main(['export', '--run-dir', str(run), '--selection-id', 'chosen']) == 0
    parent = np.load(run / 'data/raw_landscape.npz')
    source = np.load(ref if domain == 'ref' else mov)
    matched = parent[f'{domain}_source_row_id']
    mask = np.isin(source['alignments3D/class'][matched], [0, 2])
    for export_domain, path in [('ref', ref), ('mov', mov)]:
        exported = np.load(run / 'exports/chosen' / export_domain / f'selected_{export_domain}.cs')
        expected = np.load(path)[parent[f'{export_domain}_source_row_id'][mask]]
        np.testing.assert_array_equal(exported, expected)
    assert hashes == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}


def test_split_and_conflict_prevalidation(cs_run):
    run, _, _ = cs_run
    target = run / 'selections/chosen_1'
    target.mkdir(parents=True)
    (target / 'keep.txt').write_text('keep')
    with pytest.raises(FileExistsError):
        create_selection(request(run, split_by_value=True))
    assert not (run / 'selections/chosen_0').exists()
    result = create_selection(request(run, split_by_value=True, overwrite=True))
    assert {p.name for p in result.selection_dirs} == {'chosen_0', 'chosen_1', 'chosen_2'}
    assert sum(result.selected_counts) == 70


@pytest.mark.parametrize('field', ['missing', 'score', 'alignments3D/pose'])
def test_unsupported_field_creates_no_selection(cs_run, field):
    run, _, _ = cs_run
    with pytest.raises(ValueError):
        create_selection(replace(request(run, metadata_value='1'), metadata_column=field))
    assert not list((run / 'selections').iterdir())


def test_split_collision_and_changed_source_fail_before_writing(cs_run):
    run, ref, _ = cs_run
    with pytest.raises(ValueError, match='collid'):
        create_selection(replace(request(run, split_by_value=True, overwrite=True), metadata_column='label'))
    assert not list((run / 'selections').iterdir())
    data = np.load(ref)
    data['alignments3D/class'][0] = 9
    with ref.open('wb') as handle:
        np.save(handle, data)
    with pytest.raises(ValueError, match='SHA-256'):
        create_selection(request(run, metadata_value='1'))
    assert not list((run / 'selections').iterdir())


@pytest.mark.parametrize('rows', [np.array([-1]), np.array([2]), np.array([0.5]), np.array([True])])
def test_invalid_source_row_ids_fail(rows):
    with pytest.raises(ValueError, match='source-row'):
        CsMetadataColumn.from_source(np.array([0, 1]), rows)


def test_missing_utf8_and_split_limit():
    with pytest.raises(ValueError, match='UTF-8'):
        CsMetadataColumn(np.array([b'\xff']))
    column = CsMetadataColumn(np.array(['', '01', '1', '01']))
    assert column.groups() == ('01', '1')
    assert column.missing.sum() == 1
    assert len(CsMetadataColumn(np.arange(100)).groups()) == 100
    with pytest.raises(ValueError, match='limit is 100'):
        CsMetadataColumn(np.arange(101)).groups()


@pytest.mark.parametrize('field,value', [('big_id', str(2**64 - 30)), ('flag', 'true'), ('label', '01')])
def test_cli_typed_values_and_persisted_missing_count(cs_run, field, value):
    run, ref, _ = cs_run
    assert main(['select', '--run-dir', str(run), '--selection-id', 'typed', '--mode', 'metadata',
        '--metadata-domain', 'ref', '--metadata-column', field, '--metadata-value', value]) == 0
    summary = json.loads((run / 'selections/typed/selection_summary.json').read_text())
    assert summary['metadata_values'] == [value]
    assert summary['selected_count'] == {'big_id': 1, 'flag': 35, 'label': 14}[field]
    assert summary['metadata_missing_count'] == (14 if field == 'label' else 0)
    selection = json.loads((run / 'selections/typed/selection.json').read_text())
    assert selection['active_policy']['metadata_source_details']['comparison'] == 'typed_exact'


def test_invalid_overwrite_preserves_existing_selection(cs_run):
    run, _, _ = cs_run
    result = create_selection(request(run, metadata_value='1'))
    hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in result.output_dir.iterdir() if p.is_file()}
    with pytest.raises(ValueError, match='decimal integer'):
        create_selection(replace(request(run, metadata_value='1.5', overwrite=True)))
    assert hashes == {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in hashes}



def test_cs_requires_recorded_identity(cs_run):
    run, _, _ = cs_run
    for name in ('run_summary.json', 'run_manifest.json'):
        path = run / name
        payload = json.loads(path.read_text())
        payload.pop('source_identities', None)
        path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match='recorded source SHA-256'):
        create_selection(request(run, metadata_value='1'))
    assert not list((run / 'selections').iterdir())


def test_split_recompute_prerequisites_precede_any_output(cs_run):
    run, _, _ = cs_run
    # big_id makes singleton groups, which cannot support recomputed subset SLD.
    with pytest.raises(ValueError, match='at least two'):
        create_selection(replace(request(run, split_by_value=True, write_selected_landscape=True,
            recompute_sld=True), metadata_column='big_id'))
    assert not list((run / 'selections').iterdir())


def test_metadata_help_explains_cs_types_and_limit(capsys):
    with pytest.raises(SystemExit) as result:
        main(['select', '--help'])
    assert result.value.code == 0
    help_text = capsys.readouterr().out
    assert 'CryoSPARC CS' in help_text
    assert 'scalar integers, booleans' in help_text
    assert '100 groups' in help_text
    assert 'alignments3D/class' in help_text
