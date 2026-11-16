"""Real selected UIDs and source bytes must survive export failure boundaries."""
from pathlib import Path
import json
import os

import numpy as np
import pytest

from cryorole import api
from cryorole.errors import CryoroleError
from cryorole.export import metadata_subset as module
from cryorole.export.cs_subset import write_cryosparc_cs_subset
from cryorole.export.star_subset import write_relion_star_subset


def cs(path, *, shape=6, seed=8):
    data = np.zeros(shape, dtype=[('uid', '<u8'), ('alignments3D/pose', '<f8', (3,))])
    data['uid'] = np.arange(data.size).reshape(data.shape) + 100
    data['alignments3D/pose'] = np.random.default_rng(seed).normal(0, .2, (*data.shape, 3))
    with path.open('wb') as stream:
        np.save(stream, data)
    return data


def star(path):
    path.write_text('data_particles\n\nloop_\n_rlnImageName #1\n1@stack.mrcs\n2@stack.mrcs\n3@stack.mrcs\n')


def subset(kind, source, output, rows, **kwargs):
    writer = write_cryosparc_cs_subset if kind == 'cs' else write_relion_star_subset
    return writer(source, output, rows, row_id_field='source_row_id', **kwargs)


@pytest.mark.parametrize('kind', ['cs', 'star'])
@pytest.mark.parametrize('rows', [[1.2, 1.8], [True], [False, 1], [[1]], [float('nan')], ['1']])
def test_subset_rejects_noninteger_row_ids_without_output(tmp_path, kind, rows):
    source, output = tmp_path / ('input.' + kind), tmp_path / ('out.' + kind)
    (cs if kind == 'cs' else star)(source)
    before = source.read_bytes()
    with pytest.raises(ValueError, match='integer'):
        subset(kind, source, output, rows)
    assert source.read_bytes() == before
    assert not output.exists()


@pytest.mark.parametrize('kind', ['cs', 'star'])
@pytest.mark.parametrize('alias', ['self', 'symlink', 'hardlink'])
def test_subset_rejects_input_alias_as_output(tmp_path, kind, alias):
    source = tmp_path / ('input.' + kind)
    (cs if kind == 'cs' else star)(source)
    before = source.read_bytes()
    output = source if alias == 'self' else tmp_path / ('output.' + kind)
    if alias == 'symlink':
        output.symlink_to(source)
    elif alias == 'hardlink':
        os.link(source, output)
    with pytest.raises(ValueError, match='input source'):
        subset(kind, source, output, [1], overwrite=True)
    assert source.read_bytes() == before


@pytest.mark.parametrize('target', ['source', 'output', 'source_link', 'output_link'])
def test_uid_output_collision_is_rejected_before_any_write(tmp_path, target):
    source, output = tmp_path / 'input.cs', tmp_path / 'out.cs'
    cs(source)
    output.write_bytes(b'old output')
    uid = {'source': source, 'output': output}.get(target, tmp_path / 'alias.txt')
    if target.endswith('_link'):
        uid.symlink_to(source if target.startswith('source') else output)
    before_source, before_output = source.read_bytes(), output.read_bytes()
    with pytest.raises(ValueError, match='input source|distinct'):
        subset('cs', source, output, [1], overwrite=True, uid_output_path=uid)
    assert source.read_bytes() == before_source
    assert output.read_bytes() == before_output


def test_existing_uid_file_keeps_main_output_untouched(tmp_path):
    source, output, uid = tmp_path / 'input.cs', tmp_path / 'out.cs', tmp_path / 'uids.txt'
    cs(source)
    uid.write_text('previous UID list')
    with pytest.raises(FileExistsError):
        subset('cs', source, output, [1], uid_output_path=uid)
    assert not output.exists()
    assert uid.read_text() == 'previous UID list'


def test_cs_requires_one_dimensional_particle_table(tmp_path):
    source, output = tmp_path / 'input.cs', tmp_path / 'out.cs'
    cs(source, shape=(2, 3))
    with pytest.raises(ValueError, match='one-dimensional'):
        subset('cs', source, output, [1])
    assert not output.exists()


@pytest.fixture
def selected_run(tmp_path):
    ref, mov, run = tmp_path / 'ref.cs', tmp_path / 'mov.cs', tmp_path / 'run'
    original = cs(ref)
    cs(mov, seed=9)
    api.run(ref, mov, output_dir=run, k_neighbors=2, no_visualize=True, no_raw_csv=True, quiet=True)
    api.select(run, selection_id='sample', mode='random', fraction=.5, seed=12)
    return ref, original, run


@pytest.mark.parametrize('existing', [False, True])
@pytest.mark.parametrize('change_at', ['before_read', 'after_read', 'preserve_existing'])
def test_public_export_rechecks_source_before_publication(tmp_path, monkeypatch, selected_run, existing, change_at):
    ref, original, run = selected_run
    output = tmp_path / 'export'
    if existing:
        output.mkdir()
        (output / 'old.txt').write_text('previous complete output')
    before = {p.relative_to(output): p.read_bytes() for p in output.rglob('*') if p.is_file()} if existing else {}
    write_export = module._write_metadata_export
    keep = module._keep_existing_files

    def change():
        with ref.open('wb') as stream:
            np.save(stream, original[::-1])

    def staged(*args, **kwargs):
        if change_at == 'before_read':
            change()
        result = write_export(*args, **kwargs)
        if change_at == 'after_read' or (change_at == 'preserve_existing' and not existing):
            change()
        return result

    def preserve(*args):
        keep(*args)
        change()

    monkeypatch.setattr(module, '_write_metadata_export', staged)
    if change_at == 'preserve_existing':
        monkeypatch.setattr(module, '_keep_existing_files', preserve)
    with pytest.raises((CryoroleError, ValueError), match='changed|source'):
        api.export(run, selection_id='sample', domain='ref', output_dir=output, overwrite=existing)
    if existing:
        assert {p.relative_to(output): p.read_bytes() for p in output.rglob('*') if p.is_file()} == before
    else:
        assert not output.exists()
    assert not list(tmp_path.glob('.export.*'))


@pytest.mark.parametrize('point,existing', [(p,e) for p in ('backup_after','install_before','install_after')
                                         for e in (False,True) if p!='backup_after' or e])
def test_interrupted_export_retains_previous_complete_directory(tmp_path, monkeypatch, selected_run, point, existing):
    _, _, run = selected_run
    output = tmp_path / 'export'
    if existing:
        output.mkdir()
        (output / 'old.txt').write_text('previous complete output')
    replace, failed = module.os.replace, False

    def interrupt(source, destination):
        nonlocal failed
        source, destination = Path(source), Path(destination)
        hit = (point == 'backup_after' and source == output) or (point.startswith('install') and destination == output)
        if hit and not failed:
            failed = True
            if point.endswith('after'):
                replace(source, destination)
            raise KeyboardInterrupt
        replace(source, destination)

    monkeypatch.setattr(module.os, 'replace', interrupt)
    with pytest.raises(KeyboardInterrupt):
        api.export(run, selection_id='sample', domain='ref', output_dir=output, overwrite=existing)
    if existing:
        assert (output / 'old.txt').read_text() == 'previous complete output'
        assert list(output.iterdir()) == [output / 'old.txt']
    else:
        assert not output.exists()
    assert not list(tmp_path.glob('.export.*'))


@pytest.mark.parametrize('mutate', [False, True])
def test_explicit_legacy_trust_binds_current_bytes_without_claiming_historical_proof(tmp_path, monkeypatch, selected_run, mutate):
    ref, original, run = selected_run
    for name in ('run_summary.json','run_manifest.json'):
        path = run / name
        payload = json.loads(path.read_text())
        payload.pop('source_identities', None)
        path.write_text(json.dumps(payload))
    output = tmp_path / 'legacy-export'
    if mutate:
        write = module._write_metadata_export

        def change(*args, **kwargs):
            with ref.open('wb') as stream:
                np.save(stream, original[::-1])
            return write(*args, **kwargs)

        monkeypatch.setattr(module, '_write_metadata_export', change)
        with pytest.raises(CryoroleError, match='changed'):
            api.export(run, selection_id='sample', domain='ref', output_dir=output, allow_unverified_source=True)
        assert not output.exists()
    else:
        report = api.export(run, selection_id='sample', domain='ref', output_dir=output, allow_unverified_source=True)
        verification = report['source_verification']['ref']
        assert verification['verification_status'] == 'legacy_unverified_allowed'
        assert not verification['sha256_verified']
        assert len(verification['export_source_sha256']) == 64


def test_failed_export_recovery_retains_backup_and_lock(tmp_path, monkeypatch, selected_run):
    _, _, run = selected_run
    output = tmp_path / 'export'
    output.mkdir()
    (output / 'old.txt').write_text('previous complete output')
    replace = module.os.replace

    def fail(source, destination):
        source, destination = Path(source), Path(destination)
        if destination == output and source.name in {'export', 'previous'}:
            raise OSError('cannot install or restore')
        return replace(source, destination)

    monkeypatch.setattr(module.os, 'replace', fail)
    with pytest.raises(RuntimeError, match='recover preserved files'):
        api.export(run, selection_id='sample', domain='ref', output_dir=output, overwrite=True)
    # The publish lock '.export.export-publish' matches the same pattern; the staging folder is the other one.
    stage = next(path for path in tmp_path.glob('.export.export-*') if path.name != '.export.export-publish')
    assert (stage / 'previous/old.txt').read_text() == 'previous complete output'
    assert (stage / 'export/export_report.json').is_file()
    assert (tmp_path / '.export.export-publish').is_dir()
