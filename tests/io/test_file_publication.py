"""Interrupted file-set writes restore old bytes, including optional CS UID lists."""
from pathlib import Path

import numpy as np
import pytest

from cryorole.export.cs_subset import write_cryosparc_cs_subset
from cryorole.io.writers import publication as module


@pytest.mark.parametrize('point,existing', [(p,e) for p in ('first_backup_after','first_install_after','second_backup_after','second_install_before','second_install_after')
                                         for e in (False,True) if 'backup' not in p or e])
def test_file_set_interruption_rolls_back_every_destination(tmp_path, monkeypatch, point, existing):
    outputs = {'metadata': tmp_path / 'main' / 'subset.cs', 'uids': tmp_path / 'lists' / 'uids.txt'}
    for path in outputs.values():
        path.parent.mkdir()
        if existing:
            path.write_bytes(b'old ' + path.name.encode())
    previous = {name: path.read_bytes() for name, path in outputs.items()} if existing else {}
    replace, failed = module.os.replace, False

    def interrupt(source, destination):
        nonlocal failed
        source, destination = Path(source), Path(destination)
        chosen = outputs['metadata' if point.startswith('first') else 'uids']
        hit = ('backup' in point and source == chosen) or ('install' in point and destination == chosen)
        if hit and not failed:
            failed = True
            if point.endswith('after'):
                replace(source, destination)
            raise KeyboardInterrupt
        replace(source, destination)

    monkeypatch.setattr(module.os, 'replace', interrupt)
    with pytest.raises(KeyboardInterrupt):
        with module.ArtifactFiles(outputs, overwrite=existing) as publication:
            for name, path in publication.paths.items():
                path.write_bytes(b'new ' + name.encode())
    for name, path in outputs.items():
        if existing:
            assert path.read_bytes() == previous[name]
        else:
            assert not path.exists()
        assert not list(path.parent.glob('.*'))


def test_failure_while_preparing_uid_keeps_existing_files(tmp_path, monkeypatch):
    source, output, uid = (tmp_path / name for name in ('input.cs', 'out.cs', 'uids.txt'))
    array = np.array([(123,)], dtype=[('uid', '<u8')])
    with source.open('wb') as stream:
        np.save(stream, array)
    output.write_bytes(b'old cs')
    uid.write_text('old uid')
    open_path = Path.open

    def fail_uid(path, *args, **kwargs):
        if path.name == 'prepared' and path.parent.name.startswith('.uids.txt.stage-'):
            raise OSError('UID writing failed')
        return open_path(path, *args, **kwargs)

    monkeypatch.setattr(Path, 'open', fail_uid)
    with pytest.raises(OSError, match='UID writing failed'):
        write_cryosparc_cs_subset(source, output, [0], row_id_field='row', overwrite=True, uid_output_path=uid)
    assert output.read_bytes() == b'old cs'
    assert uid.read_text() == 'old uid'
    assert not list(tmp_path.glob('.*'))


def test_file_set_excludes_overlapping_cooperative_writer(tmp_path):
    common = tmp_path / 'common'
    with module.ArtifactFiles({'one': common}) as first:
        with pytest.raises(FileExistsError, match='Another writer'):
            with module.ArtifactFiles({'other': tmp_path / 'other', 'common': common}):
                pytest.fail('overlapping writer entered')
        first.paths['one'].write_text('complete')
    assert common.read_text() == 'complete'
    assert not list(tmp_path.glob('.*'))


def test_failed_file_recovery_keeps_backup_and_lock(tmp_path, monkeypatch):
    output = tmp_path / 'out.cs'
    output.write_bytes(b'old')
    replace = module.os.replace

    def fail(source, destination):
        if Path(source).name in {'prepared', 'previous'}:
            raise OSError('install and recovery blocked')
        replace(source, destination)

    monkeypatch.setattr(module.os, 'replace', fail)
    with pytest.raises(RuntimeError, match='recover preserved files'):
        with module.ArtifactFiles({'main': output}, overwrite=True) as publication:
            publication.paths['main'].write_bytes(b'new')
    stage = next(tmp_path.glob('.out.cs.stage-*'))
    assert (stage / 'previous').read_bytes() == b'old'
    assert (stage / 'prepared').read_bytes() == b'new'
    assert (tmp_path / '.out.cs.file-publish').is_dir()


def test_file_cleanup_failure_preserves_new_output(tmp_path, monkeypatch):
    output = tmp_path / 'out.cs'
    output.write_bytes(b'old')
    rmtree = module.shutil.rmtree

    def fail(path, *args, **kwargs):
        if Path(path).name.startswith('.out.cs.stage-'):
            raise OSError('cleanup blocked')
        return rmtree(path, *args, **kwargs)

    monkeypatch.setattr(module.shutil, 'rmtree', fail)
    with pytest.warns(RuntimeWarning, match='Artifacts published'):
        with module.ArtifactFiles({'main': output}, overwrite=True) as publication:
            publication.paths['main'].write_bytes(b'new')
    assert output.read_bytes() == b'new'
    assert (next(tmp_path.glob('.out.cs.stage-*')) / 'previous').read_bytes() == b'old'
