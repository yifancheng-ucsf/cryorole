"""A caught publication failure must retain the previous complete artifact group."""
from pathlib import Path
import json

import pytest

from cryorole.io.writers import publication as module


@pytest.mark.parametrize('point', ['backup', 'install_before', 'install_after'])
def test_interrupted_rename_restores_old_group(tmp_path, monkeypatch, point):
    output = tmp_path / 'out'
    output.mkdir()
    (output / 'old.txt').write_text('old')
    replace = module.os.replace
    failed = False

    def interrupt(source, destination):
        nonlocal failed
        source, destination = Path(source), Path(destination)
        hit = (point == 'backup' and destination.name == 'previous') or (
            point.startswith('install') and source.name == 'artifacts')
        if hit and not failed:
            failed = True
            if point != 'install_before':
                replace(source, destination)
            raise KeyboardInterrupt
        replace(source, destination)

    monkeypatch.setattr(module.os, 'replace', interrupt)
    with pytest.raises(KeyboardInterrupt):
        with module.ArtifactDirectory(output, overwrite=True) as publication:
            (publication.path / 'new.txt').write_text('new')
    assert (output / 'old.txt').read_text() == 'old'
    assert not (output / 'new.txt').exists()
    assert not list(tmp_path.glob('.out.*'))


def test_concurrent_publication_is_excluded(tmp_path):
    output = tmp_path / 'out'
    with module.ArtifactDirectory(output) as first:
        with pytest.raises(FileExistsError, match='Another writer'):
            with module.ArtifactDirectory(output):
                pytest.fail('second writer entered')
        (first.path / 'report.json').write_text(json.dumps({'path': str(first.path / 'value.txt')}))
        (first.path / 'value.txt').write_text('value')
    assert json.loads((output / 'report.json').read_text())['path'] == str(output / 'value.txt')


def test_failed_rollback_retains_the_original_backup_and_lock(tmp_path, monkeypatch):
    output = tmp_path / 'out'
    output.mkdir()
    (output / 'old.txt').write_text('original')
    replace = module.os.replace

    def fail_install_and_restore(source, destination):
        if Path(source).name in {'artifacts', 'previous'}:
            raise OSError('cannot install or restore')
        replace(source, destination)

    monkeypatch.setattr(module.os, 'replace', fail_install_and_restore)
    with pytest.raises(RuntimeError, match='recover preserved files'):
        with module.ArtifactDirectory(output, overwrite=True) as publication:
            (publication.path / 'new.txt').write_text('new')
    stage = next(tmp_path.glob('.out.stage-*'))
    assert (stage / 'previous/old.txt').read_text() == 'original'
    assert (stage / 'artifacts/new.txt').read_text() == 'new'
    assert (tmp_path / '.out.publication-lock').is_dir()


def test_cleanup_error_after_commit_does_not_undo_new_group(tmp_path, monkeypatch):
    output = tmp_path / 'out'
    output.mkdir()
    (output / 'old.txt').write_text('original')
    rmtree = module.shutil.rmtree

    def fail_cleanup(path, *args, **kwargs):
        if Path(path).name.startswith('.out.stage-'):
            raise OSError('cannot remove staging directory')
        return rmtree(path, *args, **kwargs)

    monkeypatch.setattr(module.shutil, 'rmtree', fail_cleanup)
    with pytest.warns(RuntimeWarning, match='Artifacts published'):
        with module.ArtifactDirectory(output, overwrite=True) as publication:
            (publication.path / 'new.txt').write_text('new')
    assert (output / 'new.txt').read_text() == 'new'
    assert (output / 'old.txt').read_text() == 'original'
    stage = next(tmp_path.glob('.out.stage-*'))
    assert (stage / 'previous/old.txt').read_text() == 'original'


@pytest.mark.parametrize('already_exists', [False, True])
def test_invalid_new_report_keeps_previous_group(tmp_path, already_exists):
    output = tmp_path / 'out'
    if already_exists:
        output.mkdir()
        (output / 'old.txt').write_text('original')
    with pytest.raises(ValueError):
        with module.ArtifactDirectory(output, overwrite=already_exists) as publication:
            (publication.path / 'report.json').write_text('{incomplete')
    if already_exists:
        assert (output / 'old.txt').read_text() == 'original'
        assert not (output / 'report.json').exists()
    else:
        assert not output.exists()
    assert not list(tmp_path.glob('.out.*'))
