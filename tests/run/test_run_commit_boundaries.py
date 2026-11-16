"""Public run interruption must restore the exact previous bundle bytes."""
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from cryorole import api
from cryorole.run_bundle import RunBundleWriter
from cryorole.run_bundle import writer as module


def native_pair(tmp_path):
    paths = tmp_path / 'ref.cs', tmp_path / 'mov.cs'
    for path, seed in zip(paths, (2, 3)):
        data = np.zeros(6, dtype=[('uid', '<u8'), ('alignments3D/pose', '<f8', (3,))])
        data['uid'] = np.arange(100, 106)
        data['alignments3D/pose'] = np.random.default_rng(seed).normal(0, .2, (6, 3))
        with path.open('wb') as file:
            np.save(file, data)
    return paths


def files(path):
    return {str(x.relative_to(path)): sha256(x.read_bytes()).hexdigest() for x in path.rglob('*') if x.is_file()}


@pytest.mark.parametrize('existing,point', [(e,p) for e in (False, True) for p in ('backup_after','install_before','install_after')
                                         if e or p != 'backup_after'])
def test_public_run_interruption_restores_previous_bundle(tmp_path, monkeypatch, existing, point):
    ref, mov = native_pair(tmp_path)
    target = tmp_path / 'run'
    request = dict(output_dir=target, k_neighbors=2, no_visualize=True, no_raw_csv=True, quiet=True)
    if existing:
        api.run(ref, mov, **request)
    prior = files(target) if existing else {}
    source_bytes = (ref.read_bytes(), mov.read_bytes())
    replace, triggered = module.os.replace, False

    def interrupt(source, destination):
        nonlocal triggered
        source, destination = Path(source), Path(destination)
        hit = (point == 'backup_after' and source == target) or (point.startswith('install') and destination == target)
        if hit and not triggered:
            triggered = True
            if point.endswith('after'):
                replace(source, destination)
            raise KeyboardInterrupt('interrupt after filesystem operation')
        return replace(source, destination)

    monkeypatch.setattr(module.os, 'replace', interrupt)
    with pytest.raises(KeyboardInterrupt, match='filesystem operation'):
        api.run(ref, mov, overwrite=existing, **request)
    assert files(target) == prior if existing else not target.exists()
    assert (ref.read_bytes(), mov.read_bytes()) == source_bytes
    assert not list(tmp_path.glob('.run.previous-*'))
    assert not (tmp_path / '.run.publish-lock').exists()


def minimum_bundle(path):
    (path / 'data').mkdir()
    (path / 'reports').mkdir()
    np.savez(path / 'data/raw_landscape.npz', particle_key=np.asarray(['p1']))
    (path / 'data/match_table.csv').write_text('particle_key\np1\n')
    for name in ('reports/match_report.json','reports/density_report.json','run_summary.json','run_manifest.json'):
        (path / name).write_text('{}\n')
    (path / 'run_report.md').write_text('New completed run\n')


def test_failed_recovery_keeps_old_backup_prepared_stage_and_lock(tmp_path, monkeypatch):
    target = tmp_path / 'run'
    target.mkdir()
    (target / 'old.txt').write_text('old complete data')
    replace = module.os.replace
    staging = None

    def fail_install_and_restore(source, destination):
        if Path(destination) == target:
            raise OSError('install and restore blocked')
        return replace(source, destination)

    monkeypatch.setattr(module.os, 'replace', fail_install_and_restore)
    with pytest.raises(RuntimeError, match='recover preserved'):
        with RunBundleWriter(target, overwrite=True, run_id='recover') as publication:
            staging = publication.path
            minimum_bundle(staging)
            publication.commit()
    assert (tmp_path / '.run.previous-recover/old.txt').read_text() == 'old complete data'
    assert staging.is_dir()
    assert (staging / '.cryorole_bundle_complete').is_file()
    assert not (staging / 'failure_report.json').exists()
    assert not (tmp_path / '.run.failed-recover').exists()
    assert (tmp_path / '.run.publish-lock').is_dir()
