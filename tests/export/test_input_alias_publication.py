"""Directory publication preserves recorded input names as well as source bytes."""
import json
import os

import numpy as np
import pytest

from cryorole import api
from cryorole.errors import CryoroleError
from cryorole.export import metadata_subset as module
from test_metadata_export_boundaries import cs


def snapshot(folder):
    return {str(p.relative_to(folder)): ('link', os.readlink(p)) if p.is_symlink()
            else ('file', p.read_bytes()) for p in folder.rglob('*') if p.is_file() or p.is_symlink()}


def prepare(tmp_path, kind='direct', relocated=False, source_domain='ref'):
    ref, mov, run, output = (tmp_path/name for name in ('ref.cs', 'mov.cs', 'run', 'export'))
    cs(ref, shape=12)
    cs(mov, shape=12, seed=9)
    (output/source_domain).mkdir(parents=True)
    (output/'previous.txt').write_text('previous output')
    target = ref if source_domain == 'ref' else mov
    name = 'input.cs' if kind == 'other-name' else 'selected_'+source_domain+'.cs'
    alias = output/source_domain/name
    alias.symlink_to(target)
    used = alias
    if kind == 'parent-alias':
        view = tmp_path/'view'
        view.symlink_to(output, target_is_directory=True)
        used = view/source_domain/name
    elif kind == 'nested-parent':
        real = tmp_path/'source-dir'
        real.mkdir()
        moved = real/target.name
        target.rename(moved)
        target.symlink_to(moved)
        alias.unlink()
        (output/'bridge').symlink_to(real, target_is_directory=True)
        used = output/'bridge'/target.name
        alias = output/'bridge'
    options = {source_domain: used} if not relocated else {}
    api.run(options.get('ref', ref), options.get('mov', mov), output_dir=run,
            k_neighbors=2, no_visualize=True, no_raw_csv=True, quiet=True)
    api.select(run, selection_id='sample', mode='random', fraction=.5, seed=12)
    return ref, mov, run, output, alias, used


@pytest.mark.parametrize('kind', ['direct', 'other-name', 'parent-alias', 'nested-parent'])
@pytest.mark.parametrize('relocated', [False, True])
def test_public_export_refuses_original_or_relocated_input_alias(tmp_path, kind, relocated):
    ref, mov, run, output, alias, used = prepare(tmp_path, kind, relocated)
    before = snapshot(output), ref.read_bytes(), mov.read_bytes(), os.readlink(alias)
    options = {'relocated_ref': used} if relocated else {}
    with pytest.raises((ValueError, CryoroleError), match='input source|separate|source metadata'):
        api.export(run, selection_id='sample', output_dir=output, overwrite=True, **options)
    assert (snapshot(output), ref.read_bytes(), mov.read_bytes(), os.readlink(alias)) == before
    assert not list(tmp_path.glob('.export.*'))


def test_unrequested_domain_input_name_is_still_protected(tmp_path):
    ref, mov, run, output, alias, _ = prepare(tmp_path, source_domain='mov')
    before = snapshot(output), ref.read_bytes(), mov.read_bytes()
    with pytest.raises((ValueError, CryoroleError), match='input source|separate|source metadata'):
        api.export(run, selection_id='sample', domain='ref', output_dir=output, overwrite=True)
    assert alias.is_symlink() and (snapshot(output), ref.read_bytes(), mov.read_bytes()) == before


def test_relative_input_uses_saved_absolute_path_after_chdir(tmp_path, monkeypatch):
    ref, mov, run, output, alias, _ = prepare(tmp_path)
    # A new normal run records a relative spelling and its absolute runtime base.
    monkeypatch.chdir(tmp_path)
    run2 = tmp_path/'relative-run'
    api.run(alias.relative_to(tmp_path), mov, output_dir=run2, k_neighbors=2,
            no_visualize=True, no_raw_csv=True, quiet=True)
    api.select(run2, selection_id='sample', mode='random', fraction=.5, seed=12)
    monkeypatch.chdir(run)
    with pytest.raises((ValueError, CryoroleError), match='input source|separate|source metadata'):
        api.export(run2, selection_id='sample', output_dir=output, overwrite=True)
    assert alias.is_symlink() and len(np.load(ref)) == 12


@pytest.mark.parametrize('relocated', [False, True])
def test_separate_alias_and_normal_overwrite_are_legal(tmp_path, relocated):
    ref, mov = tmp_path/'ref.cs', tmp_path/'mov.cs'
    cs(ref, shape=12)
    cs(mov, shape=12, seed=9)
    alias = tmp_path/'input-link.cs'
    alias.symlink_to(ref)
    view = tmp_path/'parent-link'
    view.symlink_to(tmp_path, target_is_directory=True)
    run, output = tmp_path/'run', tmp_path/'export'
    api.run(alias, mov, output_dir=run, k_neighbors=2, no_visualize=True, quiet=True)
    api.select(run, selection_id='sample', mode='random', fraction=.5, seed=12)
    before = ref.read_bytes(), mov.read_bytes(), os.readlink(alias)
    options = {'relocated_ref': view/alias.name} if relocated else {}
    for overwrite in (False, True):
        report = api.export(run, selection_id='sample', output_dir=output, overwrite=overwrite, **options)
        assert len(np.load(output/'ref/selected_ref.cs')) == 6
        assert json.loads((output/'ref/export_ref_report.json').read_text())['source_file_unchanged']
        assert report['selected_count'] == 6
    assert before == (ref.read_bytes(), mov.read_bytes(), os.readlink(alias))


def test_relocated_alias_retarget_during_preparation_preserves_previous_output(tmp_path, monkeypatch):
    ref, mov = tmp_path/'ref.cs', tmp_path/'mov.cs'
    cs(ref, shape=12)
    cs(mov, shape=12, seed=9)
    alias = tmp_path/'relocated.cs'
    alias.symlink_to(ref)
    run, output = tmp_path/'run', tmp_path/'export'
    api.run(ref, mov, output_dir=run, k_neighbors=2, no_visualize=True, quiet=True)
    api.select(run, selection_id='sample', mode='random', fraction=.5, seed=12)
    api.export(run, selection_id='sample', output_dir=output)
    before = snapshot(output), ref.read_bytes(), mov.read_bytes()
    write = module._write_metadata_export

    def retarget(*args, **kwargs):
        result = write(*args, **kwargs)
        alias.unlink()
        alias.symlink_to(mov)
        return result

    monkeypatch.setattr(module, '_write_metadata_export', retarget)
    with pytest.raises((ValueError, CryoroleError), match='changed|source'):
        api.export(run, selection_id='sample', output_dir=output, overwrite=True, relocated_ref=alias)
    assert before == (snapshot(output), ref.read_bytes(), mov.read_bytes())
