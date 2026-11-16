"""Public draft payloads are owned by callers; confirmations retain exact science."""
from hashlib import sha256
import shlex
import subprocess
import sys

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole import api
from cryorole.export import read_selection_json


@pytest.fixture
def native_view(tmp_path):
    ref, mov = tmp_path / 'ref.cs', tmp_path / 'mov.cs'
    for path, seed in ((ref, 0), (mov, 1)):
        data = np.zeros(6, dtype=[('uid', '<u8'), ('alignments3D/pose', '<f8', (3,))])
        data['uid'] = np.arange(100, 106)
        data['alignments3D/pose'] = np.random.default_rng(seed).normal(0, .3, (6, 3))
        with path.open('wb') as stream:
            np.save(stream, data)
    run = tmp_path / 'SELECTION_ID run'
    api.run(ref, mov, output_dir=run, k_neighbors=2, no_visualize=True, no_raw_csv=True, quiet=True)
    return api.open_landscape(run)


@pytest.mark.parametrize('change', ['all_true', 'all_false', 'invert'])
def test_caller_mask_changes_cannot_change_confirmed_particle_keys(native_view, change):
    view = native_view
    draft = view.evaluate_radius(view.coordinates('rotvec', [0])[0], representation='rotvec', radius_deg=1.)
    expected = view.session.arrays.particle_key[draft.selected_mask].tolist()
    assert len(expected) == draft.selected_count == 1
    raw = view.run_dir / 'data/raw_landscape.npz'
    before = sha256(raw.read_bytes()).hexdigest()
    if change == 'invert':
        np.logical_not(draft.selected_mask, out=draft.selected_mask)
    else:
        draft.selected_mask[:] = change == 'all_true'
    view.confirm('kept', draft)
    selected = read_selection_json(view.run_dir / 'selections/kept/selection.json')
    assert list(selected.selected_particle_keys) == list(map(str, expected))
    assert selected.selected_count == draft.selected_count
    assert sha256(raw.read_bytes()).hexdigest() == before


def test_returned_report_center_edits_do_not_change_confirmed_geometry(native_view):
    view = native_view
    center = view.coordinates('rotvec', [0])[0]
    report = view.session.evaluate(center=list(center), representation='rotvec', radius_deg=1.)
    report['center_rv'][:] = [np.pi, 0., 0.]
    report['center_input'][:] = [np.pi, 0., 0.]
    report['center_euler'][:] = [180., 0., 0.]
    view.session.confirm('kept')
    selected = read_selection_json(view.run_dir / 'selections/kept/selection.json')
    assert np.allclose(selected.center_evaluated, center)
    assert np.allclose(selected.center_input, center)
    rows = np.flatnonzero(np.isin(view.session.arrays.particle_key, selected.selected_particle_keys))
    distances = (Rotation.from_rotvec(center).inv() * Rotation.from_rotvec(view.coordinates('rotvec', rows))).magnitude()
    assert np.all(distances <= np.deg2rad(selected.radius) + 1e-10)


@pytest.mark.parametrize('payload', ['session', 'draft'])
def test_payload_display_filter_edits_cannot_rewrite_confirmation_provenance(native_view, payload):
    view = native_view
    draft = view.evaluate_radius(view.coordinates('rotvec', [0])[0], representation='rotvec', radius_deg=1.)
    if payload == 'session':
        report = view.session.session_payload()
    else:
        report = view.session.evaluate(center=list(draft.center_rv), representation='rotvec', radius_deg=1.)
    report['display_filter']['mode'] = 'caller-edited'
    view.session.confirm('kept')
    selected = read_selection_json(view.run_dir / 'selections/kept/selection.json')
    assert selected.interaction_provenance['display_filter'] == {'mode': 'all', 'field': 'sld_display'}


def test_exact_gui_command_preserves_placeholder_words_in_run_path(native_view):
    view = native_view
    draft = view.evaluate_radius(view.coordinates('rotvec', [0])[0], representation='rotvec', radius_deg=1.)
    parts = shlex.split(draft.command('cli-copy'))
    assert parts[parts.index('--run-dir') + 1] == str(view.run_dir)
    assert parts[parts.index('--selection-id') + 1] == 'cli-copy'
    completed = subprocess.run([sys.executable, '-B', '-m', 'cryorole.cli.main', *parts[1:]], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    view.confirm('gui-copy', draft)
    root = view.run_dir / 'selections'
    cli = read_selection_json(root / 'cli-copy/selection.json')
    gui = read_selection_json(root / 'gui-copy/selection.json')
    assert cli.selected_particle_keys == gui.selected_particle_keys
