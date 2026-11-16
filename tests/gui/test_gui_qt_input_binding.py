"""Input labels, actual run sources and readiness must follow the latest choice."""
from pathlib import Path
from threading import Event
from types import SimpleNamespace

import numpy as np
import pytest
from PySide6.QtCore import QThreadPool, QTimer
from shiboken6 import isValid

from cryorole import api
from cryorole.gui.pages.run_setup import InputsPage, RunPage, RunSetup


@pytest.fixture
def inputs(tmp_path):
    paths = [tmp_path / (name + '.cs') for name in ('old', 'new', 'moving')]
    for path in paths:
        data = np.zeros(6, dtype=[('uid', '<u8'), ('alignments3D/pose', '<f8', (3,))])
        data['uid'] = np.arange(100, 106)
        with path.open('wb') as stream:
            np.save(stream, data)
    return [(path, api.inspect_input(path)) for path in paths]


def finish_threads(qtbot):
    qtbot.waitUntil(lambda: QThreadPool.globalInstance().activeThreadCount() == 0, timeout=5000)
    qtbot.wait(30)  # deliver completion signals already queued by the worker


@pytest.mark.parametrize('old_outcome', ['success', 'failure'])
def test_two_input_choices_cannot_restore_the_older_file(qtbot, monkeypatch, tmp_path, inputs, old_outcome):
    setup = RunSetup()
    setup.output_dir = tmp_path / 'run'
    page = InputsPage(setup)
    qtbot.addWidget(page)
    box = page.boxes['ref']
    old, new = inputs[:2]
    started, release = Event(), Event()

    def inspect(path):
        if Path(path) == old[0]:
            started.set()
            assert release.wait(3)
            if old_outcome == 'failure':
                raise ValueError('old inspection failed')
            return old[1]
        return new[1]

    monkeypatch.setattr(api, 'inspect_input', inspect)
    box.set_path(old[0])
    qtbot.waitUntil(started.is_set)
    box.set_path(new[0])
    QTimer.singleShot(100, release.set)
    finish_threads(qtbot)
    assert box.path_edit.text() == str(new[0])
    assert setup.chosen['ref'] == new[0]
    assert setup.inspections['ref'].chosen.path == new[0]


def test_choosing_a_new_file_immediately_invalidates_the_old_source(qtbot, monkeypatch, tmp_path, inputs):
    setup = RunSetup()
    setup.output_dir = tmp_path / 'run'
    page, run = InputsPage(setup), RunPage(setup)
    qtbot.addWidget(page)
    qtbot.addWidget(run)
    old, new, mov = inputs
    setup.set_inspection('ref', old[1])
    setup.set_inspection('mov', mov[1])
    setup._timer.stop()
    setup.preflight = {'readiness': 'READY'}
    started, release = Event(), Event()

    def inspect(path):
        started.set()
        assert release.wait(3)
        return new[1]

    monkeypatch.setattr(api, 'inspect_input', inspect)
    try:
        page.boxes['ref'].set_path(new[0])
        qtbot.waitUntil(started.is_set)
        assert setup.chosen['ref'] is None
        assert not setup.inputs_ready
        assert setup.preflight is None
        assert not run.can_start()[0]
    finally:
        release.set()
        finish_threads(qtbot)
        setup._timer.stop()


@pytest.mark.parametrize('change', ['input_removed', 'options_changed'])
def test_old_preflight_cannot_become_ready_between_new_requests(qtbot, monkeypatch, tmp_path, inputs, change):
    setup = RunSetup()
    setup.output_dir = tmp_path / 'run'
    page = InputsPage(setup)
    qtbot.addWidget(page)
    setup.set_inspection('ref', inputs[0][1])
    setup.set_inspection('mov', inputs[2][1])
    setup._timer.stop()
    started, release = Event(), Event()

    def preflight(*args, **kwargs):
        started.set()
        assert release.wait(3)
        return SimpleNamespace(report={'readiness': 'READY', 'request': 'old'})

    monkeypatch.setattr(api, 'preflight', preflight)
    setup._run_preflight()
    qtbot.waitUntil(started.is_set)
    if change == 'input_removed':
        setup.set_inspection('ref', None)
    else:
        setup.k_neighbors = 7
        setup.options_changed.emit()
    setup._timer.stop()  # hold the next debounced request until old completion arrives
    release.set()
    finish_threads(qtbot)
    assert setup.preflight is None
    assert setup.readiness is None
    assert setup.preflight_error is None


def test_incomplete_inputs_cancel_a_pending_preflight_timer(qtbot, monkeypatch, inputs):
    setup = RunSetup()
    page = InputsPage(setup)
    qtbot.addWidget(page)
    setup.set_inspection('ref', inputs[0][1])
    setup.set_inspection('mov', inputs[2][1])
    assert setup._timer.isActive()
    setup.set_inspection('ref', None)
    try:
        assert not setup._timer.isActive()
        assert not setup.preflight_running
    finally:
        setup._timer.stop()


def test_deleted_input_page_discards_its_pending_callback(qtbot, monkeypatch, inputs):
    setup = RunSetup()
    page = InputsPage(setup)
    started, release = Event(), Event()

    def inspect(path):
        started.set()
        assert release.wait(3)
        return inputs[0][1]

    monkeypatch.setattr(api, 'inspect_input', inspect)
    page.boxes['ref'].set_path(inputs[0][0])
    qtbot.waitUntil(started.is_set)
    try:
        page.deleteLater()
        qtbot.waitUntil(lambda: not isValid(page))
    finally:
        release.set()
        finish_threads(qtbot)
        if isValid(page):
            page.deleteLater()
            qtbot.waitUntil(lambda: not isValid(page))
    assert setup.chosen['ref'] is None
