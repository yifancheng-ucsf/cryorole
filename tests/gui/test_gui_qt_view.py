"""View page on a real (small) run bundle, headless: 2D and 1D views, no OpenGL."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole import api
from cryorole.gui.app import MainWindow
from cryorole.gui.glprobe import ProbeResult
from cryorole.gui.pages.view import ViewPage

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run"))
from test_core_productionization import _write_cs  # noqa: E402

N = 300
NO_GL = ProbeResult(platform_ok=True, gl_status="failed", platform_name="offscreen",
                    reason="the Qt platform 'offscreen' has no OpenGL widgets", remedy="2D and 1D work")


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    root = tmp_path_factory.mktemp("gui_view")
    rng = np.random.default_rng(5)
    uids = list(range(100, 100 + N))
    ref = rng.normal(0.0, 0.5, (N, 3))
    mov = (Rotation.from_rotvec(ref) * Rotation.from_rotvec(rng.normal(0.0, 0.1, (N, 3)))).as_rotvec()
    _write_cs(root / "ref.cs", uids, ref)
    _write_cs(root / "mov.cs", uids, mov)
    api.run(root / "ref.cs", root / "mov.cs", output_dir=root / "run", k_neighbors=10, no_visualize=True)
    api.canonicalize(root / "run", no_visualize=True)
    return root / "run"


@pytest.fixture
def page(qtbot, run_dir):
    widget = ViewPage(gl_available=False, gl_reason="no OpenGL here")
    qtbot.addWidget(widget)
    widget.resize(1200, 700)
    widget.show()
    widget.set_landscape(api.open_landscape(run_dir))
    return widget


def _drawn_pixels(page) -> list[int]:
    return [int(np.count_nonzero(image.image[..., 3])) for image in page.view2d.images]


def test_without_opengl_the_3d_window_is_unavailable_and_explained(page):
    assert page.view3d is None
    assert [page.tabs.tabText(i) for i in range(page.tabs.count())] == ["2D", "1D"]
    assert not page.open3d_button.isEnabled()
    assert "no OpenGL here" in page.open3d_button.toolTip()
    page.open_3d()  # ignored without OpenGL
    assert page.view3d is None
    assert "3D:" not in page.status.text()


def test_loaded_landscape_fills_all_2d_panels_and_1d_histograms(page):
    assert all(count > 0 for count in _drawn_pixels(page))
    for plot in (*page.view1d.coordinate_plots, page.view1d.sld_plot):
        assert plot.listDataItems()
    assert page.status.text().endswith("display only")
    shown = page.frame.filtered_count
    assert f"2D and 1D: all {shown:,} particles with SLD >= 1 · full landscape {N:,}" in page.status.text()
    assert [page.space_combo.itemText(i) for i in range(page.space_combo.count())] == ["Raw", "Canonical: default"]


def test_default_filter_is_sld_at_least_1_like_cryorole_visualize(page):
    assert page.filter_combo.currentData() == "threshold" and page.filter_value.value() == 1.0
    sld = page.landscape.sld_display
    assert page.frame.filtered_count == int(np.count_nonzero(sld >= 1.0))
    assert "particles with SLD >= 1 ·" in page.status.text()


def test_top_percent_and_all_change_counts_but_not_the_landscape(page):
    page.filter_combo.setCurrentIndex(page.filter_combo.findData("top_percent"))
    assert page.filter_value.value() == 40.0 and page.filter_value.suffix() == " %"
    page.filter_value.setValue(10.0)
    page.refresh(reset_view=False)
    assert page.frame.filtered_count == int(np.ceil(N * 0.1))  # continuous SLD: no ties here
    assert "top 10% by SLD" in page.status.text()
    page.filter_combo.setCurrentIndex(page.filter_combo.findData("all"))
    assert page.frame.filtered_count == N and not page.filter_value.isEnabled()
    assert page.landscape.n_points == N


def _view_range(page, panel):
    return page.view2d.plots[panel].getViewBox().viewRange()


def _assert_range(actual, expected, tolerance=0.02):
    span = expected[1] - expected[0]
    assert abs(actual[0] - expected[0]) <= tolerance * span and abs(actual[1] - expected[1]) <= tolerance * span


def _settle(qtbot):
    qtbot.wait(50)


def test_standard_axis_ranges_with_equal_units_are_shown_exactly(page, qtbot):
    _settle(qtbot)
    assert page.axis_mode == "standard" and page.equal_units.isChecked()
    assert page.axis_ranges == [(-180.0, 180.0)] * 3
    for panel in range(3):
        x_range, y_range = _view_range(page, panel)
        _assert_range(x_range, (-180.0, 180.0))
        _assert_range(y_range, (-180.0, 180.0))  # not stretched past the real angle range


def test_custom_ranges_drive_2d_and_1d_and_keep_equal_units(page, qtbot):
    page.set_custom_ranges([(-40.0, 40.0), (-30.0, 30.0), (-50.0, 50.0)])
    _settle(qtbot)
    assert page.axis_mode == "custom" and page.custom_box.isVisibleTo(page)
    x_range, y_range = _view_range(page, 0)  # alpha vs beta
    _assert_range(x_range, (-40.0, 40.0))
    _assert_range(y_range, (-30.0, 30.0))
    assert page.view2d.boxes[0].aspect == pytest.approx(80.0 / 60.0)
    for plot, expected in zip(page.view1d.coordinate_plots, page.axis_ranges):
        _assert_range(plot.getViewBox().viewRange()[0], expected)
    page.equal_units.setChecked(False)
    assert page.view2d.boxes[0].aspect is None


def test_invalid_custom_range_is_reported_and_the_view_is_kept(page, qtbot):
    page.set_custom_ranges([(-40.0, 40.0), (-30.0, 30.0), (-50.0, 50.0)])
    before = list(page.axis_ranges)
    low, high = page.custom_spins[1]
    low.setValue(10.0)
    high.setValue(5.0)
    page.refresh(reset_view=True)
    assert page.axis_ranges == before
    assert "lower < upper for: beta" in page.status.text()


def test_rotation_vector_uses_its_own_standard_range(page):
    page.representation_combo.setCurrentIndex(page.representation_combo.findData("rotvec"))
    assert page.axis_ranges == [(-3.14, 3.14)] * 3
    assert [label.text() for label in page.custom_labels] == ["x", "y", "z"]


def test_fit_modes_frame_the_data_without_filtering(page):
    page.axis_combo.setCurrentIndex(page.axis_combo.findData("core"))
    core = page.axis_ranges
    page.axis_combo.setCurrentIndex(page.axis_combo.findData("all"))
    every = page.axis_ranges
    assert all(e[0] <= c[0] and e[1] >= c[1] for c, e in zip(core, every))
    assert page.frame.filtered_count == int(np.count_nonzero(page.landscape.sld_display >= 1.0))


def test_switching_to_canonical_loads_in_the_background(page, qtbot):
    with qtbot.waitSignal(page.landscape_loaded, timeout=10000):
        page.space_combo.setCurrentIndex(1)
        page.space_combo.activated.emit(1)
    assert page.landscape.space == "canonical"
    assert "(canonical default)" in page.status.text()
    # The rebuilt Space box must name the loaded frame, not fall back to "Raw".
    assert page.space_combo.currentText() == "Canonical: default"
    assert page.space_combo.currentData() == ("canonical", "default")


def test_main_window_opens_a_run_and_reports_a_wrong_folder(qtbot, run_dir, tmp_path):
    window = MainWindow(NO_GL)
    qtbot.addWidget(window)
    with qtbot.waitSignal(window.view_page.landscape_loaded, timeout=10000):
        window.open_run(run_dir)
    assert window.pages.currentWidget() is window.view_page
    assert str(run_dir) in window.windowTitle()
    window.open_run(tmp_path)
    qtbot.waitUntil(lambda: window.view_page.findChild(object, "error_input_invalid") is not None, timeout=10000)
    assert window.view_page.isEnabled()
    assert window.view_page.landscape.run_dir == Path(run_dir).resolve()  # the loaded run is kept


def test_newer_landscape_request_supersedes_pending_load(page, qtbot):
    loaded = []
    page.landscape_loaded.connect(lambda landscape: loaded.append(landscape.space))
    with qtbot.waitSignal(page.landscape_loaded, timeout=5000):
        page.open_run(page.landscape.run_dir)
        page.open_run(page.landscape.run_dir, space="canonical", canonical_id="default")
    assert loaded == ["canonical"]
    assert page.isEnabled()


def test_canonicalize_does_not_override_a_newer_open(page, qtbot, monkeypatch):
    monkeypatch.setattr(page._canonicalizer, "start", lambda request: None)
    run_dir = page.landscape.run_dir
    page.canonicalize(canonical_id="old-job", fit_top_fraction=0.5, visualize=False)
    with qtbot.waitSignal(page.landscape_loaded, timeout=5000):
        page.open_run(run_dir)
    page._canonicalized({"output_dir": str(run_dir / "canonical" / "old-job")})
    assert page.isEnabled()
    assert not page._load_runner.busy


def _two_runs(root):
    def make(name, n, spread, seed):
        rng = np.random.default_rng(seed)
        ref = rng.normal(0.0, 0.6, (n, 3))
        mov = (Rotation.from_rotvec(ref) * Rotation.from_rotvec(rng.normal(0.0, spread, (n, 3)))).as_rotvec()
        _write_cs(root / f"{name}_ref.cs", list(range(n)), ref)
        _write_cs(root / f"{name}_mov.cs", list(range(n)), mov)
        api.run(root / f"{name}_ref.cs", root / f"{name}_mov.cs", output_dir=root / name, k_neighbors=10,
                no_visualize=True)
        return root / name

    return make("A", 400, 0.02, 1), make("B", 150, 0.9, 2)  # A: high SLD; B: spread out, low SLD


def test_a_run_whose_display_is_empty_never_shows_the_previous_run(qtbot, tmp_path):
    a, b = _two_runs(tmp_path)
    max_a, max_b = (float(api.open_landscape(x).sld_display.max()) for x in (a, b))
    assert max_a > max_b
    page = ViewPage(gl_available=False)
    qtbot.addWidget(page)
    page.resize(1200, 800)
    page.show()
    page.open_run(a)
    qtbot.waitUntil(lambda: page.landscape is not None, timeout=30000)
    page.filter_value.setValue((max_a + max_b) / 2)  # keeps a few of A, none of B
    page.refresh(reset_view=False)
    assert page.frame is not None and page.frame.full_count == 400 and page.stack.currentIndex() == 1
    page.open_run(b)
    qtbot.waitUntil(lambda: page.landscape is not None and page.landscape.run_dir.name == "B", timeout=30000)
    assert page.frame is None and page.stack.currentIndex() == 0  # the message, not A's pictures
    assert "400" not in page.status.text()
    page.picking = True
    page._pick_2d(0, 0.0, 0.0)  # was an IndexError on A's rows
    page.filter_value.setValue(0.0)  # a filter that keeps B's particles brings the views back
    page.refresh(reset_view=True)
    assert page.frame.full_count == 150 and page.stack.currentIndex() == 1


def test_canonicalize_result_does_not_replace_another_run(qtbot, tmp_path):
    a, b = _two_runs(tmp_path)
    page = ViewPage(gl_available=False)
    qtbot.addWidget(page)
    page.open_run(a)
    qtbot.waitUntil(lambda: page.landscape is not None, timeout=30000)
    page._canonicalize_run_dir, page._pending_frame, page._canonicalize_command = a, "f1", "cryorole canonicalize"
    page._canonicalize_generation = page._load_generation
    page.open_run(b)  # the user opens another run while A is canonicalized
    qtbot.waitUntil(lambda: page.landscape.run_dir.name == "B", timeout=30000)
    page._canonicalized({"output_dir": str(a / "canonical" / "f1")})
    qtbot.wait(200)
    assert page.landscape.run_dir.name == "B" and "f1 of A is ready" in page.status.text()


def test_quitting_during_canonicalize_asks_and_removes_the_unfinished_frame(qtbot, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QMessageBox

    a, _ = _two_runs(tmp_path)
    window = MainWindow(NO_GL)
    qtbot.addWidget(window)
    window.show()
    with qtbot.waitSignal(window.view_page.landscape_loaded, timeout=10000):
        window.open_run(a)
    page = window.view_page
    page.canonicalize(canonical_id="late", fit_top_fraction=0.2, visualize=False)
    (a / "canonical" / "late").mkdir(parents=True, exist_ok=True)  # as if the worker had begun writing
    (a / "canonical" / "late" / "canonical_landscape.npz").write_bytes(b"partial")
    answers = []

    def ask(*_args):
        answers.append(_args[-1])
        return QMessageBox.StandardButton.No if len(answers) == 1 else QMessageBox.StandardButton.Yes

    monkeypatch.setattr(QMessageBox, "question", ask)
    assert page.canonicalizing
    window.close()  # "No": the window stays, the worker keeps going
    assert window.isVisible() and page.canonicalizing
    window.close()  # "Yes": stopped, and the frame it started is removed
    assert not window.isVisible() and not page.canonicalizing
    assert not (a / "canonical" / "late").exists()
    assert len(answers) == 2 and "canonical frame is being written" in answers[0]
