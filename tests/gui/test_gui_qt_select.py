"""Select page, headless: pick, live exact count, highlight, Confirm, and a form mode."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cryorole import api
from cryorole.export import read_selection_json
from cryorole.gui.app import MainWindow, PAGES
from cryorole.gui.glprobe import ProbeResult

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "run"))
from test_core_productionization import _write_cs  # noqa: E402

N = 300
NO_GL = ProbeResult(platform_ok=True, gl_status="failed", platform_name="offscreen",
                    reason="no OpenGL", remedy="2D and 1D work")


@pytest.fixture
def run_dir(tmp_path):
    rng = np.random.default_rng(9)
    uids = list(range(500, 500 + N))
    ref = rng.normal(0.0, 0.4, (N, 3))
    mov = (Rotation.from_rotvec(ref) * Rotation.from_rotvec(rng.normal(0.0, 0.1, (N, 3)))).as_rotvec()
    _write_cs(tmp_path / "ref.cs", uids, ref)
    _write_cs(tmp_path / "mov.cs", uids, mov)
    api.run(tmp_path / "ref.cs", tmp_path / "mov.cs", output_dir=tmp_path / "run", k_neighbors=10, no_visualize=True)
    return tmp_path / "run"


@pytest.fixture
def window(qtbot, run_dir):
    win = MainWindow(NO_GL)
    qtbot.addWidget(win)
    win.resize(1500, 900)
    win.show()
    with qtbot.waitSignal(win.view_page.landscape_loaded, timeout=10000):
        win.open_run(run_dir)
    win.nav.setCurrentRow(PAGES.index("Select"))
    qtbot.wait(50)
    return win


def _click_on_particle(window, row: int) -> np.ndarray:
    page = window.view_page
    a, b = page.landscape.projections[0]
    coords = page.landscape.coordinates(page.representation, [row])[0]
    page._pick_2d(0, float(coords[a]), float(coords[b]))
    return coords


def _wait_for_count(qtbot, panel, radius=None):
    qtbot.waitUntil(
        lambda: panel.draft is not None and not panel._runner.busy
        and (radius is None or panel.draft.radius_deg == radius),
        timeout=10000,
    )


def test_select_page_shows_the_panel_and_enables_picking(window):
    page = window.view_page
    assert window.pages.currentWidget() is page
    assert page.selection_panel.isVisible() and page.picking
    assert "Space: <b>raw</b>" in page.selection_panel.space_label.text()
    window.nav.setCurrentRow(PAGES.index("View"))
    assert not page.selection_panel.isVisible() and not page.picking


def test_click_sets_centre_counts_exactly_and_highlights(window, qtbot):
    page, panel = window.view_page, window.view_page.selection_panel
    page.filter_combo.setCurrentIndex(page.filter_combo.findData("all"))
    row = int(np.argmax(page.landscape.sld_display))
    coords = _click_on_particle(window, row)
    _wait_for_count(qtbot, panel)
    assert [spin.value() for spin in panel.center_spins] == pytest.approx(list(coords), abs=1e-4)
    expected = page.landscape.evaluate_radius(coords, representation="euler", radius_deg=10.0)
    assert panel.draft.selected_count == expected.selected_count
    assert f"<b>{expected.selected_count:,}</b> of {N:,}" in panel.count_label.text()
    assert page.highlight_mask is panel.draft.selected_mask
    assert page.view2d.markers[0].data.size == 1
    assert page.view1d._overlays  # selection outline on the 1D histograms
    panel.radius_spin.setValue(25.0)
    _wait_for_count(qtbot, panel, radius=25.0)
    assert panel.draft.selected_count >= expected.selected_count
    assert panel.radius_slider.value() == 250


def test_confirm_writes_a_gui_selection_and_never_overwrites(window, qtbot, run_dir):
    page, panel = window.view_page, window.view_page.selection_panel
    _click_on_particle(window, 0)
    _wait_for_count(qtbot, panel)
    count = panel.draft.selected_count
    panel.name_edit.setText("region_01")
    assert panel.confirm_button.isEnabled()
    assert "--selection-id region_01" in panel.command_label.text()
    with qtbot.waitSignal(panel.selection_created, timeout=10000):
        panel.confirm_button.click()
    selection = read_selection_json(run_dir / "selections" / "region_01" / "selection.json")
    assert selection.selected_count == count
    assert selection.interaction_provenance["source"] == "cryorole_gui"
    assert "cryorole select" in window.commands_view.toPlainText()
    qtbot.waitUntil(lambda: panel.saved_list.count() == 1, timeout=10000)
    assert panel.saved_list.item(0).text() == f"region_01  ({count:,} particles)"

    panel.radius_spin.setValue(30.0)
    _wait_for_count(qtbot, panel, radius=30.0)
    panel.name_edit.setText("region_01")
    panel.confirm_button.click()
    qtbot.waitUntil(lambda: panel.findChild(object, "error_output_exists") is not None, timeout=10000)
    assert read_selection_json(run_dir / "selections" / "region_01" / "selection.json").selected_count == count


def test_threshold_form_creates_a_selection_like_the_cli(window, qtbot, run_dir):
    panel = window.view_page.selection_panel
    panel.mode_combo.setCurrentIndex(panel.mode_combo.findData("threshold"))
    assert not window.view_page.picking and window.view_page.highlight_mask is None
    panel.sld_min.setValue(2.0)
    panel.name_edit.setText("high_sld")
    assert "--mode threshold" in panel.command_label.text() and "--sld-min 2" in panel.command_label.text()
    with qtbot.waitSignal(panel.selection_created, timeout=10000):
        panel.confirm_button.click()
    sld_raw = window.view_page.landscape.session.arrays.sld_raw
    selection = read_selection_json(run_dir / "selections" / "high_sld" / "selection.json")
    assert selection.selected_count == int(np.count_nonzero(sld_raw >= 2.0))


def test_switching_coordinates_keeps_the_centre_and_the_count(window, qtbot):
    page, panel = window.view_page, window.view_page.selection_panel
    _click_on_particle(window, 5)
    _wait_for_count(qtbot, panel)
    count, center_rv = panel.draft.selected_count, panel.draft.center_rv
    page.representation_combo.setCurrentIndex(page.representation_combo.findData("rotvec"))
    qtbot.waitUntil(lambda: panel.draft is not None and panel.draft.representation == "rotvec"
                    and not panel._runner.busy, timeout=10000)
    assert panel.draft.selected_count == count
    assert panel.draft.center_rv == pytest.approx(center_rv, abs=1e-4)


def test_3d_pick_finds_the_point_under_the_cursor(qtbot):
    pytest.importorskip("pyqtgraph.opengl")
    from cryorole.gui.views.landscape3d import Landscape3DView

    view = Landscape3DView()
    qtbot.addWidget(view)
    view.resize(400, 400)
    coords = np.array([[0.0, 0.0, 0.0], [10.0, 10.0, 10.0], [-10.0, -10.0, -10.0]])
    view.set_points(coords, np.full((3, 4), 255, dtype=np.ubyte), ("a", "b", "c"),
                    ((-10.0, 10.0), (-10.0, 10.0), (-10.0, 10.0)))
    centre = (view.view.width() / 2.0, view.view.height() / 2.0)
    assert view.pick(*centre) == 0  # the camera looks at the box centre
    assert view.pick(1.0, 1.0) is None


def test_confirm_waits_for_the_count_of_the_current_radius(window, qtbot, run_dir):
    panel = window.view_page.selection_panel
    _click_on_particle(window, 0)
    _wait_for_count(qtbot, panel)
    panel.name_edit.setText("must_match")
    assert panel.confirm_button.isEnabled()
    panel.radius_spin.setValue(33.0)  # the shown count is for 10 deg until the new one arrives
    assert not panel.confirm_button.isEnabled()
    assert panel.count_label.text() == "Counting..."
    panel.confirm()  # ignored while out of date
    assert not (run_dir / "selections" / "must_match").exists()
    _wait_for_count(qtbot, panel, radius=33.0)
    qtbot.waitUntil(panel.confirm_button.isEnabled, timeout=10000)
    with qtbot.waitSignal(panel.selection_created, timeout=10000):
        panel.confirm_button.click()
    assert read_selection_json(run_dir / "selections" / "must_match" / "selection.json").radius == 33.0


def test_switching_coordinates_converts_the_typed_centre(window, qtbot):
    page, panel = window.view_page, window.view_page.selection_panel
    assert page.representation == "euler"
    for spin, value in zip(panel.center_spins, (30.0, 40.0, 50.0)):  # an Euler centre, in degrees
        spin.setValue(value)
    page.representation_combo.setCurrentIndex(page.representation_combo.findData("rotvec"))  # before the count
    _wait_for_count(qtbot, panel)
    expected_rv = page.landscape.rotvec_of((30.0, 40.0, 50.0), "euler")
    assert [spin.value() for spin in panel.center_spins] == pytest.approx(list(expected_rv), abs=1e-4)
    assert panel.draft.representation == "rotvec"
    assert list(panel.draft.center_rv) == pytest.approx(list(expected_rv), abs=1e-4)


def test_confirmed_selection_records_what_the_gui_showed(window, qtbot, run_dir):
    page, panel = window.view_page, window.view_page.selection_panel
    page.filter_combo.setCurrentIndex(page.filter_combo.findData("threshold"))
    page.filter_value.setValue(1.0)
    page.refresh(reset_view=False)
    shown = page.frame.filtered_count
    _click_on_particle(window, int(np.argmax(page.landscape.sld_display)))
    _wait_for_count(qtbot, panel)
    panel.name_edit.setText("shown_01")
    with qtbot.waitSignal(panel.selection_created, timeout=10000):
        panel.confirm_button.click()
    record = read_selection_json(run_dir / "selections" / "shown_01" / "selection.json").interaction_provenance
    assert record["displayed_point_count"] == shown < N
    assert record["display_filter"] == {"mode": "threshold", "field": "sld_display", "threshold": 1.0}
    assert record["full_candidate_count"] == N and record["selection_evaluation"] == "exact_full_parent_landscape"
