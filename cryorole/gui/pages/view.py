"""View page: 2D and 1D views of one landscape plus a separate 3D window, with display-only controls.

The 3D view (an OpenGL widget) lives in its own top-level window: on NVIDIA/X11
(measured 2026-10-04) an OpenGL widget anywhere in the main window made every
update of the 2D images stall for ~2.9 s, versus 6 ms with the OpenGL widget in a
separate window.

Defaults follow `cryorole visualize`: display filter sld_display >= 1, the
standard axis ranges (Euler +/-180 deg, rotation vector +/-3.14 rad) with equal
axis units, legacy colours. The 3D view draws at most "3D max points"
(deterministic sample); the 2D and 1D views use every particle that passes the
display filter. One axis range drives all three views.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shlex
import shutil

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSpinBox,
    QSplitter,
    QStackedWidget,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from cryorole import api
from cryorole.errors import CryoroleError
from cryorole.gui.pages.select import SelectionPanel
from cryorole.gui.process import WorkerProcess
from cryorole.gui.tasks import LatestWinsRunner
from cryorole.gui.views.colors import colormap_lut, values_to_rgba
from cryorole.gui.views.histogram1d import HistogramPanels
from cryorole.gui.views.legend import ColorLegend
from cryorole.gui.views.projection2d import ProjectionPanels
from cryorole.gui.widgets.errors import show_error
from cryorole.gui.widgets.save_figures_dialog import SaveFiguresDialog
from cryorole.logs import notify_user

REPRESENTATIONS = (("euler", "Euler angles (deg)"), ("rotvec", "Rotation vector (rad)"))
FILTERS = (("threshold", "SLD at least"), ("top_percent", "Top SLD %"), ("all", "All particles"))
DEFAULT_SLD_THRESHOLD = 1.0  # `cryorole visualize` default: sld_display >= 1
DEFAULT_TOP_PERCENT = 40.0
AXIS_MODES = (
    ("standard", "Standard (as figures)"),
    ("core", "Fit core (central 99%)"),
    ("all", "Fit all particles"),
    ("custom", "Custom"),
)
CORE_COVERAGE = 0.99  # "Fit core": view range only; outliers stay in the data and in the plots
SETTINGS_DELAY_MS = 150
PICK_RADIUS_PX = 15  # a click farther than this from every displayed particle picks nothing
DIM_RGBA = np.array([215, 215, 215, 255], dtype=np.ubyte)  # unselected particles while a selection is shown
Range = tuple[float, float]


class ViewPage(QWidget):
    """View one API landscape; explicit Save and Confirm call its standard writers."""

    landscape_loaded = Signal(object)
    center_picked = Signal(int)
    command_issued = Signal(str)  # landscape row of the displayed particle nearest to a click
    figures_saved = Signal(object)
    overlay_loaded = Signal(object)

    def __init__(self, *, gl_available: bool, gl_reason: str | None = None, parent=None) -> None:
        super().__init__(parent)
        self.gl_available = gl_available
        self.landscape = None
        self.frame = None
        self.axis_ranges: list[Range] = []
        self.highlight_mask: np.ndarray | None = None  # exact full-landscape mask of the shown selection
        self.highlight_center_rv: tuple[float, float, float] | None = None
        self.dim_unselected = True
        self.picking = False
        self._load_runner = LatestWinsRunner(self.set_landscape, self._load_failed)
        self.destroyed.connect(self._load_runner.discard)
        self._load_generation = 0
        self._overlay_runner = LatestWinsRunner(self._overlay_ready, self._overlay_failed)
        self.destroyed.connect(self._overlay_runner.discard)
        self._saved_overlay = None
        self._figure_task = None
        self._run_dir: Path | None = None
        self._canonicalize_run_dir: Path | None = None
        self._canonicalize_target: Path | None = None  # the frame folder this canonicalize creates
        self._canonicalize_pending = False  # started, outcome not reported yet

        # What is shown (display only).
        self.space_combo = QComboBox()
        self.space_combo.setObjectName("space")
        self.representation_combo = QComboBox()
        for key, label in REPRESENTATIONS:
            self.representation_combo.addItem(label, key)
        self.filter_combo = QComboBox()
        for key, label in FILTERS:
            self.filter_combo.addItem(label, key)
        self.filter_value = QDoubleSpinBox()
        self.filter_value.setDecimals(3)
        self.vmax_spin = QDoubleSpinBox()
        self.vmax_spin.setObjectName("color_vmax")
        self.vmax_spin.setRange(0.0, 1e9)
        self.vmax_spin.setDecimals(2)
        self.vmax_spin.setSpecialValueText("auto")  # 0 shows "auto": vmax = min(max SLD, 100)
        self.vmax_spin.setValue(0.0)
        self.vmax_spin.setToolTip(
            "Top of the SLD colour scale (display only, like `cryorole visualize --vmax`). "
            "auto = min(max SLD, 100), as in the standard figures."
        )
        self.legend = ColorLegend()
        self.save_figures_button = QPushButton("Save figures...")
        self.save_figures_button.setObjectName("save_figures")
        self.save_figures_button.setEnabled(False)
        self.save_figures_button.clicked.connect(self.ask_save_figures)
        self.canonicalize_button = QPushButton("Canonicalize...")
        self.canonicalize_button.setObjectName("canonicalize")
        self.canonicalize_button.setToolTip("Fit a canonical frame (cryorole canonicalize) and open it.")
        self._canonicalizer = WorkerProcess(self)
        self._canonicalizer.succeeded.connect(self._canonicalized)
        self._canonicalizer.failed.connect(self._canonicalize_failed)
        self._pending_frame: str | None = None

        # How it is framed (display only).
        self.axis_combo = QComboBox()
        self.axis_combo.setObjectName("axis_range")
        for key, label in AXIS_MODES:
            self.axis_combo.addItem(label, key)
        self.axis_combo.setToolTip("Axis ranges of the 2D, 1D and 3D views. View only: no particle is filtered out.")
        self.custom_labels: list[QLabel] = []
        self.custom_spins: list[tuple[QDoubleSpinBox, QDoubleSpinBox]] = []
        custom_row = QHBoxLayout()
        custom_row.setContentsMargins(0, 0, 0, 0)
        for _ in range(3):
            label = QLabel()
            low, high = QDoubleSpinBox(), QDoubleSpinBox()
            for spin in (low, high):
                spin.setRange(-1e6, 1e6)
                spin.setDecimals(2)
                spin.setKeyboardTracking(False)
            custom_row.addWidget(label)
            custom_row.addWidget(low)
            custom_row.addWidget(QLabel("to"))
            custom_row.addWidget(high)
            self.custom_labels.append(label)
            self.custom_spins.append((low, high))
        self.custom_box = QWidget()
        self.custom_box.setLayout(custom_row)
        self.custom_box.setVisible(False)
        self.custom_box.setToolTip("Like `cryorole visualize --axis-limit AXIS:LOWER:UPPER`: view only, no filtering.")
        self.equal_units = QCheckBox("Equal axis units")
        self.equal_units.setChecked(True)
        self.equal_units.setToolTip("Same scale on both axes of each 2D panel, as in the standard figures.")
        self.reset_button = QPushButton("Reset view")
        self.reset_button.setToolTip("Undo mouse zoom and pan: return to the chosen axis range.")
        self.cap_spin = QSpinBox()
        self.cap_spin.setRange(1_000, 5_000_000)
        self.cap_spin.setSingleStep(50_000)
        self.cap_spin.setGroupSeparatorShown(True)
        self.cap_spin.setValue(100_000)
        self.open3d_button = QPushButton("3D view...")
        self.open3d_button.setObjectName("open_3d")
        if gl_available:
            self.open3d_button.setToolTip("Open the 3D view in its own window; it follows these settings.")
        else:
            self.open3d_button.setEnabled(False)
            self.open3d_button.setToolTip(
                f"3D view is off: {gl_reason or 'OpenGL is not available'}. The 2D and 1D views do not need OpenGL."
            )
        self.cap_spin.setToolTip(
            "Most points drawn in the 3D view (deterministic sample, seed 0). "
            "The 2D and 1D views always use every particle that passes the display filter."
        )

        shown = QHBoxLayout()
        for widget in (
            QLabel("Space"), self.space_combo, self.canonicalize_button,
            QLabel("Coordinates"), self.representation_combo,
            QLabel("Display filter"), self.filter_combo, self.filter_value,
            QLabel("SLD colour max"), self.vmax_spin,
        ):
            shown.addWidget(widget)
        shown.addStretch(1)
        shown.addWidget(self.legend)
        framed = QHBoxLayout()
        for widget in (QLabel("Axis range"), self.axis_combo, self.custom_box, self.equal_units, self.reset_button):
            framed.addWidget(widget)
        framed.addStretch(1)
        framed.addWidget(self.save_figures_button)
        framed.addWidget(self.open3d_button)
        framed.addWidget(QLabel("3D max points"))
        framed.addWidget(self.cap_spin)

        # Views (the 3D window is created on first use).
        self.view3d = None
        self.tabs = QTabWidget()
        self.view2d = ProjectionPanels()
        self.view1d = HistogramPanels()
        self.tabs.addTab(self.view2d, "2D")
        self.tabs.addTab(self.view1d, "1D")

        self.status = QLabel()
        self.status.setObjectName("view_status")
        self.status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.addLayout(shown)
        content_layout.addLayout(framed)
        overlay_row = QHBoxLayout()
        self.overlay_combo = QComboBox()
        self.overlay_combo.setObjectName("saved_selection_overlay")
        self.overlay_combo.addItem("No saved selection", None)
        self.overlay_combo.setToolTip("Review an existing Selection. This does not evaluate or confirm a new selection.")
        self.overlay_combo.currentIndexChanged.connect(self._overlay_chosen)
        self.overlay_status = QLabel("Read-only overlay")
        self.overlay_status.setObjectName("overlay_status")
        self.overlay_refresh = QPushButton("Refresh list")
        self.overlay_refresh.clicked.connect(self.refresh_overlays)
        overlay_row.addWidget(QLabel("Show saved selection"))
        overlay_row.addWidget(self.overlay_combo)
        overlay_row.addWidget(self.overlay_refresh)
        overlay_row.addWidget(self.overlay_status, 1)
        content_layout.addLayout(overlay_row)
        self.figure_status = QLabel()
        self.figure_status.setObjectName("figure_export_status")
        self.figure_status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.figure_status.setWordWrap(True)
        content_layout.addWidget(self.figure_status)
        content_layout.addWidget(self.tabs, 1)
        content_layout.addWidget(self.status)
        self.empty = QLabel("Open a run bundle (File > Open run...) to view its landscape.")
        self.empty.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.stack = QStackedWidget()
        self.stack.addWidget(self.empty)
        self.stack.addWidget(content)
        self.selection_panel = SelectionPanel(self)
        self.selection_panel.setVisible(False)
        self.selection_panel.selection_created.connect(self.refresh_overlays)
        splitter = QSplitter()
        splitter.addWidget(self.stack)
        splitter.addWidget(self.selection_panel)
        splitter.setStretchFactor(0, 1)
        splitter.setCollapsible(1, False)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter)

        self._set_filter_mode("threshold")

        self._settings_timer = QTimer(self)
        self._settings_timer.setSingleShot(True)
        self._settings_timer.setInterval(SETTINGS_DELAY_MS)
        self._settings_timer.timeout.connect(lambda: self.refresh(reset_view=False))
        self._range_timer = QTimer(self)
        self._range_timer.setSingleShot(True)
        self._range_timer.setInterval(SETTINGS_DELAY_MS)
        self._range_timer.timeout.connect(lambda: self.refresh(reset_view=True))
        self.representation_combo.currentIndexChanged.connect(self._representation_changed)
        self.filter_combo.currentIndexChanged.connect(lambda _: self._filter_mode_changed())
        self.filter_value.valueChanged.connect(lambda _: self._settings_timer.start())
        self.cap_spin.valueChanged.connect(lambda _: self._settings_timer.start())
        self.vmax_spin.valueChanged.connect(lambda _: self._settings_timer.start())
        self.axis_combo.currentIndexChanged.connect(lambda _: self._axis_mode_changed())
        for low, high in self.custom_spins:
            low.valueChanged.connect(lambda _: self._range_timer.start())
            high.valueChanged.connect(lambda _: self._range_timer.start())
        self.equal_units.toggled.connect(lambda _: self.refresh(reset_view=True))
        self.reset_button.clicked.connect(lambda: self.refresh(reset_view=True))
        self.space_combo.activated.connect(self._space_chosen)
        self.view2d.point_clicked.connect(self._pick_2d)
        self.open3d_button.clicked.connect(self.open_3d)
        self.canonicalize_button.clicked.connect(self.ask_canonicalize)

    # --- explicit display exports / saved-selection review ----------------------

    @property
    def saving_figures(self) -> bool:
        return self._figure_task is not None

    def _figure_snapshot(self):
        reset = self._range_timer.isActive()
        self._settings_timer.stop()
        self._range_timer.stop()
        self.refresh(reset_view=reset)
        if self.landscape is None or self.frame is None:
            raise CryoroleError("input_invalid", "Show a non-empty landscape before saving figures.")
        return self.landscape, self.frame, self.representation, tuple(self.axis_ranges)

    def ask_save_figures(self) -> None:
        if self.saving_figures:
            return
        try:
            snapshot = self._figure_snapshot()
        except CryoroleError as exc:
            show_error(self, exc)
            return
        landscape, _, representation, _ = snapshot
        description = f"{landscape.run_dir.name} · {landscape.space} · {representation}"
        if landscape.canonical_id:
            description += f" · {landscape.canonical_id}"
        dialog = SaveFiguresDialog(description, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self._start_figure_export(snapshot, dialog.visual_id, dialog.views)

    def save_figures(self, visual_id: str, *, views: str = "2d,1d") -> None:
        """Explicit action; snapshot controls before the background writer starts."""

        if not self.saving_figures:
            self._start_figure_export(self._figure_snapshot(), visual_id, views)

    def _start_figure_export(self, snapshot, visual_id: str, views: str) -> None:
        """Save figures in a child process (`cryorole.gui.worker`), never in a GUI-process thread.

        Matplotlib running in a thread beside Qt painting can crash or hang the GUI
        (PySide6 on Python 3.10/3.11), so the GUI only captures the display state here.
        """

        landscape, frame, representation, ranges = snapshot
        try:
            request = landscape.figure_request(visual_id, frame=frame, representation=representation,
                                               axis_ranges=ranges, views=views)
        except CryoroleError as exc:
            self._figures_failed(exc)
            return
        self.save_figures_button.setEnabled(False)
        self.figure_status.setText(f"Saving {visual_id} from {landscape.run_dir.name} ({landscape.space})...")
        worker = WorkerProcess(self)  # one per figure set: a finished worker may still be exiting
        worker.succeeded.connect(lambda message: self._figures_saved(
            FigureResult(Path(message["output_dir"]), message["report"])))
        worker.failed.connect(lambda message: self._figures_failed(CryoroleError(
            message.get("code") or "internal", message.get("message") or "Saving figures failed.",
            message.get("remedy"))))
        worker.stopped.connect(worker.deleteLater)
        self._figure_task = worker
        worker.start({"command": "save_figures", **request})

    def _figures_saved(self, result) -> None:
        self._figure_task = None
        self.save_figures_button.setEnabled(self.frame is not None)
        self.figure_status.setText(f"Saved figures: {result.output_dir}")
        self.command_issued.emit(result.report["gui_display_snapshot"]["equivalent_cli_command"])
        self.figures_saved.emit(result)

    def _figures_failed(self, exc: BaseException) -> None:
        self._figure_task = None
        self.save_figures_button.setEnabled(self.frame is not None)
        self.figure_status.setText("Figures were not saved. Choose a new name or correct the error and retry.")
        if isinstance(exc, FileExistsError) or getattr(exc, "code", None) == "output_exists":
            exc = CryoroleError("output_exists", "A figure set already uses this name.",
                                "Choose a new figure name; existing files were kept.")
        show_error(self, exc, title="cryoROLE: cannot save figures")

    def refresh_overlays(self, *_args, preferred: str | None = None) -> None:
        current = self.overlay_combo.currentData() if preferred is None else preferred
        self.overlay_combo.blockSignals(True)
        self.overlay_combo.clear()
        self.overlay_combo.addItem("No saved selection", None)
        if self.landscape is not None:
            for name, count in self.landscape.saved_selections():
                label = name if count is None else f"{name} ({count:,})"
                self.overlay_combo.addItem(label, name)
        self.overlay_combo.setCurrentIndex(max(0, self.overlay_combo.findData(current)))
        self.overlay_combo.blockSignals(False)
        self._overlay_chosen()

    def _overlay_chosen(self, *_args) -> None:
        self._overlay_runner.discard()
        self._saved_overlay = None
        self._update_overlay_status()
        if self.landscape is not None and self.frame is not None:
            self._draw(reset_view=False)
        name, landscape = self.overlay_combo.currentData(), self.landscape
        if name is not None and landscape is not None:
            self.overlay_status.setText(f"Loading {name} (read only)...")
            self._overlay_runner.submit(lambda: (landscape, landscape.saved_selection_overlay(name)))

    def _overlay_ready(self, result) -> None:
        landscape, overlay = result
        if self.landscape is not landscape or self.overlay_combo.currentData() != overlay.selection_id:
            return
        self._saved_overlay = overlay
        self._update_overlay_status()
        if self.frame is not None:
            self._draw(reset_view=False)
        self.overlay_loaded.emit(overlay)

    def _overlay_failed(self, exc: BaseException) -> None:
        self._saved_overlay = None
        self.overlay_status.setText("Saved selection could not be shown; no selection changed.")
        show_error(self, exc, title="cryoROLE: cannot show saved selection")

    def _display_highlight(self):
        if self.highlight_mask is not None:
            return self.highlight_mask
        return None if self._saved_overlay is None else self._saved_overlay.selected_mask

    def _update_overlay_status(self) -> None:
        overlay = self._saved_overlay
        if overlay is None:
            self.overlay_status.setText("Read-only overlay")
            return
        visible = 0 if self.frame is None else int(overlay.selected_mask[self.frame.filtered_indices].sum())
        text = f"{overlay.selection_id}: {overlay.selected_count:,} of {overlay.total_count:,}; {visible:,} visible · read only"
        if self.highlight_mask is not None:
            text += " · draft highlight shown"
        self.overlay_status.setText(text)

    # --- canonicalize ------------------------------------------------------------

    def ask_canonicalize(self) -> None:
        if self.landscape is None or self._canonicalizer.running:
            return
        from cryorole.gui.widgets.canonicalize_dialog import CanonicalizeDialog

        dialog = CanonicalizeDialog(self, existing=self.landscape.canonical_ids)
        if dialog.exec():
            self.canonicalize(**dialog.values())

    def canonicalize(self, *, canonical_id: str, fit_top_fraction: float, visualize: bool) -> None:
        """Run `cryorole canonicalize` in the worker; open the new frame when it is done."""

        run_dir = self.landscape.run_dir
        self._canonicalize_generation = self._load_generation
        self._canonicalize_run_dir = run_dir
        self._pending_frame = canonical_id
        target = Path(run_dir) / "canonical" / canonical_id
        self._canonicalize_target = None if target.exists() else target
        self._canonicalize_pending = True
        self.canonicalize_button.setEnabled(False)
        self.status.setText(f"Canonicalizing ({canonical_id}, fit to the top {fit_top_fraction:.0%} by SLD)...")
        request = {"command": "canonicalize", "run_dir": str(run_dir), "canonical_id": canonical_id,
                   "fit_top_fraction": fit_top_fraction}
        if not visualize:
            request["no_visualize"] = True
        self._canonicalize_command = shlex.join(
            ["cryorole", "canonicalize", "--run-dir", str(run_dir), "--canonical-id", canonical_id,
             "--fit-top", f"{fit_top_fraction:g}"] + ([] if visualize else ["--no-visualize"])
        )
        self._canonicalizer.start(request)

    @property
    def canonicalizing(self) -> bool:  # the worker reports before it exits: a reported frame is complete
        return self._canonicalize_pending and self._canonicalizer.running

    def stop_canonicalize(self) -> None:
        """Stop the canonicalize worker and remove the frame folder it started (it would be incomplete)."""

        if not self.canonicalizing:
            return
        self._canonicalizer.blockSignals(True)  # quitting: no error dialog for the stopped worker
        self._canonicalizer.kill()
        self._canonicalizer.process.waitForFinished(5000)
        if self._canonicalize_target is not None:
            shutil.rmtree(self._canonicalize_target, ignore_errors=True)

    def _canonicalized(self, result: dict) -> None:
        self._canonicalize_pending = False
        self.canonicalize_button.setEnabled(True)
        notify_user(f"Canonical frame {self._pending_frame} written: {result.get('output_dir')}")
        self.command_issued.emit(self._canonicalize_command)
        run_dir = self._canonicalize_run_dir
        same_request = self._load_generation == self._canonicalize_generation
        if same_request and self.landscape is not None and self.landscape.run_dir == run_dir:
            self.open_run(run_dir, space="canonical", canonical_id=self._pending_frame)
        else:
            self.status.setText(f"Canonical frame {self._pending_frame} of {Path(run_dir).name} is ready "
                                "(File > Open run... and choose it under Space).")

    def _canonicalize_failed(self, error: dict) -> None:
        self._canonicalize_pending = False
        self.canonicalize_button.setEnabled(True)
        self.status.setText(self.status_text())
        show_error(self, CryoroleError(error.get("code", "internal"), error.get("message", ""), error.get("remedy")),
                   title="cryoROLE: canonicalize failed")

    # --- 3D window ---------------------------------------------------------------

    def open_3d(self) -> None:
        """Show the 3D window (created on first use) with the current display."""

        if not self.gl_available:
            return
        if self.view3d is None:
            from cryorole.gui.views.landscape3d import Landscape3DView

            self.view3d = Landscape3DView()
            self.view3d.resize(1000, 850)
            self.view3d.point_clicked.connect(self._pick_3d)
        self.view3d.setWindowTitle(self._3d_title())
        self.view3d.show()
        self.view3d.raise_()
        self.view3d.activateWindow()
        if self.frame is not None:
            self._draw_3d(reset_view=True)

    def close_3d(self) -> None:
        if self.view3d is not None:
            self.view3d.close()

    def _3d_window_shown(self) -> bool:
        return self.view3d is not None and self.view3d.isVisible()

    def _3d_title(self) -> str:
        landscape = self.landscape
        if landscape is None:
            return "cryoROLE 3D"
        where = "raw" if landscape.space == "raw" else f"canonical {landscape.canonical_id}"
        return f"cryoROLE 3D - {landscape.run_dir.name} ({where})"

    def _draw_3d(self, *, reset_view: bool) -> None:
        landscape, frame, representation = self.landscape, self.frame, self.representation
        lut = colormap_lut(landscape.style["color_map"])
        rows = frame.sample_indices
        rgba = values_to_rgba(landscape.sld_display[rows], frame.vmin, frame.vmax, lut)
        mask = self._display_highlight()
        dimmed = ~mask[rows] if mask is not None and self.dim_unselected else None
        self.view3d.set_points(
            landscape.coordinates(representation, rows), rgba, landscape.axis_labels(representation),
            self.axis_ranges if reset_view else None, dimmed=dimmed,
        )
        center = None
        if self.highlight_center_rv is not None:
            center = np.asarray(landscape.represent(self.highlight_center_rv, representation))
        self.view3d.set_marker(center)
        self.view3d.setWindowTitle(self._3d_title())

    # --- selection support (display only; counts come from the exact evaluator) ---

    def set_selecting(self, selecting: bool) -> None:
        """Show or hide the Select panel next to the views."""

        self.selection_panel.setVisible(selecting)
        self.selection_panel.on_shown(selecting)

    def set_highlight(self, mask: np.ndarray | None, center_rv=None) -> None:
        """Show an exact selection mask: selected particles in SLD colours, the rest dimmed."""

        self.highlight_mask = mask
        self.highlight_center_rv = center_rv
        self._update_overlay_status()
        if self.landscape is not None and self.frame is not None:
            self._draw(reset_view=False)

    def set_dim_unselected(self, dim: bool) -> None:
        self.dim_unselected = dim
        if self.highlight_mask is not None:
            self.set_highlight(self.highlight_mask, self.highlight_center_rv)

    def highlighted_display_count(self) -> int:
        """How many particles of the shown selection pass the display filter (and are drawn)."""

        if self.highlight_mask is None or self.frame is None:
            return 0
        return int(np.count_nonzero(self.highlight_mask[self.frame.filtered_indices]))

    def _pick_2d(self, panel: int, x: float, y: float) -> None:
        if not self.picking or self.frame is None:
            return
        a, b = self.landscape.projections[panel]
        rows = self.frame.filtered_indices
        if rows.size == 0:
            return
        coords = self.landscape.coordinates(self.representation, rows)
        px, py = self.view2d.pixel_size(panel)
        distance = np.hypot((coords[:, a] - x) / px, (coords[:, b] - y) / py)  # in screen pixels
        best = int(np.argmin(distance))
        if distance[best] > PICK_RADIUS_PX:
            self.status.setText("No displayed particle near the click; click on the landscape.")
            return
        self.center_picked.emit(int(rows[best]))

    def _pick_3d(self, index: int) -> None:
        if self.picking and self.frame is not None and self.view3d is not None:
            self.center_picked.emit(int(self.frame.sample_indices[index]))

    # --- loading ---------------------------------------------------------------

    def open_run(self, run_dir: str | Path, *, space: str = "raw", canonical_id: str | None = None) -> None:
        """Load a landscape in the background; errors are shown, nothing is written."""

        self.status.setText(f"Loading {run_dir} ({space}{' ' + canonical_id if canonical_id else ''}) ...")
        self.setEnabled(False)
        self._load_generation += 1

        def load():
            return api.open_landscape(run_dir, space=space, canonical_id=canonical_id)

        self._load_runner.submit(load)

    def _load_failed(self, exc: BaseException) -> None:
        self.setEnabled(True)
        self.status.setText(self.status_text())
        show_error(self, exc, title="cryoROLE: cannot open run")

    def set_landscape(self, landscape) -> None:
        self.setEnabled(True)
        previous = self.overlay_combo.currentData() if self.landscape is not None and self.landscape.run_id == landscape.run_id else None
        self._overlay_runner.discard()
        self._saved_overlay = None
        self.landscape = landscape
        self.highlight_mask = None  # a selection belongs to one landscape (space)
        self.highlight_center_rv = None
        self._run_dir = landscape.run_dir
        self.space_combo.blockSignals(True)
        self.space_combo.clear()
        self.space_combo.addItem("Raw", ("raw", None))
        for frame_id in landscape.canonical_ids:
            self.space_combo.addItem(f"Canonical: {frame_id}", ("canonical", frame_id))
        # findData() cannot match Python tuples (PySide6 compares the wrapped objects by identity).
        key = (landscape.space, landscape.canonical_id)
        current = next((i for i in range(self.space_combo.count()) if self.space_combo.itemData(i) == key), 0)
        self.space_combo.setCurrentIndex(current)
        self.space_combo.blockSignals(False)
        self._fill_custom_ranges()
        self.refresh(reset_view=True)
        self.refresh_overlays(preferred=previous)
        self.landscape_loaded.emit(landscape)

    def _space_chosen(self, _index: int) -> None:
        space, frame_id = self.space_combo.currentData()
        if self.landscape is not None and (space, frame_id) != (self.landscape.space, self.landscape.canonical_id):
            self.open_run(self._run_dir, space=space, canonical_id=frame_id)

    # --- display settings ------------------------------------------------------

    def _set_filter_mode(self, mode: str) -> None:
        self.filter_value.blockSignals(True)
        self.filter_combo.blockSignals(True)
        self.filter_combo.setCurrentIndex(self.filter_combo.findData(mode))
        if mode == "threshold":
            self.filter_value.setSuffix("")
            self.filter_value.setRange(0.0, 1e9)
            self.filter_value.setSingleStep(0.5)
            self.filter_value.setValue(DEFAULT_SLD_THRESHOLD)
        elif mode == "top_percent":
            self.filter_value.setSuffix(" %")
            self.filter_value.setRange(0.1, 100.0)
            self.filter_value.setSingleStep(5.0)
            self.filter_value.setValue(DEFAULT_TOP_PERCENT)
        self.filter_value.setEnabled(mode != "all")
        self.filter_combo.blockSignals(False)
        self.filter_value.blockSignals(False)

    def _filter_mode_changed(self) -> None:
        self._set_filter_mode(self.filter_combo.currentData())
        self.refresh(reset_view=self.axis_mode in ("core", "all"))

    def display_options(self) -> dict[str, float | int | None]:
        mode = self.filter_combo.currentData()
        return {
            "threshold": self.filter_value.value() if mode == "threshold" else None,
            "top_fraction": self.filter_value.value() / 100.0 if mode == "top_percent" else None,
            "max_points_3d": self.cap_spin.value(),
            "color_vmax": self.vmax_spin.value() or None,
        }

    def filter_description(self) -> str:
        mode = self.filter_combo.currentData()
        if mode == "threshold":
            return f"with SLD >= {self.filter_value.value():g}"
        if mode == "top_percent":
            return f"in the top {self.filter_value.value():g}% by SLD"
        return "(no display filter)"

    @property
    def representation(self) -> str:
        return self.representation_combo.currentData()

    @property
    def axis_mode(self) -> str:
        return self.axis_combo.currentData()

    def _representation_changed(self, _index: int) -> None:
        self._fill_custom_ranges()
        self.refresh(reset_view=True)

    def _axis_mode_changed(self) -> None:
        self.custom_box.setVisible(self.axis_mode == "custom")
        self.refresh(reset_view=True)

    def _fill_custom_ranges(self) -> None:
        """Start custom ranges from the standard ones of the current coordinates."""

        if self.landscape is None:
            return
        names = self.landscape.axis_names(self.representation)
        standard = self.landscape.axis_limits(self.representation)
        for label, name, (low, high), (lo, hi) in zip(self.custom_labels, names, self.custom_spins, standard):
            label.setText(name)
            for spin, value in ((low, lo), (high, hi)):
                spin.blockSignals(True)
                spin.setValue(value)
                spin.blockSignals(False)

    def set_custom_ranges(self, ranges: list[Range]) -> None:
        """Switch to Custom with these (low, high) ranges, one per axis."""

        for (low, high), (lo, hi) in zip(self.custom_spins, ranges):
            for spin, value in ((low, lo), (high, hi)):
                spin.blockSignals(True)
                spin.setValue(value)
                spin.blockSignals(False)
        if self.axis_mode == "custom":
            self.refresh(reset_view=True)
        else:
            self.axis_combo.setCurrentIndex(self.axis_combo.findData("custom"))

    def _resolve_axis_ranges(self, filtered: np.ndarray) -> list[Range]:
        mode = self.axis_mode
        if mode == "standard":
            return [tuple(limits) for limits in self.landscape.axis_limits(self.representation)]
        if mode == "custom":
            ranges = [(low.value(), high.value()) for low, high in self.custom_spins]
            bad = [label.text() for label, (lo, hi) in zip(self.custom_labels, ranges) if not hi > lo]
            if bad:
                raise CryoroleError(
                    "input_invalid",
                    f"Custom axis range needs lower < upper for: {', '.join(bad)}.",
                    "Fix the range, or choose another axis range.",
                )
            return ranges
        return [_fitted_range(filtered[:, axis], mode) for axis in range(3)]

    # --- drawing ---------------------------------------------------------------

    def refresh(self, *, reset_view: bool) -> None:
        landscape = self.landscape
        if landscape is None:
            return
        try:
            frame = landscape.apply(**self.display_options())
            representation = self.representation
            filtered = landscape.coordinates(representation, frame.filtered_indices)
            axis_ranges = self._resolve_axis_ranges(filtered) if reset_view or not self.axis_ranges else self.axis_ranges
        except CryoroleError as exc:
            # Nothing to show (e.g. the display filter keeps no particle of this landscape): clear the
            # views, so no picture, count or click target of the previous landscape or filter remains.
            message = f"{exc.message} {exc.remedy or ''}".strip()
            self.frame, self._base_cache = None, None
            self.save_figures_button.setEnabled(False)
            self._update_overlay_status()
            self.empty.setText(message)
            self.stack.setCurrentIndex(0)
            if self._3d_window_shown():
                self.view3d.set_points(np.zeros((0, 3)), np.zeros((0, 4), dtype=np.ubyte),
                                       landscape.axis_labels(self.representation))
            self.status.setText(message)
            return
        self.frame = frame
        self.save_figures_button.setEnabled(not self.saving_figures)
        self._update_overlay_status()
        self.axis_ranges = axis_ranges
        self.stack.setCurrentIndex(1)
        self._draw(reset_view=reset_view)

    def _base(self) -> dict:
        """Per (display settings, coordinates) data that a highlight change does not alter."""

        frame, representation = self.frame, self.representation
        key = (id(self.landscape), id(frame), representation)
        if getattr(self, "_base_cache", None) is None or self._base_cache["key"] != key:
            landscape = self.landscape
            coords = landscape.coordinates(representation, frame.filtered_indices)
            self._base_cache = {
                "key": key,
                "coords": coords,
                "sld": landscape.sld_display[frame.filtered_indices],
                "histograms": [landscape.histogram_of(coords[:, axis], frame) for axis in range(3)],
                "sld_histogram": landscape.sld_histogram(frame),
            }
        return self._base_cache

    def _draw(self, *, reset_view: bool) -> None:
        landscape, frame, axis_ranges = self.landscape, self.frame, self.axis_ranges
        representation = self.representation
        lut = colormap_lut(landscape.style["color_map"])
        self.legend.set_scale(lut, frame.vmin, frame.vmax)
        labels = landscape.axis_labels(representation)
        base = self._base()
        coords, sld = base["coords"], base["sld"]
        mask = self._display_highlight()
        center = None
        if self.highlight_center_rv is not None:
            center = np.asarray(landscape.represent(self.highlight_center_rv, representation))

        if self._3d_window_shown():
            self._draw_3d(reset_view=reset_view)

        chosen = mask[frame.filtered_indices] if mask is not None else None
        if chosen is not None:
            layers = [(coords[~chosen], sld[~chosen]), (coords[chosen], sld[chosen])]
        dim = self.dim_unselected

        def render(index, x_range, y_range, width, height):
            a, b = landscape.projections[index]

            def image(points, values):
                return landscape.raster(points[:, a], points[:, b], values, x_range=x_range, y_range=y_range,
                                        width=width, height=height)

            if chosen is None:
                return values_to_rgba(image(coords, sld), frame.vmin, frame.vmax, lut)
            (other_xyz, other_sld), (chosen_xyz, chosen_sld) = layers
            others = image(other_xyz, other_sld)
            if dim:
                rgba = np.zeros(others.shape + (4,), dtype=np.ubyte)
                rgba[np.isfinite(others)] = DIM_RGBA
            else:
                rgba = values_to_rgba(others, frame.vmin, frame.vmax, lut)
            top = image(chosen_xyz, chosen_sld)
            drawn = np.isfinite(top)
            rgba[drawn] = values_to_rgba(top[drawn], frame.vmin, frame.vmax, lut)
            return rgba

        self.view2d.set_source(
            render,
            [(labels[a], labels[b]) for a, b in landscape.projections],
            [(axis_ranges[a], axis_ranges[b]) for a, b in landscape.projections],
            equal_units=self.equal_units.isChecked(),
            reset_view=reset_view,
        )
        self.view2d.set_markers(
            None if center is None else [(center[a], center[b]) for a, b in landscape.projections]
        )
        histograms = base["histograms"]
        self.view1d.set_data(histograms, labels, base["sld_histogram"], axis_ranges)
        if chosen is not None:
            self.view1d.set_overlay([
                landscape.histogram_of(coords[chosen, axis], frame, edges=edges)
                for axis, (_counts, edges) in enumerate(histograms)
            ])
        self.status.setText(self.status_text())

    def status_text(self) -> str:
        frame, landscape = self.frame, self.landscape
        if frame is None or landscape is None:
            return ""
        where = "raw" if landscape.space == "raw" else f"canonical {landscape.canonical_id}"
        parts = [f"{landscape.run_dir.name} ({where})"]
        if self.gl_available:
            sampled = " (display sample, seed 0)" if frame.downsampled else ""
            parts.append(f"3D: showing {frame.sample_count:,} of {frame.filtered_count:,}{sampled}")
        parts.append(f"2D and 1D: all {frame.filtered_count:,} particles {self.filter_description()}")
        parts.append(f"full landscape {frame.full_count:,}")
        return " · ".join(parts) + " · display only"


def _fitted_range(values: np.ndarray, mode: str, margin: float = 0.05) -> Range:
    if mode == "core" and values.size > 100:
        tail = (1.0 - CORE_COVERAGE) / 2.0 * 100.0
        low, high = (float(v) for v in np.percentile(values, [tail, 100.0 - tail]))
    else:
        low, high = float(np.min(values)), float(np.max(values))
    span = high - low if high > low else 1.0
    return low - margin * span, high + margin * span


@dataclass(frozen=True)
class FigureResult:
    """A saved figure set, as the figure worker reported it (like ``VisualizationResult``)."""

    output_dir: Path
    report: dict
