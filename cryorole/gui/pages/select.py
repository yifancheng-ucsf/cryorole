"""Select panel, shown next to the views on the Select page.

Radius mode is interactive: click a displayed particle (2D or 3D) or type a
centre, move the radius, and the exact count over the full landscape updates
(the CLI's SO(3) evaluator, in the view's space). Nothing is written until
Confirm, which creates a standard Selection and never overwrites. The other
modes (SLD threshold, coordinate range, random, metadata) call `api.select`.
Display filters only change what is drawn; no selection uses them.
"""

from __future__ import annotations

import shlex

from PySide6.QtCore import QRegularExpression, Qt, QTimer, Signal
from PySide6.QtGui import QRegularExpressionValidator
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QPushButton,
    QSlider,
    QSpinBox,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from cryorole import api
from cryorole.errors import CryoroleError
from cryorole.gui.tasks import LatestWinsRunner, run_in_background
from cryorole.gui.widgets.errors import show_error
from cryorole.logs import notify_user

MODES = (
    ("radius", "Radius around a centre (interactive)"),
    ("threshold", "SLD threshold"),
    ("range", "Coordinate range"),
    ("random", "Random fraction"),
    ("metadata", "Source metadata values"),
)
DEFAULT_RADIUS_DEG = 10.0
RADIUS_SLIDER_STEPS_PER_DEG = 10
EVALUATE_DELAY_MS = 80
NAME_PATTERN = r"[A-Za-z0-9][A-Za-z0-9._-]*"


class SelectionPanel(QWidget):
    command_issued = Signal(str)  # equivalent CLI command of a confirmed action
    selection_created = Signal(str)

    def __init__(self, view_page) -> None:
        super().__init__()
        self.view = view_page
        self.draft = None
        self._busy_writing = False
        self._stale = False  # centre or radius changed; the shown count is not the current one yet
        self._task = None
        self._runner = LatestWinsRunner(self._evaluated, self._evaluation_failed)
        self.setMinimumWidth(380)
        self.setMaximumWidth(520)

        intro = QLabel(
            "Selections use every particle of the landscape (exact). "
            "Display filters and the 3D sample only change what is drawn."
        )
        intro.setWordWrap(True)
        self.space_label = QLabel()
        self.space_label.setObjectName("selection_space")
        self.name_edit = QLineEdit()
        self.name_edit.setObjectName("selection_name")
        self.name_edit.setPlaceholderText("selection name, e.g. region_01")
        self.name_edit.setValidator(QRegularExpressionValidator(QRegularExpression(NAME_PATTERN)))
        self.confirm_button = QPushButton("Confirm selection")
        self.confirm_button.setObjectName("confirm_selection")
        self.confirm_button.setToolTip("Write a standard Selection under selections/NAME (never overwrites).")
        self.command_label = QLabel()
        self.command_label.setObjectName("selection_command")
        self.command_label.setWordWrap(True)
        self.command_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.command_label.setStyleSheet("font-family: monospace; color: #444;")
        self.saved_list = QListWidget()
        self.saved_list.setObjectName("saved_selections")
        self.saved_list.setMaximumHeight(140)

        self.mode_combo = QComboBox()
        for key, label in MODES:
            self.mode_combo.addItem(label, key)
        self.forms = QStackedWidget()
        self.forms.addWidget(self._radius_form())
        self.forms.addWidget(self._threshold_form())
        self.forms.addWidget(self._range_form())
        self.forms.addWidget(self._random_form())
        self.forms.addWidget(self._metadata_form())

        name_row = QHBoxLayout()
        name_row.addWidget(self.name_edit, 1)
        name_row.addWidget(self.confirm_button)
        layout = QVBoxLayout(self)
        layout.addWidget(intro)
        layout.addWidget(self.space_label)
        layout.addWidget(self.mode_combo)
        layout.addWidget(self.forms)
        layout.addLayout(name_row)
        layout.addWidget(QLabel("Equivalent command:"))
        layout.addWidget(self.command_label)
        layout.addWidget(QLabel("Selections in this run:"))
        layout.addWidget(self.saved_list)
        layout.addStretch(1)

        self._evaluate_timer = QTimer(self)
        self._evaluate_timer.setSingleShot(True)
        self._evaluate_timer.setInterval(EVALUATE_DELAY_MS)
        self._evaluate_timer.timeout.connect(self.evaluate)
        self.mode_combo.currentIndexChanged.connect(self._mode_changed)
        self.name_edit.textChanged.connect(lambda _: (self._update_texts_command(), self._update_buttons()))
        self.confirm_button.clicked.connect(self.confirm)
        self.view.center_picked.connect(self._center_picked)
        self.view.landscape_loaded.connect(self._landscape_loaded)
        self.view.representation_combo.currentIndexChanged.connect(lambda _: self._representation_changed())
        self._update_buttons()

    # --- forms -----------------------------------------------------------------

    def _radius_form(self) -> QWidget:
        box = QGroupBox("Radius (SO(3) geodesic)")
        form = QFormLayout(box)
        hint = QLabel("Click a particle in the 2D or 3D view to set the centre (the nearest displayed particle), "
                      "or type it.")
        hint.setWordWrap(True)
        form.addRow(hint)
        self.center_labels: list[QLabel] = []
        self.center_spins: list[QDoubleSpinBox] = []
        for _ in range(3):
            label, spin = QLabel(), QDoubleSpinBox()
            spin.setRange(-1e4, 1e4)
            spin.setDecimals(4)
            spin.setKeyboardTracking(False)
            spin.valueChanged.connect(lambda _: self._center_typed())
            self.center_labels.append(label)
            self.center_spins.append(spin)
            form.addRow(label, spin)
        self.radius_spin = QDoubleSpinBox()
        self.radius_spin.setObjectName("radius")
        self.radius_spin.setRange(0.1, 180.0)
        self.radius_spin.setDecimals(1)
        self.radius_spin.setSuffix(" deg")
        self.radius_spin.setValue(DEFAULT_RADIUS_DEG)
        self.radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.radius_slider.setRange(1, 180 * RADIUS_SLIDER_STEPS_PER_DEG)
        self.radius_slider.setValue(int(DEFAULT_RADIUS_DEG * RADIUS_SLIDER_STEPS_PER_DEG))
        self.radius_spin.valueChanged.connect(self._radius_spin_changed)
        self.radius_slider.valueChanged.connect(self._radius_slider_changed)
        form.addRow("Radius", self.radius_spin)
        form.addRow(self.radius_slider)
        self.count_label = QLabel("Set a centre to see the exact count.")
        self.count_label.setObjectName("selection_count")
        self.count_label.setWordWrap(True)
        form.addRow(self.count_label)
        self.dim_check = QCheckBox("Dim unselected particles")
        self.dim_check.setChecked(True)
        self.dim_check.toggled.connect(lambda on: self.view.set_dim_unselected(on))
        self.clear_button = QPushButton("Clear")
        self.clear_button.clicked.connect(self.clear_draft)
        row = QHBoxLayout()
        row.addWidget(self.dim_check)
        row.addStretch(1)
        row.addWidget(self.clear_button)
        form.addRow(row)
        self._center_set = False
        self._spin_representation = self.view.representation  # the coordinates the centre boxes hold
        return box

    def _threshold_form(self) -> QWidget:
        box = QGroupBox("SLD threshold (sld_raw, inclusive)")
        form = QFormLayout(box)
        self.sld_min_check, self.sld_min = QCheckBox("Minimum"), QDoubleSpinBox()
        self.sld_max_check, self.sld_max = QCheckBox("Maximum"), QDoubleSpinBox()
        for check, spin, value in ((self.sld_min_check, self.sld_min, 2.0), (self.sld_max_check, self.sld_max, 100.0)):
            spin.setRange(0.0, 1e9)
            spin.setDecimals(3)
            spin.setValue(value)
            check.toggled.connect(lambda _: self._form_changed())
            spin.valueChanged.connect(lambda _: self._form_changed())
            form.addRow(check, spin)
        self.sld_min_check.setChecked(True)
        note = QLabel("Uses the scientific sld_raw values, independent of the display filter.")
        note.setWordWrap(True)
        form.addRow(note)
        return box

    def _range_form(self) -> QWidget:
        box = QGroupBox("Coordinate range (inclusive, in the view's coordinates)")
        form = QFormLayout(box)
        self.range_rows: list[tuple[QCheckBox, QDoubleSpinBox, QDoubleSpinBox]] = []
        for _ in range(3):
            check, low, high = QCheckBox(), QDoubleSpinBox(), QDoubleSpinBox()
            for spin in (low, high):
                spin.setRange(-1e4, 1e4)
                spin.setDecimals(3)
                spin.valueChanged.connect(lambda _: self._form_changed())
            check.toggled.connect(lambda _: self._form_changed())
            row = QHBoxLayout()
            row.addWidget(low)
            row.addWidget(QLabel("to"))
            row.addWidget(high)
            form.addRow(check, row)
            self.range_rows.append((check, low, high))
        note = QLabel("Bounds of several axes are intersected. Euler bounds with lower > upper wrap across ±180°.")
        note.setWordWrap(True)
        form.addRow(note)
        return box

    def _random_form(self) -> QWidget:
        box = QGroupBox("Random fraction of all particles")
        form = QFormLayout(box)
        self.random_percent = QDoubleSpinBox()
        self.random_percent.setRange(0.01, 100.0)
        self.random_percent.setDecimals(2)
        self.random_percent.setSuffix(" %")
        self.random_percent.setValue(10.0)
        self.seed_check, self.seed_spin = QCheckBox("Seed"), QSpinBox()
        self.seed_spin.setRange(0, 2**31 - 1)
        self.seed_check.setToolTip("Without a seed, one is generated and recorded in the selection.")
        for widget in (self.random_percent, self.seed_spin):
            widget.valueChanged.connect(lambda _: self._form_changed())
        self.seed_check.toggled.connect(lambda _: self._form_changed())
        form.addRow("Fraction", self.random_percent)
        form.addRow(self.seed_check, self.seed_spin)
        return box

    def _metadata_form(self) -> QWidget:
        box = QGroupBox("Source metadata values (run-time ref/mov files)")
        form = QFormLayout(box)
        self.meta_domain = QComboBox()
        self.meta_domain.addItems(["ref", "mov"])
        self.meta_column = QLineEdit()
        self.meta_column.setPlaceholderText("e.g. alignments3D/class or rlnClassNumber")
        self.meta_values = QLineEdit()
        self.meta_values.setPlaceholderText("comma-separated, e.g. 1,3")
        self.meta_split = QCheckBox("One selection per value (NAME_VALUE)")
        for widget in (self.meta_column, self.meta_values):
            widget.textChanged.connect(lambda _: self._form_changed())
        self.meta_domain.currentIndexChanged.connect(lambda _: self._form_changed())
        self.meta_split.toggled.connect(lambda on: (self.meta_values.setEnabled(not on), self._form_changed()))
        form.addRow("Domain", self.meta_domain)
        form.addRow("Column", self.meta_column)
        form.addRow("Values", self.meta_values)
        form.addRow(self.meta_split)
        return box

    # --- state -----------------------------------------------------------------

    @property
    def mode(self) -> str:
        return self.mode_combo.currentData()

    def on_shown(self, shown: bool) -> None:
        self.view.picking = shown and self.mode == "radius"
        if shown:
            self.refresh_saved()
            self._update_texts()
        if not shown or self.mode != "radius":
            self.view.set_highlight(None)
        elif self.draft is not None:
            self.view.set_highlight(self.draft.selected_mask, self.draft.center_rv)

    def _landscape_loaded(self, landscape) -> None:
        self._spin_representation = self.view.representation
        self.clear_draft()
        self._update_texts()
        self.refresh_saved()

    def _mode_changed(self, _index: int) -> None:
        self.forms.setCurrentIndex(self.mode_combo.currentIndex())
        self.on_shown(self.isVisible())
        self._update_buttons()

    def _representation_changed(self) -> None:
        landscape = self.view.landscape
        old, new = self._spin_representation, self.view.representation
        self._spin_representation = new
        if landscape is None:
            return
        self._update_texts()
        if self._center_set and old != new:
            # the centre in the boxes (typed, clicked or counted), in the coordinates it was entered in,
            # shown in the new ones; re-evaluated so Confirm records them
            values = tuple(spin.value() for spin in self.center_spins)
            self._set_center(landscape.represent(landscape.rotvec_of(values, old), new))
        elif not self._center_set:
            self._set_center((0.0, 0.0, 0.0), mark_set=False)

    def _update_texts(self) -> None:
        landscape = self.view.landscape
        if landscape is None:
            self.space_label.setText("Open a run to select particles.")
            return
        where = "raw" if landscape.space == "raw" else f"canonical frame {landscape.canonical_id}"
        self.space_label.setText(f"Space: <b>{where}</b> (selections use the space shown)")
        names = landscape.axis_names(self.view.representation)
        unit = "deg" if self.view.representation == "euler" else "rad"
        for label, name in zip(self.center_labels, names):
            label.setText(f"Centre {name} ({unit})")
        for (check, low, high), name, (lo, hi) in zip(
            self.range_rows, names, landscape.axis_limits(self.view.representation)
        ):
            check.setText(f"{name} ({unit})")
            for spin, value in ((low, lo), (high, hi)):
                spin.blockSignals(True)
                spin.setValue(value)
                spin.blockSignals(False)
        self._form_changed()

    def _update_buttons(self) -> None:
        has_name = bool(self.name_edit.text()) and self.name_edit.hasAcceptableInput()
        if self.mode == "radius":
            ready = self.draft is not None and not self._runner.busy and not self._stale
        else:
            ready = self.view.landscape is not None and self._form_options() is not None
        self.confirm_button.setEnabled(has_name and ready and not self._busy_writing)
        self.confirm_button.setText("Confirm selection" if self.mode == "radius" else "Create selection")

    # --- radius ------------------------------------------------------------------

    def _set_center(self, values, *, mark_set: bool = True) -> None:
        for spin, value in zip(self.center_spins, values):
            spin.blockSignals(True)
            spin.setValue(float(value))
            spin.blockSignals(False)
        if mark_set:
            self._center_set = True
            self._request_evaluation()

    def _center_picked(self, row: int) -> None:
        if self.mode != "radius" or self.view.landscape is None:
            return
        coords = self.view.landscape.coordinates(self.view.representation, [row])[0]
        self._set_center(coords)

    def _center_typed(self) -> None:
        self._center_set = True
        self._request_evaluation()

    def _request_evaluation(self) -> None:
        """Mark the shown count out of date at once (Confirm waits), then count after a short pause."""

        self._stale = True
        self.count_label.setText("Counting...")
        self._update_buttons()
        self._evaluate_timer.start()

    def _radius_spin_changed(self, value: float) -> None:
        self.radius_slider.blockSignals(True)
        self.radius_slider.setValue(int(round(value * RADIUS_SLIDER_STEPS_PER_DEG)))
        self.radius_slider.blockSignals(False)
        if self._center_set:
            self._request_evaluation()

    def _radius_slider_changed(self, value: int) -> None:
        self.radius_spin.setValue(value / RADIUS_SLIDER_STEPS_PER_DEG)

    def evaluate(self) -> None:
        landscape = self.view.landscape
        if landscape is None or self._busy_writing or not self._center_set:
            return
        center = tuple(spin.value() for spin in self.center_spins)
        representation = self.view.representation
        radius = self.radius_spin.value()
        self.count_label.setText("Counting...")
        self._runner.submit(
            lambda: (landscape, landscape.evaluate_radius(center, representation=representation, radius_deg=radius))
        )
        self._update_buttons()

    def _evaluated(self, result) -> None:
        landscape, draft = result
        if landscape is not self.view.landscape:
            return  # the run or space changed while counting
        self.draft = draft
        self._stale = self._evaluate_timer.isActive()  # a newer change is already waiting
        shown = 0
        if self.mode == "radius" and self.isVisible():
            self.view.set_highlight(draft.selected_mask, draft.center_rv)
            shown = self.view.highlighted_display_count()
        self.count_label.setText(
            f"<b>{draft.selected_count:,}</b> of {draft.total_count:,} particles "
            f"({draft.selected_fraction:.1%}) within {draft.radius_deg:g}° - exact, full landscape, "
            f"SO(3) geodesic. {shown:,} of them pass the display filter and are highlighted."
        )
        self._update_texts_command()
        self._update_buttons()

    def _evaluation_failed(self, exc: BaseException) -> None:
        self.draft = None
        self._stale = False
        message = exc.message if isinstance(exc, CryoroleError) else str(exc)
        self.count_label.setText(f"Cannot count: {message}")
        self.view.set_highlight(None)
        self._update_buttons()

    def clear_draft(self) -> None:
        self._runner.discard()
        self._evaluate_timer.stop()
        self.draft = None
        self._stale = False
        self._center_set = False
        self._set_center((0.0, 0.0, 0.0), mark_set=False)
        self.count_label.setText("Set a centre to see the exact count.")
        self.view.set_highlight(None)
        self._update_texts_command()
        self._update_buttons()

    # --- other modes -----------------------------------------------------------------

    def _form_options(self) -> dict | None:
        """``api.select`` options for the current non-radius form, or None while incomplete."""

        mode = self.mode
        if mode == "threshold":
            options = {}
            if self.sld_min_check.isChecked():
                options["sld_min"] = self.sld_min.value()
            if self.sld_max_check.isChecked():
                options["sld_max"] = self.sld_max.value()
            return options or None
        if mode == "range":
            bounds = [
                f"{check.text().split()[0]}:{low.value():g}:{high.value():g}"
                for check, low, high in self.range_rows
                if check.isChecked()
            ]
            return {"range_bound": tuple(bounds)} if bounds else None
        if mode == "random":
            options = {"fraction": self.random_percent.value() / 100.0}
            if self.seed_check.isChecked():
                options["seed"] = self.seed_spin.value()
            return options
        if mode == "metadata":
            column = self.meta_column.text().strip()
            if not column:
                return None
            options = {"metadata_domain": self.meta_domain.currentText(), "metadata_column": column}
            if self.meta_split.isChecked():
                options["split_by_value"] = True
            elif self.meta_values.text().strip():
                options["metadata_value"] = self.meta_values.text().strip()
            else:
                return None
            return options
        return None

    def _form_changed(self) -> None:
        self._update_texts_command()
        self._update_buttons()

    def _cli_command(self, name: str) -> str:
        landscape = self.view.landscape
        if landscape is None:
            return ""
        if self.mode == "radius":
            return self.draft.command(name) if self.draft is not None else ""
        options = self._form_options()
        if options is None:
            return ""
        parts = ["cryorole", "select", "--run-dir", str(landscape.run_dir), "--selection-id", name,
                 "--mode", self.mode, "--space", landscape.space]
        if landscape.space == "canonical":
            parts += ["--canonical-id", landscape.canonical_id]
        for key, value in options.items():
            flag = "--" + key.replace("_", "-")
            if value is True:
                parts.append(flag)
            elif isinstance(value, tuple):
                for item in value:
                    parts += [flag, item]
            else:
                parts += [flag, f"{value:g}" if isinstance(value, float) else str(value)]
        return shlex.join(parts)

    def _update_texts_command(self) -> None:
        self.command_label.setText(self._cli_command(self.name_edit.text() or "NAME"))

    # --- writing -----------------------------------------------------------------

    def confirm(self) -> None:
        landscape = self.view.landscape
        name = self.name_edit.text()
        if landscape is None or not name:
            return
        if self.mode == "radius":
            draft = self.draft
            if draft is None or self._runner.busy or self._stale:
                return
            command = draft.command(name)
            frame = self.view.frame  # what was shown: recorded as display provenance

            def write():
                return landscape.confirm(name, draft, frame=frame)
        else:
            options = self._form_options()
            if options is None:
                return
            if name in landscape.selection_ids():
                show_error(self, CryoroleError(
                    "output_exists", f"A selection named {name} already exists in this run.",
                    "Choose another name; the GUI never overwrites a selection.",
                ))
                return
            command = self._cli_command(name)
            mode, space, frame = self.mode, landscape.space, landscape.canonical_id

            def write():
                return api.select(landscape.run_dir, selection_id=name, mode=mode, space=space,
                                  canonical_id=frame, **options)

        self._busy_writing = True
        self._set_inputs_enabled(False)
        self._task = run_in_background(
            write,
            lambda result: self._written(name, command, result),
            self._write_failed,
        )

    def _written(self, name: str, command: str, result) -> None:
        self._busy_writing = False
        self._task = None
        self._set_inputs_enabled(True)
        if isinstance(result, dict):  # radius confirm returns the selection summary
            counts = [int(result.get("selected_count", 0))]
        else:
            counts = list(result.selected_counts)
        total = sum(counts)
        notify_user(f"Selection {name}: {total:,} particles written ({len(counts)} selection(s)).")
        self.command_issued.emit(command)
        self.selection_created.emit(name)
        self.name_edit.clear()
        self.refresh_saved()
        self._update_buttons()

    def _write_failed(self, exc: BaseException) -> None:
        self._busy_writing = False
        self._task = None
        self._set_inputs_enabled(True)
        self._update_buttons()
        show_error(self, exc, title="cryoROLE: selection not written")

    def _set_inputs_enabled(self, enabled: bool) -> None:
        for widget in (self.mode_combo, self.forms, self.name_edit):
            widget.setEnabled(enabled)
        self.view.picking = enabled and self.mode == "radius" and self.isVisible()

    def refresh_saved(self) -> None:
        landscape = self.view.landscape
        if landscape is None:
            self.saved_list.clear()
            return
        self.saved_list.clear()
        for name, count in landscape.saved_selections():
            self.saved_list.addItem(name + (f"  ({count:,} particles)" if count is not None else ""))
