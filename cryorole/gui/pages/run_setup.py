"""Inputs, Parameters and Run pages: from two refinement outputs to a run bundle.

* Inputs: a reference and a moving body, each a .cs/.star file or a job folder
  (`api.inspect_input` picks the particle file only when unambiguous), then the
  same preflight as `cryorole preflight` (nothing is written).
* Parameters: output folder and the few run options, each with its CLI default.
* Run: `cryorole run` in a child process (`cryorole.gui.worker`) with stage
  progress and Cancel; a cancelled run publishes nothing. The new bundle opens
  in View.
"""

from __future__ import annotations

from pathlib import Path
import shlex
import time

from PySide6.QtCore import QObject, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from cryorole import api
from cryorole.errors import CryoroleError
from cryorole.gui.process import WorkerProcess
from cryorole.gui.tasks import LatestWinsRunner
from cryorole.gui.widgets.errors import show_error
from cryorole.logs import notify_user

DOMAINS = (("ref", "Reference body (ref)"), ("mov", "Moving body (mov)"))
# Keep the default equal to cryorole.core.density.DEFAULT_SLD_METRIC (the GUI imports only cryorole.api).
DEFAULT_SLD_METRIC = "so3_geodesic"
SLD_METRICS = (("so3_geodesic", "SO(3) geodesic (default)"),
               ("rotvec_euclidean", "Euclidean, rotation-vector space (faster for very large data)"))
DEFAULT_K = 50
PREFLIGHT_DELAY_MS = 400
CANCEL_GRACE_S = 20
READINESS_STYLE = {"READY": "#1e7b34", "READY_WITH_WARNINGS": "#a15c00", "BLOCKED": "#b00020"}
WARNING_TEXT = {  # preflight warning codes that are not sentences yet
    "large_input_consider_--no-visualize": (
        "Large input: consider turning off the quick-look figures (Parameters) to save time and memory."
    ),
}


class RunSetup(QObject):
    """Shared state of the three pages."""

    inputs_changed = Signal()
    options_changed = Signal()
    preflight_changed = Signal()

    def __init__(self) -> None:
        super().__init__()
        self.inspections: dict[str, object | None] = {"ref": None, "mov": None}
        self.chosen: dict[str, Path | None] = {"ref": None, "mov": None}
        self.output_dir = Path.cwd() / "cryorole_outputs"
        self.sld_metric = DEFAULT_SLD_METRIC
        self.k_neighbors = DEFAULT_K
        self.visualize = True
        self.raw_csv = True
        self.point_style = "scatter"
        self.preflight: dict | None = None
        self.preflight_error: str | None = None
        self.preflight_running = False
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(PREFLIGHT_DELAY_MS)
        self._timer.timeout.connect(self._run_preflight)
        self._runner = LatestWinsRunner(self._preflight_done, self._preflight_failed)
        self.inputs_changed.connect(self.request_preflight)
        self.options_changed.connect(self.request_preflight)

    # --- inputs ------------------------------------------------------------------

    def set_inspection(self, domain: str, inspection) -> None:
        self.inspections[domain] = inspection
        chosen = inspection.chosen if inspection is not None else None
        self.chosen[domain] = chosen.path if chosen is not None else None
        self.inputs_changed.emit()

    def choose(self, domain: str, path: Path) -> None:
        self.chosen[domain] = path
        self.inputs_changed.emit()

    def software(self, domain: str) -> str | None:
        inspection = self.inspections[domain]
        if inspection is None or self.chosen[domain] is None:
            return None
        for candidate in inspection.candidates:
            if candidate.path == self.chosen[domain]:
                return candidate.software
        return None

    @property
    def mixed_software(self) -> bool:
        kinds = {self.software("ref"), self.software("mov")}
        return None not in kinds and len(kinds) == 2

    @property
    def inputs_ready(self) -> bool:
        return self.chosen["ref"] is not None and self.chosen["mov"] is not None and not self.mixed_software

    # --- options -----------------------------------------------------------------

    def run_options(self) -> dict:
        options = {"k_neighbors": self.k_neighbors, "sld_metric": self.sld_metric}
        if not self.visualize:
            options["no_visualize"] = True
        if not self.raw_csv:
            options["no_raw_csv"] = True
        if self.point_style != "scatter":
            options["point_style"] = self.point_style
        return options

    def cli_command(self) -> str:
        if not self.inputs_ready:
            return ""
        parts = ["cryorole", "run", "--ref", str(self.chosen["ref"]), "--mov", str(self.chosen["mov"]),
                 "--output-dir", str(self.output_dir)]
        if self.sld_metric != DEFAULT_SLD_METRIC:
            parts += ["--sld-metric", self.sld_metric]
        if self.k_neighbors != DEFAULT_K:
            parts += ["--k-neighbors", str(self.k_neighbors)]
        if not self.visualize:
            parts.append("--no-visualize")
        if not self.raw_csv:
            parts.append("--no-raw-csv")
        if self.point_style != "scatter":
            parts += ["--point-style", self.point_style]
        return shlex.join(parts)

    # --- preflight -----------------------------------------------------------------

    def request_preflight(self) -> None:
        self._timer.stop()
        self._runner.discard()
        self.preflight, self.preflight_error = None, None
        self.preflight_running = self.inputs_ready
        if self.inputs_ready:
            self._timer.start()
        else:  # nothing to check; a check still running answers for older inputs and is dropped
            self.preflight_running = False
            self._timer.stop()
        self.preflight_changed.emit()

    def _preflight_key(self) -> tuple:
        """What a preflight result answers: the inputs and every option it checks."""

        return (self.chosen["ref"], self.chosen["mov"], str(self.output_dir), self.k_neighbors, self.sld_metric,
                self.visualize, self.raw_csv)

    def _run_preflight(self) -> None:
        if not self.inputs_ready:
            self.preflight_running = False
            self.preflight_changed.emit()
            return
        ref, mov = self.chosen["ref"], self.chosen["mov"]
        options = {"output_dir": str(self.output_dir), "k_neighbors": self.k_neighbors,
                   "sld_metric": self.sld_metric, "visualize": self.visualize, "raw_csv": self.raw_csv}
        key = self._preflight_key()

        def check():
            try:
                return key, api.preflight(ref, mov, **options).report, None
            except Exception as exc:  # noqa: BLE001 - reported with the key it belongs to
                return key, None, exc

        self._runner.submit(check)

    def _stale(self, key: tuple) -> bool:
        # an answer for other inputs or options, or one a newer check (still waiting) will replace
        return key != self._preflight_key() or self._timer.isActive() or self._runner.busy

    def _preflight_done(self, result) -> None:
        key, report, exc = result
        if self._stale(key):
            return
        self.preflight_running = False
        if exc is not None:
            message = exc.message if isinstance(exc, CryoroleError) else str(exc)
            self.preflight, self.preflight_error = None, message
        else:
            self.preflight, self.preflight_error = report, None
        self.preflight_changed.emit()

    def _preflight_failed(self, exc: BaseException) -> None:  # check() returns its errors; this is a bug guard
        if self._timer.isActive() or self._runner.busy:
            return
        self.preflight_running = False
        message = exc.message if isinstance(exc, CryoroleError) else str(exc)
        self.preflight, self.preflight_error = None, message
        self.preflight_changed.emit()

    @property
    def readiness(self) -> str | None:
        return None if self.preflight is None else self.preflight.get("readiness")


# --- Inputs ----------------------------------------------------------------------------


class _InputBox(QGroupBox):
    def __init__(self, page: "InputsPage", domain: str, title: str) -> None:
        super().__init__(title)
        self.page, self.domain = page, domain
        self.path_edit = QLineEdit()
        self.path_edit.setReadOnly(True)
        self.path_edit.setPlaceholderText("Choose a .cs / .star file or a job folder, or drop one here")
        file_button, folder_button = QPushButton("File..."), QPushButton("Job folder...")
        file_button.clicked.connect(self.choose_file)
        folder_button.clicked.connect(self.choose_folder)
        self.candidates = QComboBox()
        self.candidates.setVisible(False)
        self.candidates.activated.connect(self._candidate_chosen)
        self.info = QLabel()
        self.info.setWordWrap(True)
        self.info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        row = QHBoxLayout()
        row.addWidget(self.path_edit, 1)
        row.addWidget(file_button)
        row.addWidget(folder_button)
        layout = QVBoxLayout(self)
        layout.addLayout(row)
        layout.addWidget(self.candidates)
        layout.addWidget(self.info)
        self._inspection_runner = LatestWinsRunner(self._inspected, self._failed)
        self.destroyed.connect(self._inspection_runner.discard)

    def choose_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, self.title(), str(Path.cwd()),
                                              "Particle metadata (*.cs *.star);;All files (*)")
        if path:
            self.set_path(path)

    def choose_folder(self) -> None:
        path = QFileDialog.getExistingDirectory(self, self.title(), str(Path.cwd()))
        if path:
            self.set_path(path)

    def set_path(self, path: str | Path) -> None:
        self.path_edit.setText(str(path))
        self.info.setText("Looking for the particle file...")
        self.candidates.clear()
        self.candidates.setVisible(False)
        self.page.setup.set_inspection(self.domain, None)
        self._inspection_runner.submit(lambda: api.inspect_input(path))

    def _inspected(self, inspection) -> None:
        self.candidates.clear()
        for candidate in inspection.candidates:
            self.candidates.addItem(f"{candidate.path.name} - {candidate.reason}", str(candidate.path))
        many = len(inspection.candidates) > 1
        self.candidates.setVisible(many)
        if inspection.chosen is not None:
            self.candidates.setCurrentIndex(self.candidates.findData(str(inspection.chosen.path)))
        else:
            self.candidates.setCurrentIndex(-1)
            self.candidates.setPlaceholderText(inspection.note or "Choose one")
        self.page.setup.set_inspection(self.domain, inspection)
        self.show_info()

    def _failed(self, exc: BaseException) -> None:
        self.candidates.setVisible(False)
        self.page.setup.set_inspection(self.domain, None)
        message = exc.message if isinstance(exc, CryoroleError) else str(exc)
        self.info.setText(f"<span style='color:#b00020'>{message}</span>")

    def _candidate_chosen(self, index: int) -> None:
        self.page.setup.choose(self.domain, Path(self.candidates.itemData(index)))
        self.show_info()

    def show_info(self) -> None:
        setup = self.page.setup
        inspection, chosen = setup.inspections[self.domain], setup.chosen[self.domain]
        if inspection is None:
            return
        if chosen is None:
            self.info.setText(inspection.note or "Choose one of the particle outputs above.")
            return
        candidate = next(c for c in inspection.candidates if c.path == chosen)
        parts = [f"<b>{candidate.path}</b>", f"{'CryoSPARC' if candidate.software == 'cryosparc' else 'RELION'}"
                 + (f", job {candidate.job}" if candidate.job else "")]
        if candidate.particle_count is not None:
            parts.append(f"{candidate.particle_count:,} particles")
        parts.append(f"chosen: {candidate.reason}")
        if candidate.group is not None:
            parts.append(f"result group {candidate.group.name}")
        self.info.setText(" · ".join(parts))


class InputsPage(QWidget):
    def __init__(self, setup: RunSetup) -> None:
        super().__init__()
        self.setup = setup
        self.setAcceptDrops(True)
        intro = QLabel(
            "Choose the two refinement outputs. The reference body defines the frame of the relative "
            "orientation (RO = R_ref⁻¹ R_mov); swapping them asks a different question."
        )
        intro.setWordWrap(True)
        self.boxes = {domain: _InputBox(self, domain, title) for domain, title in DOMAINS}
        self.readiness = QLabel()
        self.readiness.setObjectName("readiness")
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        recheck = QPushButton("Check again")
        recheck.clicked.connect(setup.request_preflight)
        check_box = QGroupBox("Preflight (the same checks as cryorole preflight; nothing is written)")
        check_layout = QVBoxLayout(check_box)
        top = QHBoxLayout()
        top.addWidget(self.readiness, 1)
        top.addWidget(recheck)
        check_layout.addLayout(top)
        check_layout.addWidget(self.summary)
        layout = QVBoxLayout(self)
        layout.addWidget(intro)
        for box in self.boxes.values():
            layout.addWidget(box)
        layout.addWidget(check_box)
        layout.addStretch(1)
        setup.preflight_changed.connect(self.show_preflight)
        self.show_preflight()

    def dragEnterEvent(self, event) -> None:  # noqa: N802 - Qt override
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:  # noqa: N802 - Qt override
        urls = [url.toLocalFile() for url in event.mimeData().urls() if url.isLocalFile()]
        if not urls:
            return
        point = event.position().toPoint()
        target = next((domain for domain, box in self.boxes.items() if box.geometry().contains(point)), None)
        if target is None:
            target = "ref" if self.setup.chosen["ref"] is None else "mov"
        self.boxes[target].set_path(urls[0])
        if len(urls) > 1 and target == "ref":
            self.boxes["mov"].set_path(urls[1])

    def show_preflight(self) -> None:
        setup = self.setup
        if setup.mixed_software:
            self._badge("BLOCKED", "Both inputs must come from the same software (CryoSPARC or RELION) in this version.")
            self.summary.setText("")
            return
        if not setup.inputs_ready:
            self.readiness.setText("Choose both inputs to run the checks.")
            self.summary.setText("")
            return
        if setup.preflight_error:
            self._badge("BLOCKED", setup.preflight_error)
            self.summary.setText("")
            return
        if setup.preflight is None:
            self.readiness.setText("Checking (reading both inputs)...")
            self.summary.setText("")
            return
        report = setup.preflight
        self._badge(report["readiness"], {"READY": "Ready to run.", "READY_WITH_WARNINGS": "Ready, with warnings.",
                                          "BLOCKED": "Blocked: fix the errors below."}.get(report["readiness"], ""))
        self.summary.setText(preflight_summary_html(report))

    def _badge(self, readiness: str, text: str) -> None:
        color = READINESS_STYLE.get(readiness, "#444")
        self.readiness.setText(f"<b style='color:{color}'>{readiness}</b> {text}")


def preflight_summary_html(report: dict) -> str:
    matching = report.get("matching") or {}
    estimate = report.get("resource_estimate") or {}
    lines = []
    if matching:
        lines.append(
            f"Matched <b>{matching.get('matched_count', 0):,}</b> particles by "
            f"<code>{matching.get('identity_key')}</code>: {matching.get('ref_coverage', 0):.1%} of ref "
            f"({matching.get('ref_row_count', 0):,}), {matching.get('mov_coverage', 0):.1%} of mov "
            f"({matching.get('mov_row_count', 0):,})."
        )
    if estimate:
        lines.append(
            f"Estimated peak memory {estimate.get('estimated_peak_memory_mib', 0) / 1024:.1f} GB, "
            f"bundle {estimate.get('estimated_bundle_mib', 0):,.0f} MB on disk."
        )
    for error in report.get("errors") or []:
        lines.append(f"<span style='color:#b00020'>Error: {error}</span>")
    for warning in report.get("warnings") or []:
        lines.append(f"<span style='color:#a15c00'>Warning: {WARNING_TEXT.get(warning, warning)}</span>")
    return "<br>".join(lines)


# --- Parameters ------------------------------------------------------------------------


class ParametersPage(QWidget):
    def __init__(self, setup: RunSetup) -> None:
        super().__init__()
        self.setup = setup
        self.output_edit = QLineEdit(str(setup.output_dir))
        self.output_edit.setObjectName("run_output")
        browse = QPushButton("Browse...")
        browse.clicked.connect(self.choose_output)
        self.metric_combo = QComboBox()
        for key, label in SLD_METRICS:
            self.metric_combo.addItem(label, key)
        self.k_spin = QSpinBox()
        self.k_spin.setRange(2, 1000)
        self.k_spin.setValue(DEFAULT_K)
        self.k_spin.setToolTip(f"Neighbours used for the local sampling density (SLD). Default {DEFAULT_K}.")
        self.visualize_check = QCheckBox("Write quick-look figures (PNG)")
        self.visualize_check.setChecked(True)
        self.raw_csv_check = QCheckBox("Write the full raw landscape table (data/raw_landscape.csv)")
        self.raw_csv_check.setChecked(True)
        self.raw_csv_check.setToolTip("The NPZ always holds the same data; the CSV is for reading outside cryoROLE.")
        self.point_combo = QComboBox()
        self.point_combo.addItem("scatter (default)", "scatter")
        self.point_combo.addItem("density image (faster for large runs)", "density")
        self.command = QLabel()
        self.command.setObjectName("run_command")
        self.command.setWordWrap(True)
        self.command.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.command.setStyleSheet("font-family: monospace; color: #444;")
        self.output_note = QLabel()
        self.output_note.setWordWrap(True)

        output_row = QHBoxLayout()
        output_row.addWidget(self.output_edit, 1)
        output_row.addWidget(browse)
        form = QFormLayout()
        form.addRow("Output folder", output_row)
        form.addRow("", self.output_note)
        form.addRow("SLD metric", self.metric_combo)
        advanced = QGroupBox("Advanced (defaults are those of cryorole run)")
        advanced_form = QFormLayout(advanced)
        advanced_form.addRow("Neighbours k", self.k_spin)
        advanced_form.addRow(self.visualize_check)
        advanced_form.addRow(self.raw_csv_check)
        advanced_form.addRow("Figure points", self.point_combo)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(advanced)
        layout.addWidget(QLabel("Equivalent command:"))
        layout.addWidget(self.command)
        layout.addStretch(1)

        self.output_edit.editingFinished.connect(self._changed)
        for widget in (self.metric_combo, self.point_combo):
            widget.currentIndexChanged.connect(lambda _: self._changed())
        self.k_spin.valueChanged.connect(lambda _: self._changed())
        for check in (self.visualize_check, self.raw_csv_check):
            check.toggled.connect(lambda _: self._changed())
        setup.inputs_changed.connect(self.refresh)
        setup.preflight_changed.connect(self.refresh)
        self.refresh()

    def showEvent(self, event) -> None:  # noqa: N802 - Qt override; the output folder may have appeared
        super().showEvent(event)
        self.refresh()

    def choose_output(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, "Folder to create the run bundle in", str(Path.cwd()))
        if folder:
            self.output_edit.setText(str(Path(folder) / Path(self.output_edit.text() or "cryorole_outputs").name))
            self._changed()

    def _changed(self) -> None:
        setup = self.setup
        setup.output_dir = Path(self.output_edit.text()).expanduser()
        setup.sld_metric = self.metric_combo.currentData()
        setup.k_neighbors = self.k_spin.value()
        setup.visualize = self.visualize_check.isChecked()
        setup.raw_csv = self.raw_csv_check.isChecked()
        setup.point_style = self.point_combo.currentData()
        setup.options_changed.emit()
        self.refresh()

    def refresh(self) -> None:
        self.command.setText(self.setup.cli_command() or "Choose both inputs first.")
        if self.setup.output_dir.exists():
            self.output_note.setText("<span style='color:#b00020'>This folder already exists; "
                                     "choose a new one (the GUI never overwrites a run).</span>")
        else:
            self.output_note.setText("A new folder; it is created when the run is published.")


# --- Run -------------------------------------------------------------------------------


class RunPage(QWidget):
    run_finished = Signal(str)  # output folder of a published run
    command_issued = Signal(str)

    def __init__(self, setup: RunSetup) -> None:
        super().__init__()
        self.setup = setup
        self.worker = WorkerProcess(self)
        self.started_at: float | None = None
        self._run_output_dir: Path | None = None  # the running run's folder (Parameters may change meanwhile)
        self.cancel_requested_at: float | None = None
        self._command = ""
        self._stage: tuple[int, int] | None = None
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.start_button = QPushButton("Start run")
        self.start_button.setObjectName("start_run")
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.stop_button = QPushButton("Stop now")
        self.stop_button.setVisible(False)
        self.stop_button.setToolTip("Stop the worker at once. Its hidden staging folder may remain next to the output.")
        self.progress = QProgressBar()
        self.progress.setRange(0, 1000)
        self.stage_label = QLabel()
        self.stage_label.setObjectName("run_stage")
        self.elapsed_label = QLabel()
        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setObjectName("run_log")
        buttons = QHBoxLayout()
        for button in (self.start_button, self.cancel_button, self.stop_button):
            buttons.addWidget(button)
        buttons.addStretch(1)
        buttons.addWidget(self.elapsed_label)
        layout = QVBoxLayout(self)
        layout.addWidget(self.summary)
        layout.addLayout(buttons)
        layout.addWidget(self.progress)
        layout.addWidget(self.stage_label)
        layout.addWidget(self.log_view, 1)
        self._clock = QTimer(self)
        self._clock.setInterval(500)
        self._clock.timeout.connect(self._tick)

        self.start_button.clicked.connect(self.start)
        self.cancel_button.clicked.connect(self.cancel)
        self.stop_button.clicked.connect(self.worker.kill)
        self.worker.event.connect(self._event)
        self.worker.log.connect(lambda level, message: self.log_view.appendPlainText(f"{level}  {message}"))
        self.worker.succeeded.connect(self._succeeded)
        self.worker.failed.connect(self._failed)
        self.worker.cancelled.connect(self._cancelled)
        self.worker.stopped.connect(self.refresh)
        setup.inputs_changed.connect(self.refresh)
        setup.options_changed.connect(self.refresh)
        setup.preflight_changed.connect(self.refresh)
        self.refresh()

    def showEvent(self, event) -> None:  # noqa: N802 - Qt override; the output folder may have appeared
        super().showEvent(event)
        self.refresh()

    def can_start(self) -> tuple[bool, str]:
        setup = self.setup
        if self.worker.running:
            return False, "A run is in progress."
        if not setup.inputs_ready:
            return False, "Choose both inputs on the Inputs page."
        if setup.output_dir.exists():
            return False, f"{setup.output_dir} already exists; choose a new output folder on the Parameters page."
        if setup.preflight_running:
            return False, "Waiting for the preflight checks."
        if setup.preflight is None:
            return False, setup.preflight_error or "Preflight has not run."
        if setup.readiness == "BLOCKED":
            return False, "Preflight is BLOCKED; see the Inputs page."
        return True, ""

    def refresh(self) -> None:
        ok, reason = self.can_start()
        self.start_button.setEnabled(ok)
        setup = self.setup
        text = [f"ref: {setup.chosen['ref'] or '-'}", f"mov: {setup.chosen['mov'] or '-'}",
                f"output: {setup.output_dir}"]
        if setup.readiness:
            text.append(f"preflight: {setup.readiness}")
        if reason and not self.worker.running:
            text.append(f"<i>{reason}</i>")
        self.summary.setText("<br>".join(text))

    def start(self) -> None:
        ok, reason = self.can_start()
        if not ok:
            show_error(self, CryoroleError("blocked", reason))
            return
        setup = self.setup
        self._command = setup.cli_command()
        request = {"command": "run", "ref": str(setup.chosen["ref"]), "mov": str(setup.chosen["mov"]),
                   "output_dir": str(setup.output_dir), **setup.run_options()}
        self.log_view.clear()
        self.log_view.appendPlainText(f"$ {self._command}")
        self.progress.setValue(0)
        self.stage_label.setText("Starting...")
        self.started_at, self.cancel_requested_at = time.monotonic(), None
        self._run_output_dir = setup.output_dir
        self.worker.start(request)
        self.cancel_button.setEnabled(True)
        self._clock.start()
        self.refresh()

    def cancel(self) -> None:
        self.worker.cancel()
        self.cancel_requested_at = time.monotonic()
        self.cancel_button.setEnabled(False)
        self.stage_label.setText("Cancelling at the next stage boundary...")

    def _tick(self) -> None:
        if self.started_at is not None:
            self.elapsed_label.setText(f"{time.monotonic() - self.started_at:.0f} s")
        if self.cancel_requested_at is not None and time.monotonic() - self.cancel_requested_at > CANCEL_GRACE_S:
            self.stop_button.setVisible(True)

    def _event(self, event: dict) -> None:
        number, total = event.get("stage_number"), event.get("stage_total")
        if event.get("kind") == "stage_start" and number and total:
            self._stage = (number, total)
            self.stage_label.setText(f"Stage {number} of {total}: {event.get('stage')}")
            self.progress.setValue(int(1000 * (number - 1) / total))
        elif event.get("kind") == "batch" and event.get("fraction") is not None and self._stage:
            number, total = self._stage
            self.progress.setValue(int(1000 * (number - 1 + event["fraction"]) / total))
        if event.get("kind") in ("stage_start", "stage_done", "info"):  # warnings arrive through the log
            self.log_view.appendPlainText(event.get("message", ""))

    def _finish(self) -> None:
        self._clock.stop()
        self._tick()
        self.cancel_button.setEnabled(False)
        self.stop_button.setVisible(False)
        self.refresh()

    def _succeeded(self, result: dict) -> None:
        self.progress.setValue(1000)
        self.stage_label.setText(f"Done: {result.get('matched_count', 0):,} particles in {result.get('output_dir')}")
        self._finish()
        notify_user(f"Run published: {result.get('output_dir')}")
        self.command_issued.emit(self._command)
        self.run_finished.emit(result["output_dir"])

    def _failed(self, error: dict) -> None:
        self.stage_label.setText("Failed.")
        self._finish()
        message = error.get("message", "")
        if error.get("code") == "worker_stopped":
            out = self._run_output_dir or self.setup.output_dir
            leftovers = sorted(out.parent.glob(f".{out.name}.preparing-*"))
            if leftovers:
                message += " A staging folder remains and can be deleted: " + ", ".join(str(p) for p in leftovers)
        show_error(self, CryoroleError(error.get("code", "internal"), message, error.get("remedy")),
                   title="cryoROLE: run failed")

    def _cancelled(self) -> None:
        self.stage_label.setText("Cancelled; nothing was published.")
        self.progress.setValue(0)
        self._finish()
