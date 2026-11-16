"""Main window of the cryoROLE GUI."""

from __future__ import annotations

import logging
from pathlib import Path

from PySide6.QtCore import QObject, Qt, Signal
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication,
    QMessageBox,
    QDockWidget,
    QFileDialog,
    QLabel,
    QListWidget,
    QMainWindow,
    QPlainTextEdit,
    QSplitter,
    QStackedWidget,
)

from cryorole.gui.glprobe import ProbeResult
from cryorole.gui.pages.export import ExportPage
from cryorole.gui.pages.run_setup import InputsPage, ParametersPage, RunPage, RunSetup
from cryorole.gui.pages.view import ViewPage
from cryorole.logs import LOGGER_NAME

# Page order follows the workflow.
PAGES = ("Inputs", "Parameters", "Run", "View", "Select", "Export")


class _LogBridge(QObject):
    message = Signal(str)


class _QtLogHandler(logging.Handler):
    """Forward ``cryorole`` log records to the log panel (thread-safe via a queued signal)."""

    def __init__(self) -> None:
        super().__init__(level=logging.INFO)
        self.bridge = _LogBridge()
        self.setFormatter(logging.Formatter("%(levelname)s  %(message)s"))

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self.bridge.message.emit(self.format(record))
        except RuntimeError:  # the panel was already deleted
            pass


class MainWindow(QMainWindow):
    def __init__(self, probe: ProbeResult) -> None:
        super().__init__()
        self.probe = probe
        self.setWindowTitle("cryoROLE")
        self.resize(1400, 900)

        self.nav = QListWidget()
        self.nav.addItems(list(PAGES))
        self.nav.setMaximumWidth(180)
        self.pages = QStackedWidget()
        self.view_page = ViewPage(gl_available=probe.gl_available, gl_reason=probe.reason)
        self.view_page.setObjectName("page_view")
        self.view_page.landscape_loaded.connect(self._landscape_loaded)
        # View and Select share one widget: selecting needs the views; Select adds its panel.
        self.export_page = ExportPage(self.view_page)
        self.export_page.setObjectName("page_export")
        self.run_setup = RunSetup()
        self.inputs_page = InputsPage(self.run_setup)
        self.parameters_page = ParametersPage(self.run_setup)
        self.run_page = RunPage(self.run_setup)
        self.page_widgets = {
            "Inputs": self.inputs_page, "Parameters": self.parameters_page, "Run": self.run_page,
            "View": self.view_page, "Select": self.view_page, "Export": self.export_page,
        }
        for widget in (self.inputs_page, self.parameters_page, self.run_page, self.view_page, self.export_page):
            self.pages.addWidget(widget)
        self.run_page.run_finished.connect(self.open_run)
        self.nav.currentRowChanged.connect(self._nav_changed)
        self.nav.setCurrentRow(0)

        splitter = QSplitter()
        splitter.addWidget(self.nav)
        splitter.addWidget(self.pages)
        splitter.setStretchFactor(1, 1)
        self.setCentralWidget(splitter)

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setMaximumBlockCount(5000)
        dock = QDockWidget("Log", self)
        dock.setObjectName("log_dock")
        dock.setWidget(self.log_view)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, dock)
        self.resizeDocks([dock], [110], Qt.Orientation.Vertical)
        self.commands_view = QPlainTextEdit()
        self.commands_view.setReadOnly(True)
        self.commands_view.setObjectName("commands")
        self.commands_view.setPlaceholderText("Equivalent cryorole commands of confirmed actions appear here.")
        commands_dock = QDockWidget("Commands", self)
        commands_dock.setObjectName("commands_dock")
        commands_dock.setWidget(self.commands_view)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, commands_dock)
        self.tabifyDockWidget(dock, commands_dock)
        dock.raise_()
        self.view_page.selection_panel.command_issued.connect(self.add_command)
        self.view_page.selection_panel.selection_created.connect(self.export_page.refresh_selections)
        self.export_page.command_issued.connect(self.add_command)
        self.run_page.command_issued.connect(self.add_command)
        self.view_page.command_issued.connect(self.add_command)
        self._log_handler = _QtLogHandler()
        self._log_handler.bridge.message.connect(self.log_view.appendPlainText)
        self._logger = logging.getLogger(LOGGER_NAME)
        self._logger.addHandler(self._log_handler)

        self.gl_label = QLabel(probe.summary())
        self.gl_label.setObjectName("gl_status")
        if probe.remedy and probe.gl_status != "ok":
            self.gl_label.setToolTip(probe.remedy)
        self.statusBar().addPermanentWidget(self.gl_label)

        file_menu = self.menuBar().addMenu("&File")
        open_action = QAction("&Open run...", self)
        open_action.setShortcut(QKeySequence.StandardKey.Open)
        open_action.triggered.connect(self.choose_run)
        file_menu.addAction(open_action)
        file_menu.addSeparator()
        quit_action = QAction("&Quit", self)
        quit_action.setShortcut(QKeySequence.StandardKey.Quit)
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

    def _nav_changed(self, row: int) -> None:
        name = PAGES[row]
        self.pages.setCurrentWidget(self.page_widgets[name])
        self.view_page.set_selecting(name == "Select")
        if name == "Export":
            self.export_page.on_shown()

    def add_command(self, command: str) -> None:
        self.commands_view.appendPlainText(command)

    def choose_run(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, "Open a cryoROLE run bundle", str(Path.cwd()))
        if folder:
            self.open_run(folder)

    def open_run(self, run_dir: str | Path) -> None:
        if PAGES[self.nav.currentRow()] not in ("View", "Select"):
            self.nav.setCurrentRow(self.page_index("View"))
        self.view_page.open_run(run_dir)

    @staticmethod
    def page_index(name: str) -> int:
        return PAGES.index(name)

    def _landscape_loaded(self, landscape) -> None:
        self.setWindowTitle(f"cryoROLE - {landscape.run_dir}")

    @property
    def gl_available(self) -> bool:
        return self.probe.gl_available

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt override
        if self.view_page.saving_figures:
            QMessageBox.information(self, "Saving figures", "Figures are still being saved. Wait for the result before closing this window.")
            event.ignore()
            return
        busy = []  # asked about together: nothing is stopped unless every answer is yes
        if self.run_page.worker.running:
            busy.append("A run is in progress. Stop it and quit? Nothing will be published.")
        if self.view_page.canonicalizing:
            busy.append("A canonical frame is being written. Stop it and quit? The unfinished frame is removed.")
        for question in busy:
            if QMessageBox.question(self, "cryoROLE", question) != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
        self.run_page.worker.kill()
        self.view_page.stop_canonicalize()
        self.view_page.close_3d()
        self._logger.removeHandler(self._log_handler)
        super().closeEvent(event)


def run_app(probe: ProbeResult, *, run_dir: str | Path | None = None, argv: list[str] | None = None) -> int:
    app = QApplication.instance() or QApplication(argv or ["cryorole"])
    app.setApplicationName("cryoROLE")
    window = MainWindow(probe)
    window.show()
    if run_dir is not None:
        window.open_run(run_dir)
    return app.exec()
