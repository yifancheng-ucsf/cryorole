"""Export page: write a selection's source-file subsets and say how to use them.

Calls `api.export` (verbatim row subsets of the run's source STAR/CS files, by
recorded source rows). For CryoSPARC sources the source job's result group
(.csg) is mirrored so the subset can go back through Import Result Group; the
page shows the copy commands and the import path. The GUI never overwrites an
export folder.
"""

from __future__ import annotations

from html import escape
from pathlib import Path
import shlex

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from cryorole import api
from cryorole.errors import CryoroleError
from cryorole.gui.tasks import LatestWinsRunner, run_in_background
from cryorole.gui.widgets.errors import show_error
from cryorole.logs import notify_user

DOMAINS = (("both", "ref and mov"), ("ref", "ref only"), ("mov", "mov only"))
FORMAT_NAMES = {"cryosparc_cs": "CryoSPARC .cs", "relion_star": "RELION .star"}


class ExportPage(QWidget):
    command_issued = Signal(str)
    exported = Signal(object)  # the export report

    def __init__(self, view_page) -> None:
        super().__init__()
        self.view = view_page
        self.preview: dict = {}
        self.groups: dict[str, str | None] = {"ref": None, "mov": None}  # explicit .csg choices
        self._groups_run: Path | None = None  # the run the .csg choices belong to
        self._output_edited = False
        self._task = None
        self._preview_runner = LatestWinsRunner(self._preview_ready, self._preview_failed)
        self.destroyed.connect(self._preview_runner.discard)

        intro = QLabel(
            "Writes subsets of the run's source files: the selected particles' rows, unchanged. "
            "Nothing is reselected or rewritten, and an existing export folder is never overwritten."
        )
        intro.setWordWrap(True)
        self.run_label = QLabel()
        self.run_label.setObjectName("export_run")
        self.selection_combo = QComboBox()
        self.selection_combo.setObjectName("export_selection")
        self.domain_combo = QComboBox()
        for key, label in DOMAINS:
            self.domain_combo.addItem(label, key)

        sources = QGroupBox("Source files")
        sources_form = QFormLayout(sources)
        self.source_labels: dict[str, QLabel] = {}
        self.group_labels: dict[str, QLabel] = {}
        self.group_buttons: dict[str, QPushButton] = {}
        for domain in ("ref", "mov"):
            source, group = QLabel(), QLabel()
            for label in (source, group):
                label.setWordWrap(True)
                label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            button = QPushButton("Choose .csg...")
            button.setToolTip("The source job's particle result group, e.g. PROJECT/J75/J75_particles.csg.")
            button.clicked.connect(lambda _=False, d=domain: self.choose_group(d))
            row = QHBoxLayout()
            row.addWidget(group, 1)
            row.addWidget(button)
            sources_form.addRow(f"{domain}:", source)
            sources_form.addRow("", row)
            self.source_labels[domain], self.group_labels[domain], self.group_buttons[domain] = source, group, button

        self.output_edit = QLineEdit()
        self.output_edit.setObjectName("export_output")
        self.output_edit.textEdited.connect(lambda _: setattr(self, "_output_edited", True))
        browse = QPushButton("Browse...")
        browse.clicked.connect(self.choose_output)
        output_row = QHBoxLayout()
        output_row.addWidget(self.output_edit, 1)
        output_row.addWidget(browse)
        self.export_button = QPushButton("Export")
        self.export_button.setObjectName("export_button")
        self.command_label = QLabel()
        self.command_label.setObjectName("export_command")
        self.command_label.setWordWrap(True)
        self.command_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.command_label.setStyleSheet("font-family: monospace; color: #444;")
        self.result_view = QPlainTextEdit()
        self.result_view.setObjectName("export_result")
        self.result_view.setReadOnly(True)
        self.result_view.setPlaceholderText("After an export: the files written and what to do next.")
        self.result_view.setStyleSheet("font-family: monospace;")

        form = QFormLayout()
        form.addRow("Run", self.run_label)
        form.addRow("Selection", self.selection_combo)
        form.addRow("Domains", self.domain_combo)
        form.addRow("Output folder", output_row)
        layout = QVBoxLayout(self)
        layout.addWidget(intro)
        layout.addLayout(form)
        layout.addWidget(sources)
        buttons = QHBoxLayout()
        buttons.addWidget(self.export_button)
        buttons.addStretch(1)
        layout.addLayout(buttons)
        layout.addWidget(QLabel("Equivalent command:"))
        layout.addWidget(self.command_label)
        layout.addWidget(QLabel("Result and next steps:"))
        layout.addWidget(self.result_view, 1)

        self.selection_combo.currentIndexChanged.connect(lambda _: self._selection_changed())
        self.domain_combo.currentIndexChanged.connect(lambda _: self._update())
        self.output_edit.textChanged.connect(lambda _: self._update())
        self.export_button.clicked.connect(self.export)
        self.view.landscape_loaded.connect(lambda _: self.on_shown())
        self._update()

    # --- state -----------------------------------------------------------------

    @property
    def run_dir(self) -> Path | None:
        return None if self.view.landscape is None else self.view.landscape.run_dir

    @property
    def domain(self) -> str:
        return self.domain_combo.currentData()

    @property
    def selection_id(self) -> str | None:
        return self.selection_combo.currentData()

    def on_shown(self) -> None:
        """Refresh selections and source information for the open run."""

        run = self.run_dir
        if run != self._groups_run:  # .csg and output-folder choices belong to one run
            self.groups = {"ref": None, "mov": None}
            self._groups_run = run
            self._output_edited = False
            self._preview_ready({})
        self.refresh_selections()
        if run is None:
            self.run_label.setText("Open a run first (File > Open run...).")
            self._preview_runner.discard()
            self._preview_ready({})
            return
        self.run_label.setText(str(run))
        self._preview_runner.submit(lambda: api.export_preview(run))

    def refresh_selections(self, *_args) -> None:
        current = self.selection_id
        self.selection_combo.blockSignals(True)
        self.selection_combo.clear()
        if self.view.landscape is not None:
            for name, count in self.view.landscape.saved_selections():
                self.selection_combo.addItem(name + (f"  ({count:,} particles)" if count is not None else ""), name)
        index = self.selection_combo.findData(current)
        self.selection_combo.setCurrentIndex(index if index >= 0 else self.selection_combo.count() - 1)
        self.selection_combo.blockSignals(False)
        self._selection_changed()

    def _preview_ready(self, preview: dict) -> None:
        self.preview = preview or {}
        self._show_sources()
        self._update()

    def _preview_failed(self, exc: BaseException) -> None:
        self._preview_ready({})
        show_error(self, exc, title="cryoROLE: cannot preview export sources")

    def _show_sources(self) -> None:
        for domain in ("ref", "mov"):
            entry = self.preview.get(domain) or {}
            path = entry.get("path")
            kind = FORMAT_NAMES.get(entry.get("source_type"), entry.get("source_type") or "unknown format")
            self.source_labels[domain].setText(f"{path} ({kind})" if path else "not recorded in this run")
            is_cs = entry.get("source_type") == "cryosparc_cs"
            self.group_buttons[domain].setVisible(is_cs)
            if not is_cs:
                self.group_labels[domain].setText("")
            elif self.groups[domain]:
                self.group_labels[domain].setText(f"Result group: {self.groups[domain]} (chosen)")
            elif entry.get("cryosparc_group") and entry.get("cryosparc_group_usable") is not False:
                self.group_labels[domain].setText(f"Result group: {entry['cryosparc_group']} (found next to the source)")
            else:
                if entry.get("cryosparc_group_resolved_by") == "none":
                    note = (f"No result group (.csg) next to {Path(path).name}: the subset cannot go back through "
                            "Import Result Group. Choose the source job's particle .csg (e.g. J75/J75_particles.csg).")
                else:
                    note = entry.get("cryosparc_group_note") or "No result group (.csg) found."
                self.group_labels[domain].setText(f"<span style='color:#a15c00'>{escape(note, quote=False)}</span>")

    def _selection_changed(self) -> None:
        if not self._output_edited:
            run, name = self.run_dir, self.selection_id
            self.output_edit.setText(str(run / "exports" / name) if run and name else "")
        self._update()

    def _update(self) -> None:
        ready = self.run_dir is not None and self.selection_id is not None and bool(self.output_edit.text())
        self.export_button.setEnabled(ready and self._task is None)
        self.command_label.setText(self.cli_command() if ready else "")

    def choose_group(self, domain: str) -> None:
        start = str(Path((self.preview.get(domain) or {}).get("path") or Path.cwd()).parent)
        path, _ = QFileDialog.getOpenFileName(self, f"CryoSPARC result group for {domain}", start,
                                              "CryoSPARC result group (*.csg)")
        if path:
            self.set_group(domain, path)

    def set_group(self, domain: str, path: str | None) -> None:
        self.groups[domain] = path
        self._show_sources()
        self._update()

    def choose_output(self) -> None:
        start = str(self.run_dir / "exports") if self.run_dir else str(Path.cwd())
        folder = QFileDialog.getExistingDirectory(self, "Parent folder for the export", start)
        if folder and self.selection_id:
            self._output_edited = True
            self.output_edit.setText(str(Path(folder) / self.selection_id))

    # --- export ------------------------------------------------------------------

    def _options(self) -> dict:
        options = {"selection_id": self.selection_id, "domain": self.domain,
                   "output_dir": Path(self.output_edit.text()).expanduser()}
        for domain in ("ref", "mov"):
            if self.groups[domain] and self.domain in (domain, "both"):
                options[f"cryosparc_group_{domain}"] = self.groups[domain]
        return options

    def cli_command(self) -> str:
        options = self._options()
        parts = ["cryorole", "export", "--run-dir", str(self.run_dir), "--selection-id", options["selection_id"],
                 "--domain", options["domain"], "--output-dir", str(options["output_dir"])]
        for domain in ("ref", "mov"):
            if f"cryosparc_group_{domain}" in options:
                parts += [f"--cryosparc-group-{domain}", str(options[f"cryosparc_group_{domain}"])]
        return shlex.join(parts)

    def export(self) -> None:
        if self.run_dir is None or self.selection_id is None:
            return
        options = self._options()
        if options["output_dir"].exists():
            show_error(self, CryoroleError(
                "output_exists",
                f"The export folder {options['output_dir']} already exists.",
                "Choose another output folder; the GUI never overwrites an export.",
            ))
            return
        run, command = self.run_dir, self.cli_command()
        self.result_view.setPlainText("Exporting (verifying the source files, then writing subsets)...")
        self._task = run_in_background(
            lambda: api.export(run, **options),
            lambda report: self._exported(report, command),
            self._failed,
        )
        self._update()

    def _exported(self, report: dict, command: str) -> None:
        self._task = None
        lines = [f"Exported {report.get('selected_count', 0):,} particles to {report.get('output_dir')}", ""]
        for domain in ("ref", "mov"):
            path = (report.get("outputs") or {}).get(f"{domain}_metadata")
            if path:
                lines.append(f"{domain}: {path}")
        for domain, steps in (report.get("next_steps") or {}).items():
            if steps:
                lines += ["", f"[{domain}] next steps:"] + list(steps)
        if report.get("warnings"):
            lines += ["", "Warnings:"] + [f"  {warning}" for warning in report["warnings"]]
        self.result_view.setPlainText("\n".join(lines))
        notify_user(f"Export of {report.get('selection_id')}: {report.get('output_dir')}")
        self.command_issued.emit(command)
        self.exported.emit(report)
        self._output_edited = False
        self._update()

    def _failed(self, exc: BaseException) -> None:
        self._task = None
        self.result_view.setPlainText("")
        self._update()
        show_error(self, exc, title="cryoROLE: export failed")
