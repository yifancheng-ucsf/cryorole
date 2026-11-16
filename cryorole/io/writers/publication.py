"""Cooperative directory publication for complete metadata artifact groups."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any
import warnings

from cryorole.provenance.output_paths import protect_input_directory, protect_input_file


class ArtifactFiles:
    """Prepare named files, then replace their destinations as a recoverable group.

    Each file uses a stage on its destination filesystem. Cooperative writers
    take every destination lock before staging. Caught failures restore the
    previous group; failed recovery leaves the backups and locks in place.
    Readers may observe individual renames, and power-loss recovery is manual.
    """

    def __init__(self, outputs, *, overwrite=False, inputs=()):
        self.outputs = {name: Path(path).expanduser().absolute() for name, path in outputs.items()}
        self.outputs = {name: path.parent.resolve() / path.name for name, path in self.outputs.items()}
        self.inputs, self.overwrite = tuple(inputs), overwrite
        self.paths, self.backups, self.stages, self.locks = {}, {}, {}, []
        self.committed = self.recovery_required = False
        self._validate()

    def _validate(self):
        outputs = list(self.outputs.values())
        for i, path in enumerate(outputs):
            for source in self.inputs:
                protect_input_file(path, source)
            for other in outputs[:i]:
                if path.resolve() == other.resolve() or (path.exists() and other.exists() and path.samefile(other)):
                    raise ValueError('Artifact output paths must be distinct files')
            if path.exists() or path.is_symlink():
                if not path.is_file():
                    raise ValueError(f'Output path is not a regular file: {path}')
                if not self.overwrite:
                    raise FileExistsError(f'Output path already exists: {path}')

    def __enter__(self):
        try:
            for path in sorted(self.outputs.values()):
                path.parent.mkdir(parents=True, exist_ok=True)
                lock = path.parent / f'.{path.name}.file-publish'
                try:
                    lock.mkdir()
                except FileExistsError as error:
                    raise FileExistsError(f'Another writer is publishing to {path}') from error
                self.locks.append(lock)
            self._validate()
            for name, path in self.outputs.items():
                stage = Path(tempfile.mkdtemp(prefix=f'.{path.name}.stage-', dir=path.parent))
                self.stages[name] = stage
                self.paths[name], self.backups[name] = stage / 'prepared', stage / 'previous'
            return self
        except BaseException:
            self._cleanup()
            raise

    def _commit(self):
        self._validate()
        for path in self.paths.values():
            if not path.is_file():
                raise ValueError(f'Artifact was not prepared: {path}')
        try:
            for name, destination in self.outputs.items():
                if destination.exists() or destination.is_symlink():
                    os.replace(destination, self.backups[name])
                os.replace(self.paths[name], destination)
        except BaseException:
            try:
                for name, destination in reversed(list(self.outputs.items())):
                    if not self.paths[name].exists() and (destination.exists() or destination.is_symlink()):
                        destination.unlink()
                    if self.backups[name].exists() or self.backups[name].is_symlink():
                        os.replace(self.backups[name], destination)
            except BaseException as error:
                self.recovery_required = True
                locations = ', '.join(map(str, self.stages.values()))
                raise RuntimeError(f'File publication rollback failed; recover preserved files in {locations}') from error
            raise
        self.committed = True

    def _cleanup(self):
        for stage in self.stages.values():
            shutil.rmtree(stage)
        for lock in reversed(self.locks):
            lock.rmdir()

    def __exit__(self, kind, value, traceback):
        try:
            if kind is None:
                self._commit()
        finally:
            if self.committed:
                try:
                    self._cleanup()
                except OSError as error:
                    warnings.warn(f'Artifacts published; leftover stage/lock needs cleanup: {error}', RuntimeWarning)
            elif not self.recovery_required:
                self._cleanup()
        return False


class ArtifactDirectory:
    """Stage a complete artifact group; preserve existing files on caught failures.

    Unrelated files in an explicitly overwritten directory are retained. A failed
    rollback or uncatchable process kill leaves backups for manual recovery.
    """

    def __init__(self, output_dir, *, overwrite=False, inputs=(), allow_empty=False, label=None,
                 omit_existing=()):
        self.output = Path(output_dir).expanduser().resolve()
        self.input_paths = tuple(inputs)
        protect_input_directory(output_dir, self.input_paths)
        self.overwrite, self.allow_empty = overwrite, allow_empty
        self.label = label
        # Managed optional outputs absent from the new group must not survive.
        self.omit_existing = frozenset(omit_existing)
        self.lock = self.output.parent / f'.{self.output.name}.publication-lock'
        self.scratch = None
        self.committed = False
        self.recovery_required = False

    def _validate(self):
        label = f'{self.label} output' if self.label else 'Output'
        if self.output.exists() and not self.output.is_dir():
            raise ValueError(f'{label} path exists and is not a directory: {self.output}')
        if self.output.exists() and not self.overwrite and not (self.allow_empty and not any(self.output.iterdir())):
            raise FileExistsError(f'{label} directory already exists: {self.output}. Use --overwrite or another output directory.')
        protect_input_directory(self.output, self.input_paths)

    def __enter__(self):
        self._validate()
        self.output.parent.mkdir(parents=True, exist_ok=True)
        try:
            self.lock.mkdir()
        except FileExistsError as error:
            raise FileExistsError(f'Another writer is publishing to {self.output}') from error
        try:
            self._validate()
            self.scratch = Path(tempfile.mkdtemp(prefix=f'.{self.output.name}.stage-', dir=self.output.parent))
            self.path, self.backup = self.scratch / 'artifacts', self.scratch / 'previous'
            self.path.mkdir()
            return self
        except BaseException:
            self._cleanup()
            raise

    def protect(self, inputs):
        self.input_paths = (*self.input_paths, *inputs)
        protect_input_directory(self.output, self.input_paths)

    def publicize(self, value: Any) -> Any:
        if isinstance(value, Path):
            return Path(str(value).replace(str(self.path), str(self.output)))
        if isinstance(value, str):
            return value.replace(str(self.path), str(self.output))
        if isinstance(value, dict):
            return {key: self.publicize(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            converted = [self.publicize(item) for item in value]
            return tuple(converted) if isinstance(value, tuple) else converted
        return value

    def _keep_existing(self):
        for source in self.output.iterdir():
            if source.name in self.omit_existing:
                continue
            target = self.path / source.name
            if target.exists() or target.is_symlink():
                continue
            if source.is_symlink():
                target.symlink_to(os.readlink(source), target_is_directory=source.is_dir())
            elif source.is_dir():
                shutil.copytree(source, target, symlinks=True)
            else:
                shutil.copy2(source, target)

    def _commit(self):
        self._validate()
        for path in self.path.rglob('*.json'):
            payload = json.loads(path.read_text(encoding='utf-8'))
            path.write_text(json.dumps(self.publicize(payload), indent=2, sort_keys=True) + '\n', encoding='utf-8')
        if self.output.exists():
            self._keep_existing()
        try:
            if self.output.exists():
                os.replace(self.output, self.backup)
            os.replace(self.path, self.output)
        except BaseException:
            try:
                if self.backup.exists():
                    if self.output.exists():
                        shutil.rmtree(self.output)
                    os.replace(self.backup, self.output)
                elif not self.path.exists() and self.output.exists():
                    shutil.rmtree(self.output)
            except BaseException as error:
                self.recovery_required = True
                raise RuntimeError(f'Publication rollback failed; recover preserved files in {self.scratch}') from error
            raise
        self.committed = True

    def _cleanup(self):
        if self.scratch is not None and self.scratch.exists():
            shutil.rmtree(self.scratch)
        self.lock.rmdir()

    def __exit__(self, kind, value, traceback):
        try:
            if kind is None:
                self._commit()
        finally:
            if self.committed:
                try:
                    self._cleanup()
                except OSError as error:
                    warnings.warn(f'Artifacts published; leftover stage/lock needs cleanup: {error}', RuntimeWarning)
            elif not self.recovery_required:
                self._cleanup()
        return False
