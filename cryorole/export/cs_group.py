"""CryoSPARC result-group (.csg) export, so a .cs subset can go back through Import Result Group.

A CryoSPARC job describes each output group in a small YAML ``.csg`` file: one
entry per slot (``alignments3D``, ``ctf``, ``blob``, ...) naming the ``.cs``
file that holds it, with ``'>'`` meaning "relative to the .csg folder". A
refinement's particle group usually takes ``alignments3D`` and ``ctf`` from the
job's particle file and ``blob`` from its ``*_passthrough_particles.cs``.

The export mirrors the source job's own group exactly. The run's source ``.cs``
is subset by the recorded source rows (as before); every other file the group
names is subset by ``uid`` in the same order, and every selected ``uid`` must be
present exactly once. Nothing is reselected, reordered or rewritten.

Verified 2026-10-05 on CryoSPARC v5.0.6: the subset group of a refinement job
imported and reconstructed. Import Result Group only accepts a ``.csg`` located
inside the project directory; the export report says where to copy the files.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shlex
from typing import Any, Sequence

import numpy as np

GROUP_SUFFIX = ".csg"
_KEY = re.compile(r"^(?P<indent> *)(?P<key>[A-Za-z0-9_./-]+):(?: (?P<value>.*))?$")


@dataclass(frozen=True)
class CsGroup:
    """One parsed ``.csg`` file."""

    path: Path
    created: str | None
    group: dict[str, str]
    results: dict[str, dict[str, Any]]  # slot -> {"metafile", "num_items", "type"}
    version: str | None

    def metafile_name(self, slot: str) -> str:
        return _metafile_name(self.results[slot]["metafile"], self.path)

    @property
    def metafile_names(self) -> list[str]:
        """Distinct referenced file names, in slot order."""

        names: list[str] = []
        for slot in self.results:
            name = self.metafile_name(slot)
            if name not in names:
                names.append(name)
        return names


@dataclass(frozen=True)
class CsGroupExport:
    """What the group export wrote for one domain."""

    source_group: str
    resolved_by: str  # "explicit" or "only_candidate"
    group_path: str
    files: dict[str, str]  # source metafile name -> written subset path
    slots: dict[str, str]  # slot -> written file name
    selected_count: int
    version: str | None
    project_dir: str | None
    import_steps: list[str] = field(default_factory=list)


# --- reading and writing .csg ----------------------------------------------------------


def read_cs_group(path: str | Path) -> CsGroup:
    """Parse a CryoSPARC ``.csg`` file (the small, machine-written YAML subset CryoSPARC uses)."""

    group_path = Path(path)
    if not group_path.is_file():
        raise ValueError(f"CryoSPARC group file does not exist: {group_path}")
    tree = _parse_yaml_mapping(group_path.read_text(encoding="utf-8"), group_path)
    group = tree.get("group")
    results = tree.get("results")
    if not isinstance(group, dict) or not isinstance(results, dict) or not results:
        raise ValueError(f"{group_path} is not a CryoSPARC result group (needs 'group' and 'results').")
    for slot, entry in results.items():
        if not isinstance(entry, dict) or "metafile" not in entry or "num_items" not in entry:
            raise ValueError(f"{group_path}: slot {slot!r} needs 'metafile' and 'num_items'.")
        try:
            entry["num_items"] = int(entry["num_items"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{group_path}: slot {slot!r} has a non-integer num_items.") from exc
        _metafile_name(entry["metafile"], group_path)  # validates the reference
    return CsGroup(
        path=group_path,
        created=_optional_str(tree.get("created")),
        group={str(key): str(value) for key, value in group.items()},
        results=results,
        version=_optional_str(tree.get("version")),
    )


def format_cs_group(
    *,
    group: dict[str, str],
    results: dict[str, dict[str, Any]],
    version: str | None,
    created: str | None = None,
) -> str:
    """Write a ``.csg`` in CryoSPARC's own layout (sorted keys, two-space indent)."""

    created = created or datetime.now(timezone.utc).isoformat(sep=" ")
    lines = [f"created: {created}", "group:"]
    for key in sorted(group):
        lines.append(_folded(f"  {key}: ", _scalar(group[key]), indent=2))
    lines.append("results:")
    for slot in sorted(results):
        lines.append(f"  {slot}:")
        for key in sorted(results[slot]):
            value = results[slot][key]
            text = _quote(value) if key == "metafile" else _scalar(value)
            lines.append(_folded(f"    {key}: ", text, indent=4))
    if version is not None:
        lines.append(f"version: {_scalar(version)}")
    return "\n".join(lines) + "\n"


_FOLD_WIDTH = 80  # the YAML writer CryoSPARC uses (PyYAML) folds plain scalars past this column


def _folded(prefix: str, text: str, *, indent: int) -> str:
    """``prefix + text``, a long plain scalar folded as CryoSPARC's writer does: at the first
    single space after column 80, continuing two spaces deeper than the key."""

    if text[:1] in ("'", '"') or len(prefix) + len(text) <= _FOLD_WIDTH:
        return prefix + text
    out, column = prefix, len(prefix)
    words = text.split(" ")
    for number, word in enumerate(words):
        if number:
            single = word != "" and words[number - 1] != ""
            if single and column > _FOLD_WIDTH:
                out += "\n" + " " * (indent + 2)
                column = indent + 2
            else:
                out += " "
                column += 1
        out += word
        column += len(word)
    return out


# --- finding the source job's group ----------------------------------------------------


def group_candidates(source_path: str | Path) -> list[Path]:
    """``.csg`` files next to a source ``.cs`` that describe a particle group naming that file."""

    source = Path(source_path)
    candidates = []
    for path in sorted(source.parent.glob(f"*{GROUP_SUFFIX}")):
        try:
            group = read_cs_group(path)
        except (ValueError, OSError):
            continue
        if group.group.get("type") == "particle" and source.name in group.metafile_names:
            candidates.append(path)
    return candidates


def resolve_group(source_path: str | Path, explicit: str | Path | None) -> tuple[CsGroup | None, str, list[str]]:
    """The group to mirror: explicit, else the only candidate next to the source; else None.

    Returns (group, resolved_by, candidates). Several candidates are an error: never guess.
    """

    if explicit is not None:
        return read_cs_group(explicit), "explicit", [str(explicit)]
    candidates = group_candidates(source_path)
    if len(candidates) > 1:
        raise ValueError(
            f"Several CryoSPARC particle groups name {Path(source_path).name}: "
            + ", ".join(str(path) for path in candidates)
            + ". Choose one with --cryosparc-group-ref/--cryosparc-group-mov."
        )
    if not candidates:
        return None, "none", []
    return read_cs_group(candidates[0]), "only_candidate", [str(candidates[0])]


# --- checking and writing the subset group ---------------------------------------------


def check_cs_group(group: CsGroup, source_path: str | Path, source_row_count: int | None) -> None:
    """Raise ValueError if ``group`` cannot be mirrored for ``source_path``; writes nothing.

    The same checks the subset writer makes, run before an export writes its first file: a
    particle group that names the source, with the source's row count, and every other file
    it names present next to it, with a uid field and the row count the group states. (Missing
    or duplicate uids need the full files and are found while writing.)
    """

    source = Path(source_path)
    if group.group.get("type") != "particle":
        raise ValueError(f"{group.path} is a {group.group.get('type')!r} group; export needs a particle group.")
    names = group.metafile_names
    if source.name not in names:
        raise ValueError(
            f"{group.path} does not name the run's source file {source.name} "
            f"(it names {', '.join(names)}); it belongs to another job or output."
        )
    if source_row_count is None:
        source_row_count = len(_cs_header(source))
    for slot in group.results:
        if group.metafile_name(slot) == source.name and group.results[slot]["num_items"] != source_row_count:
            raise ValueError(
                f"{group.path} describes {group.results[slot]['num_items']:,} particles in {source.name}, "
                f"but the run's source has {source_row_count:,}; it is not the same file."
            )
    if "uid" not in (_cs_header(source).dtype.names or ()):
        raise ValueError(f"{source.name} has no uid field; other group files cannot be matched to it.")
    for name in names:
        if name == source.name:
            continue
        other = group.path.parent / name
        if not other.is_file():
            raise ValueError(f"CryoSPARC group file {name} is missing next to its .csg: {other}")
        array = _cs_header(other)
        if array.dtype.names is None or "uid" not in array.dtype.names:
            raise ValueError(f"{other} is not a CryoSPARC .cs file with a uid field.")
        expected = {group.results[slot]["num_items"] for slot in group.results if group.metafile_name(slot) == name}
        if expected and len(array) not in expected:
            raise ValueError(f"{other} has {len(array):,} rows; its .csg says {sorted(expected)}.")


def _cs_header(path: Path) -> np.ndarray:
    """A .cs file memory-mapped: its fields and row count without reading the rows."""

    try:
        return np.load(path, mmap_mode="r", allow_pickle=False)
    except (OSError, ValueError) as exc:
        raise ValueError(f"{path} could not be read as a CryoSPARC .cs file: {exc}") from exc


def write_cs_group_subset(
    group: CsGroup,
    *,
    resolved_by: str,
    source_path: str | Path,
    source_row_count: int,
    main_subset_path: str | Path,
    domain: str,
    domain_dir: str | Path,
    selection_id: str,
    run_id: str | None,
    overwrite: bool = False,
) -> CsGroupExport:
    """Subset every file of ``group`` like the already written main subset; write the new ``.csg``."""

    source = Path(source_path)
    out_dir = Path(domain_dir)
    check_cs_group(group, source, source_row_count)
    names = group.metafile_names
    main_subset = np.load(main_subset_path, allow_pickle=False)
    if "uid" not in (main_subset.dtype.names or ()):
        raise ValueError(f"{source.name} has no uid field; other group files cannot be matched to it.")
    uids = np.asarray(main_subset["uid"])

    written = {source.name: Path(main_subset_path)}
    for name in names:
        if name == source.name:
            continue
        other_path = group.path.parent / name
        expected = {group.results[slot]["num_items"] for slot in group.results if group.metafile_name(slot) == name}
        target = out_dir / f"selected_{domain}_{name}"
        _write_uid_subset(other_path, target, uids, expected_rows=expected, overwrite=overwrite)
        written[name] = target

    slots = {slot: written[group.metafile_name(slot)].name for slot in group.results}
    results = {
        slot: {
            "metafile": ">" + slots[slot],
            "num_items": int(len(uids)),
            "type": group.results[slot].get("type", f"particle.{slot}"),
        }
        for slot in group.results
    }
    job = group.path.parent.name
    new_group = dict(group.group)
    new_group["title"] = f"cryoROLE {selection_id}"
    new_group["description"] = (
        f"cryoROLE selection {selection_id} ({len(uids):,} particles) of {job} {group.group.get('name', '')}"
        + (f", run {run_id}" if run_id else "")
    ).strip()
    csg_path = out_dir / f"selected_{domain}{GROUP_SUFFIX}"
    if csg_path.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {csg_path}")
    csg_path.write_text(format_cs_group(group=new_group, results=results, version=group.version), encoding="utf-8")

    project_dir = _project_dir(group.path)
    files = [str(written[name]) for name in names] + [str(csg_path)]
    return CsGroupExport(
        source_group=str(group.path),
        resolved_by=resolved_by,
        group_path=str(csg_path),
        files={name: str(path) for name, path in written.items()},
        slots=slots,
        selected_count=int(len(uids)),
        version=group.version,
        project_dir=str(project_dir) if project_dir else None,
        import_steps=import_steps(files, csg_path.name, project_dir, folder=f"cryorole_{job}_{selection_id}_{domain}"),
    )


def import_steps(files: Sequence[str], csg_name: str, project_dir: Path | None, *, folder: str) -> list[str]:
    """Plain-text steps to bring the written group into CryoSPARC with Import Result Group."""

    if project_dir is None:
        return [
            "Copy these files together into a folder inside your CryoSPARC project directory "
            "(Import Result Group only accepts files inside the project):",
            "  " + " ".join(shlex.quote(path) for path in files),
            f"In CryoSPARC, run Import Result Group with the path of {csg_name} in that folder.",
        ]
    destination = project_dir / "exports" / "groups" / folder
    prefix = _run_as_owner_prefix(project_dir)
    steps = [
        "Import Result Group only accepts files inside the project directory. Copy the files there"
        + (" as the project's owner:" if prefix else ":"),
        f"  {prefix}mkdir -p {shlex.quote(str(destination))}",
        f"  {prefix}cp " + " ".join(shlex.quote(path) for path in files) + f" {shlex.quote(str(destination))}/",
        "In CryoSPARC, run Import Result Group with this path:",
        f"  {destination / csg_name}",
        "The import copies the files into its own job folder; the exports folder can be removed afterwards.",
    ]
    return steps


# --- helpers ---------------------------------------------------------------------------


def _write_uid_subset(source: Path, target: Path, uids: np.ndarray, *, expected_rows: set[int], overwrite: bool) -> None:
    if not source.is_file():
        raise ValueError(f"CryoSPARC group file {source.name} is missing next to its .csg: {source}")
    if target.exists() and not overwrite:
        raise FileExistsError(f"Output path already exists: {target}")
    array = np.load(source, allow_pickle=False)
    if array.dtype.names is None or "uid" not in array.dtype.names:
        raise ValueError(f"{source} is not a CryoSPARC .cs file with a uid field.")
    if expected_rows and len(array) not in expected_rows:
        raise ValueError(f"{source} has {len(array):,} rows; its .csg says {sorted(expected_rows)}.")
    order = np.argsort(array["uid"], kind="stable")
    sorted_uids = array["uid"][order]
    if np.any(sorted_uids[1:] == sorted_uids[:-1]):
        raise ValueError(f"{source} contains duplicate uids; rows cannot be matched by uid.")
    positions = np.searchsorted(sorted_uids, uids)
    found = (positions < len(sorted_uids)) & (sorted_uids[np.minimum(positions, len(sorted_uids) - 1)] == uids)
    if not np.all(found):
        missing = uids[~found][:5].tolist()
        raise ValueError(
            f"{int((~found).sum()):,} selected uids are missing from {source.name} (e.g. {missing}); "
            "the group's files do not match the run's source."
        )
    subset = array[order[positions]]
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("wb") as handle:
        np.save(handle, subset)


def _project_dir(group_path: Path) -> Path | None:
    """The CryoSPARC project folder when the .csg sits in one of its job folders."""

    job_dir = group_path.parent
    project = job_dir.parent
    if not (re.fullmatch(r"J\d+", job_dir.name) and (project / "project.json").is_file()):
        return None
    # Use the path CryoSPARC registered for the project: Import Result Group compares against it,
    # so a symlinked path to the same folder would be rejected.
    try:
        registered = json.loads((project / "project.json").read_text(encoding="utf-8")).get("project_dir")
    except (OSError, ValueError):
        registered = None
    if registered and Path(registered).exists() and Path(registered).resolve() == project.resolve():
        return Path(registered)
    return project.resolve()


def _run_as_owner_prefix(project_dir: Path) -> str:
    try:
        owner = project_dir.owner()
    except (KeyError, NotImplementedError, OSError):
        return ""
    if os.access(project_dir, os.W_OK):
        return ""
    return f"sudo -u {shlex.quote(owner)} "


def _metafile_name(reference: Any, group_path: Path) -> str:
    text = str(reference)
    if not text.startswith(">") or "/" in text[1:] or not text[1:]:
        raise ValueError(
            f"{group_path}: metafile {text!r} is not a file next to the .csg ('>NAME'); not supported."
        )
    return text[1:]


def _optional_str(value: Any) -> str | None:
    return None if value is None else str(value)


def _parse_yaml_mapping(text: str, path: Path) -> dict[str, Any]:
    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]
    # A scalar may continue on the following lines: CryoSPARC's YAML writer folds long values
    # (e.g. Import Particles descriptions) onto lines indented deeper than their key, and a
    # quoted value continues until its closing quote. YAML joins the pieces with single spaces.
    pending: list[Any] | None = None  # [parent, key, key indent, text, first line number]

    def finish() -> None:
        if pending is not None:
            parent_, key_, _indent, value_text, first = pending
            parent_[key_] = _parse_scalar(value_text, path, first)

    for number, raw in enumerate(text.splitlines(), start=1):
        if not raw.strip():
            continue
        line_indent = len(raw) - len(raw.lstrip(" "))
        match = _KEY.match(raw.rstrip())
        if pending is not None and (
            _open_quote(pending[3]) or (match is None and line_indent > pending[2])
        ):
            pending[3] = f"{pending[3]} {raw.strip()}"
            continue
        finish()
        pending = None
        if raw.lstrip().startswith("#"):
            continue
        if match is None:
            raise ValueError(f"{path}:{number}: unsupported .csg line {raw.strip()!r}.")
        indent = len(match["indent"])
        while stack and indent <= stack[-1][0]:
            stack.pop()
        if not stack:
            raise ValueError(f"{path}:{number}: bad indentation.")
        parent = stack[-1][1]
        value = match["value"]
        if value is None or value == "":
            child: dict[str, Any] = {}
            parent[match["key"]] = child
            stack.append((indent, child))
        else:
            pending = [parent, match["key"], indent, value.strip(), number]
    finish()
    return root


def _open_quote(value: str) -> bool:
    """True while a quoted scalar has not reached its closing quote."""

    if value.startswith("'"):
        return not value[1:].replace("''", "").endswith("'")
    if value.startswith('"'):
        return len(value) < 2 or not value.endswith('"') or value.endswith('\\"')
    return False


def _parse_scalar(value: str, path: Path, number: int) -> Any:
    value = value.strip()
    if value.startswith("'"):
        if len(value) < 2 or not value.endswith("'"):
            raise ValueError(f"{path}:{number}: unterminated quoted value.")
        return value[1:-1].replace("''", "'")
    if value.startswith('"'):
        if len(value) < 2 or not value.endswith('"'):
            raise ValueError(f"{path}:{number}: unterminated quoted value.")
        return value[1:-1]
    if value.startswith(("[", "{", "- ", "|", ">")):
        raise ValueError(f"{path}:{number}: unsupported .csg value {value!r}.")
    if re.fullmatch(r"-?\d+", value):
        return int(value)
    return value


def _quote(value: Any) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _scalar(value: Any) -> str:
    """A YAML plain scalar where CryoSPARC writes one, quoted only when YAML needs it."""

    text = str(value)
    needs_quotes = (
        text == ""
        or text != text.strip()
        or ": " in text
        or " #" in text
        or text.endswith(":")
        or text[0] in "-?:,[]{}#&*!|>'\"%@`"
    )
    return _quote(text) if needs_quotes else text
