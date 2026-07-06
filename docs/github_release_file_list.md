# GitHub Release File List

This checklist is for migrating the cleaned cryoROLE 2.0 release candidate from
the staging workspace:

```text
C:\Users\xinzh\OneDrive\Desktop\work2022\codex
```

to the real GitHub repository:

```text
C:\Users\xinzh\OneDrive\Documents\GitHub\cryorole
```

Do not use this file as a request to copy files automatically. Review the list
before migration.

## Include

Top-level project files:

```text
.gitignore
AGENTS.md
LICENSE
README.md
environment.yml
pyproject.toml
```

Package source:

```text
cryorole/
```

Tests:

```text
tests/
```

Include only source test files and fixtures that are intentionally tracked. Do
not include `__pycache__` or generated artifacts.

Documentation:

```text
docs/AGENTS.md
docs/architecture.md
docs/production_run_plan.md
docs/roadmap.md
docs/installation.md
docs/quick_start.md
docs/output_files.md
docs/migration_from_0x.md
docs/relion_workflow.md
docs/cryosparc_workflow.md
docs/github_release_file_list.md
```

## Exclude

Local or generated directories:

```text
.codex_tmp/
.idea/
.pytest_cache/
__pycache__/
cryorole.egg-info/
cryorole_outputs/
selection_out/
alignments/
New folder1/
old_files/
```

Historical archives and local bundles:

```text
*.zip
*.rar
*.7z
```

Backup and internal prompt files:

```text
AGENTS.md.bak.euler_convention_20260608
docs/architecture.md.bak.euler_convention_20260608
docs/codex_prompt_canonicalize_fit_fraction_cleanup.md
docs/codex_prompt_production_run_progress.md
```

Large local data and manuscript files:

```text
data/tomo_subvol_job122.star
docs/manusicript_cryorole_version01.docx
*.mrc
*.map
*.cs
*.star
```

Root-level development helpers that should not be copied as public entry
points:

```text
main.py
star_align_exclude.py
cryorole2_landscape_projection.py
cryorole2_plot_1d_updated.py
```

The two `cryorole2_*` plotting scripts may be reconsidered later as
`tools/` or `examples/scripts/` helpers, but they should not be copied into the
initial clean release root.

## Legacy GitHub Folder

Do not copy `Github/` wholesale.

Use it only as a source for:

```text
Github/LICENSE -> LICENSE
Github/README.md -> background text already summarized in new docs
```

The old commands are documented in `docs/migration_from_0x.md`:

```text
orientation_analysis -> cryorole run
landscape_projection -> cryorole visualize
point_select -> cryorole select
particle_backtrack -> cryorole export
```

## Pre-Migration Checks

Run these from the staging workspace:

```bash
python -m pip install -e ".[test]" --dry-run
python -m cryorole.cli.main --help
python -m cryorole.cli.main run --help
python -m pytest tests\core tests\normalize tests\match tests\io -q
python -m pytest tests\docs -q
```

If time allows, run broader CLI/export/visualization tests before copying into
the real GitHub repository.

## Version Note

The package version remains:

```text
2.0.0a1
```

Public wording may describe this as a 2.0 beta-preview or release-candidate
workflow, but the package metadata should remain `2.0.0a1` until the project is
ready to cut a true `2.0.0b1` package version.
