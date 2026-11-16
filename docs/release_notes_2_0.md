# cryoROLE 2.0 pre-release notes

Current package version: `2.0.0a1`.

## What's new

This update adds a desktop GUI, makes SO(3) geodesic SLD the default, changes how CryoSPARC poses are read, adds tools for RELION inputs that do not match directly, and makes the command line easier to use.

- **Desktop GUI.** `cryorole gui` runs the whole workflow (Inputs, Parameters, Run, View, Select, Export) and opens existing run bundles. It uses the same services and writes the same artifacts as the command line. See [Desktop GUI](gui.md).
- **SLD uses SO(3) geodesic distance by default.** This is the intrinsic metric on rotations. `--sld-metric rotvec_euclidean` keeps the earlier, faster rotation-vector approximation for very large datasets; the metric used is recorded in every run.
- **Safer export.** An export is written to a staging folder and published only when it is complete, and CryoSPARC result groups are checked before anything is written. `--no-cryosparc-group` writes the `.cs` subset alone.

- **CryoSPARC pose convention (please re-run `.cs` analyses).** `alignments3D/pose` is now read the same way as pyem `csparc2star.py`, so a `.cs` input and the same particles converted to STAR give the same relative orientations. Runs made with earlier versions from `.cs` files have the same RO angles but different RO axes; re-run them. RELION `.star` inputs are unchanged. See the [CryoSPARC workflow](cryosparc_workflow.md).
- **`cryorole align` for RELION.** Pairs particles when default matching fails, for example after signal subtraction or re-extraction: explicit key pairs (`--key-pair`), exact matching of recentred re-extractions (`--coordinate-match recentered-exact`), and chains of jobs. `--fix-subtract-coordinates` writes a verified copy of a recentred subtraction with corrected coordinates. `preflight` suggests the right command. See the [RELION workflow](relion_workflow.md#preparing-aligned-star-files).
- **Input sanity checks.** `preflight` blocks when `--ref` and `--mov` are the same file, and every run reports a relative-orientation angle summary with warnings for identical or nearly identical poses.
- **Shorter commands and clearer errors.** `--run-dir`, `--canonical-id` and, for `export`, `--selection-id` can be omitted when there is only one candidate. Errors end with a `what to do:` line. `cryorole --version` works, and `COMMAND --help` ends with examples. See [Common behaviour](cli_reference.md#common-behaviour).
- **Python API.** `from cryorole import api` gives typed functions for every step, with progress reporting and cancellation. Actionable service errors are raised as `api.CryoroleError` with a stable `code`; cancellation remains `api.CancelledError`. See [Python API](python_api.md).
- **Reproducible random selections.** A random selection without `--seed` now records the seed it used.

A demo dataset (the INO80–hexasome landscape from the paper) is described under [Demo](../README.md#2-get-the-demo).

