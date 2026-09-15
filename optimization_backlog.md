# Optimization Backlog

This document tracks optimization ideas, agreed plans, and completed work. It does not replace `AGENTS.md`, architecture documents, or formal implementation plans. Update the relevant authoritative documents and tests during implementation.

## Maintenance Rules

- Write persistent records in English by default unless the user explicitly requests Chinese.
- Organize items under Small Optimizations, Major Optimizations, and Completed. Scope and priority are independent; order does not imply priority.
- Assign unique, sequential `OPT-NNN` identifiers. Preserve them when moving items and never reuse them.
- Use statuses: **Needs discussion / Agreed, pending implementation / In progress / Completed**.
- Recording an idea does not authorize implementation. Mark the direction as agreed only after discussion, and explicitly retain unresolved decisions.
- Keep problems, proposed changes, and acceptance criteria aligned with the latest agreement. Record the rationale for material changes.
- After implementation and verification, move the entire item to Completed and record its completion date, changed files, and verification results.
- Major optimizations should additionally record open questions, affected areas, and possible approaches. Distinguish tentative approaches from agreed decisions.

## Small Optimizations

### OPT-002: Fix JavaScript newline escaping in offline interactive 3D HTML

**Status:** In progress; code and executable-script checks complete, Chrome acceptance pending.

**Recorded:** 2026-09-07

**Problem:** The offline `landscape_3d.html` viewer can appear blank because its inline JavaScript fails to parse. The Python HTML template interprets the newline escape inside the tooltip's JavaScript `.join("\n")` expression, emitting a literal line break inside a double-quoted JavaScript string. This causes `SyntaxError: Invalid or unexpected token` and prevents initialization, drawing, and interaction.

**Evidence:** Inspected `C:/Users/xinzh/OneDrive/Desktop/work2022/cryorole2_test/0907/sld_1.5/landscape_3d.html`. Its embedded JSON parses successfully and contains 50,000 points in each of the rotvec and Euler representations, with matching particle-key, SLD, and color array lengths. It has no external script dependencies. Parsing the actual inline script with Node.js reproduced the syntax error; correcting only the join-string newline in memory passed syntax validation. Corrected rendering has not yet been verified in Chrome. Neither the original HTML nor the generator was modified during diagnosis.

**Proposed changes:**

- Correct Python-to-JavaScript newline escaping in `_INTERACTIVE_3D_HTML` in `cryorole/visualize/service.py`, preserving readable multiline tooltips.
- Add regression coverage that parses the generated executable JavaScript, rather than checking only that an HTML file was written or contains expected text.
- Verify actual offline rendering and interaction in Chrome using representative data. If a corrected copy of the reported artifact is produced, preserve the original and label the copy clearly.

**Acceptance criteria:**

- Generated inline JavaScript parses successfully, including tooltip string handling.
- Opening the generated HTML locally in Chrome displays the point cloud and summary without script errors; representation switching, rotation, zoom, reset, and hover tooltips work.
- Verify rendering with a representative 50,000-point display sample, matching the reported artifact's scale.
- The viewer remains self-contained, offline, display-only, and unable to create a Selection. Scientific coordinates, SLD, sampling policy, and source artifacts remain unchanged.
- Update relevant documentation and tests; record script-validation and browser-validation results separately. Syntax validation alone is not evidence of successful rendering.

**Implementation record (2026-09-14):** Fixed template escaping, large-array max evaluation, canvas-relative hover coordinates, and pointer-cancellation handling in `cryorole/visualize/service.py`. `tests/visualize/test_offline_viewer.py` exercises generated JavaScript, tooltip escaping, rotation, zoom, reset, representation switching, and 150,000-point drawing with a canvas stub. The original 50,000-point HTML was preserved; a corrected copy is at `.codex_tmp/optimization_20260913/landscape_3d_fixed.html`. Chrome UI automation stopped because the tool could not confidently determine the local browser URL for policy enforcement. Actual Chrome rendering and interaction remain unverified; do not mark this item Completed on script tests alone.

## Major Optimizations

No items yet. Include open questions, affected areas, and possible approaches alongside the standard fields.

## Completed

### OPT-001: Clarify RO coordinate coincidence diagnostics and refine warning criteria

**Status:** Completed

**Recorded:** 2026-09-07

**Problem:** The `many_near_duplicate_ro_coordinates` warning can imply that many particles are duplicates. The diagnostic measures quantized RO coordinate coincidence, which cannot establish duplicate particle identity or images. The term `largest_cluster` can also suggest angular-distance clustering.

**Pre-update behavior:** Each RO analysis rotation-vector component is quantized with `round(coordinate / tolerance)`, using a default grid step of `1e-8 rad`. Rows with identical quantized coordinates form a group. The diagnostic counts all rows in groups of at least two and reports the largest group size. It warns when the participating row count is at least 10 and at least 1% of all rows. The trigger ignores the largest group size and does not check particle identity, source duplication, high SLD, or distance-floor activation.

**Discussion example:** 48,807 of 356,280 rows (approximately 13.70%) belong to coincident groups, with a largest group of 3. These dispersed groups of two or three triggered the previous warning. The participating count is neither a duplicate-particle count nor a count of particles to remove.

**Agreed direction:** Preserve routine coincidence statistics as informational diagnostics. Retain warnings for substantial RO coordinate concentration or supported evidence of duplicate particle information; do not remove warnings unconditionally.

**Proposed changes:**

- Preserve existing statistics and report the grid step, participating row count and percentage, largest group size, and interpretation limits.
- Distinguish dispersed small groups from large concentrations of nearly identical RO coordinates. Consider the largest group's absolute size and fraction when defining warning criteria.
- Diagnose duplicate particle information separately using explicitly defined identity keys or source provenance. State the evidence and its limitations; never infer duplicate particle identity or images from RO coincidence alone.
- Make warning wording identify coordinate concentration or particle-information duplication as appropriate. Use "group" for quantization-based statistics rather than implying SO(3) distance clustering.
- Keep diagnostics report-only; do not automatically remove particles or change scientific computations.

**Approved policy (2026-09-13):** Count rows in quantized groups of at least 10. Warn when their combined count is at least 100 or at least 1% of all landscape rows. These are explicit heuristic thresholds, not established scientific cutoffs. Particle-duplication evidence remains limited to existing input identity/matching diagnostics; existing hard failures are preserved. Row-aligned mode gains no implicit identity check. Broader SO(3) clustering or image-duplication detection is outside this item.

**Acceptance criteria:**

- Existing coincidence statistics remain numerically unchanged; RO, SLD, matching, selection, and export semantics are preserved.
- Routine dispersed small-group coincidence is informational and does not warn solely because the total participating fraction exceeds 1%.
- Substantial concentration triggers a warning under agreed explicit criteria. Supported particle-information duplication retains or receives an appropriate warning without weakening existing input/matching failures.
- Reports distinguish coordinate coincidence from particle identity duplication and state relevant thresholds or evidence. Coordinate coincidence alone does not imply a need to remove particles.
- Update relevant documentation and tests for preserved statistics, routine informational output, abnormal concentration warnings, threshold boundaries, and supported versus unsupported duplication claims.

**Example informational wording (subject to implementation review):**

> RO coordinate coincidences: 48,807/356,280 rows (13.70%) share quantized RO coordinates with another row (grid step: 1e-8 rad; largest group: 3). This measures coordinate coincidence, not duplicate particle identity. No particles were removed.

**Decision update (2026-09-07):** Replaced the initial blanket warning removal with differentiated informational and warning behavior. Abnormal concentration and supported duplication evidence must remain visible as warnings.

**Completion record (2026-09-14):** Implemented quantized-group diagnostics with the approved 10-row / 100-row / 1% heuristic, preserved legacy counts, and recorded the resolved policy and diagnostic in reports, summary, and manifest. Existing identity/matching failures remain unchanged. Changed files: `cryorole/core/density.py`, `cryorole/models/policies.py`, `cryorole/models/density_report.py`, `cryorole/workflows/run_report.py`, and `cryorole/workflows/run_service.py`. Verification: threshold boundaries, informational pairs/triples, configurable policy, array/DataFrame parity, unchanged scientific arrays, and persisted report/summary/manifest agreement in `tests/core/test_ro_coordinate_diagnostics.py` and `tests/run/test_optimization_workflow.py`. See the verification record below.

### OPT-003: Standardize user-facing color bar labels as SLD

**Status:** Completed

**Recorded:** 2026-09-13

**Problem:** Color bar labels such as `sld_raw` and `sld_display` expose implementation field names and can confuse users about the quantity being displayed.

**Proposed changes:** Use `SLD` consistently as the user-facing color bar label. Review explicit visualize outputs, run quick-look figures, and canonicalize previews for consistency. Preserve exact field names in data artifacts and reports so the actual color source and display transformations remain auditable.

**Acceptance criteria:**

- Relevant SLD color bars use the label `SLD` regardless of the underlying SLD field.
- Color values, normalization, clipping, sampling, and scientific data remain unchanged.
- Data and reports retain the actual color-source field and relevant display policies; this is a label change, not a field rename.
- Update relevant documentation and verify representative figure labels across the affected rendering paths.

**Completion record (2026-09-14):** Centralized SLD color-bar labeling in `cryorole/export/visualization.py`, preserving color-source fields and normalization. Verification: raw/display labels and unchanged color arrays in `tests/export/test_visualization.py`; a generated canonical rotvec projection was visually inspected and showed SLD. See the verification record below.

### OPT-004: Rename the public opacity option to --opacity

**Status:** Completed

**Recorded:** 2026-09-13

**Problem:** The opacity option `--alpha` shares its name with the landscape's Euler alpha coordinate. The name `--opacity` communicates its purpose more clearly.

**Proposed changes:**

- Use `--opacity` in public CLI help, examples, and validation messages.
- Preserve opacity defaults and rendering semantics; document the range from 0 (fully transparent) to 1 (fully opaque).
- Retain `--alpha` as a hidden compatibility alias for existing scripts. Public documentation should introduce only `--opacity`.
- Reject commands that specify both names, including identical values, rather than allowing argument order to determine the result.
- Keep Euler alpha coordinates and their names unchanged. An internal field rename is not required for this public CLI improvement.

**Acceptance criteria:**

- `--opacity` controls the same display behavior as the existing opacity option, with unchanged defaults and valid-value semantics.
- Legacy `--alpha` commands remain functional, while the alias is absent from public help.
- Supplying both option names produces a clear error; invalid opacity values are rejected with guidance using the public name.
- Scientific coordinates, SLD, selections, and export behavior remain unchanged.
- Update CLI documentation and relevant tests for the new name, legacy compatibility, default behavior, invalid values, and conflicting option names.

**Completion record (2026-09-14):** Added public `--opacity` and a hidden, mutually exclusive `--alpha` compatibility alias in `cryorole/cli/parsers.py`; updated validation wording in `cryorole/visualize/service.py`. Verification: aliases, defaults, boundaries, invalid values, conflicting names, help, and end-to-end visualization in `tests/cli/test_optimization_cli.py` and `tests/run/test_optimization_workflow.py`. See the verification record below.

### OPT-005: Organize selection help by mode and improve actionable guidance

**Status:** Completed

**Recorded:** 2026-09-13

**Problem:** Selection options are presented largely as a flat list. Users must infer each mode's purpose, required inputs, applicable controls, units, and conflicts. The brief selection-ID help also does not explain that users should choose a name for each saved subset.

**Proposed changes:**

- Preserve the existing `--mode` command structure. Start help with common inputs, including required `--run-dir` and `--selection-id`, followed by `--space` and `--canonical-id` with their applicability and defaults.
- Provide a concise mode overview, then a dedicated section for each mode with its purpose, required and optional parameters, defaults, units, dependencies, mutually exclusive options, and one complete command example.
- Put output controls in a separate section: `--write-selected-landscape`, `--recompute-sld`, and `--overwrite`. Explain that SLD recomputation requires a selected-derived landscape and preserves parent SLD provenance.
- Use public mode and option names in actionable errors. Missing inputs should identify the mode, required options, and a short example. Audit existing handling of mode-specific options; reject explicitly supplied options that do not apply instead of silently ignoring them. Distinguish explicit arguments from parser defaults during validation.
- Preserve explicit selection naming without an automatic default. Explain that the ID is a user-chosen name, give examples such as `region_01` and `high_sld`, and show the output location. Missing-ID errors should explain how to add the name; existing-ID errors should prioritize choosing another name and explain that `--overwrite` explicitly replaces only that selection. Successful output should show the name, selected count, location, and a subsequent command using that ID. Do not introduce an interactive prompt into ordinary CLI error handling.

**Mode coverage:**

| Mode | Purpose | Required inputs | Additional explanations and controls |
| --- | --- | --- | --- |
| `radius` | Select particles near a center orientation; SO(3) geodesic distance is the default. | `--center` and exactly one of `--radius` or `--radius-rad`. | Explain `--center-representation` and `--metric`; distinguish Euler-center degrees, rotvec-center radians, and independently specified radius units. Identify the selected raw/canonical coordinate space and recorded Euler convention. |
| `threshold` | Select by original SLD bounds. | At least one of `--sld-min` or `--sld-max`. | Explain combined bounds and boundary inclusion. State that selection uses `sld_raw`, independently of display filtering or color scaling. |
| `range` | Select by coordinate-axis bounds. | At least one `--range-bound`. | List supported axis names and units, repeated-bound syntax, open-ended bounds, multi-axis combination rules, representation restrictions, and boundary behavior. Distinguish coordinate ranges from an SO(3) radius. |
| `random` | Select an explicit random fraction from the parent landscape. | `--fraction`. | Explain valid fractions, optional `--seed`, behavior when the seed is omitted, recorded reproducibility information, and the candidate population. Distinguish scientific random selection from visualization sampling. |
| `metadata` | Select source-metadata values or split selections by value. | `--metadata-domain`, `--metadata-column`, and exactly one of `--metadata-value` or `--split-by-value`. | Explain explicit `ref`/`mov` choice, run-time source metadata and recorded row provenance, comma-separated values, and child-selection naming for split mode. |

**Example complete command:**

```bash
cryorole select --run-dir RUN --selection-id high_sld --mode threshold --sld-min 2
```

**Example missing-input guidance:**

```text
Mode 'threshold' requires --sld-min, --sld-max, or both.
Example: --mode threshold --sld-min 2
```

**Acceptance criteria:**

- `cryorole select --help` lets users identify each mode's purpose and all applicable public controls without inferring relationships from an ungrouped list. Include a complete example for every mode.
- Check documented units, defaults, range rules, seed behavior, and split naming against the implementation. Do not invent scientific semantics to simplify help.
- Required-input, conflict, and inapplicable-option errors use public CLI terminology and suggest a corrective action. Help remains available without supplying run inputs or a selection ID and has no artifact side effects.
- Selection-ID guidance explains naming and preserves explicit overwrite protection. Scientific evaluation, full-parent selection, standard Selection artifacts, and export provenance remain unchanged.
- Update relevant CLI/workflow documentation and tests for grouped help, valid examples, missing inputs, option conflicts/applicability, and selection-ID guidance. Keep validation in the appropriate shared service boundary where applicable.

**Deferred extension:** Consider focused help such as `cryorole select --mode radius --help` only if full help becomes unwieldy. A new subcommand hierarchy is not part of this item.

**Completion record (2026-09-14):** Grouped mode help and complete examples in `cryorole/cli/parsers.py`; added shared request validation and typed selection counts in `cryorole/select/service.py`; added actionable success/export guidance in `cryorole/cli/commands/downstream.py`. Removed the obsolete automatic metadata-prefix fallback. Selection IDs must be names rather than paths, including when overwrite is requested. Verification: all six help examples, missing/conflicting/inapplicable options, explicit defaults, name validation, overwrite protection, shell quoting, standard artifact counts, and source-row export in `tests/cli/test_optimization_cli.py`, `tests/cli/test_cli_workflow.py`, and `tests/run/test_optimization_workflow.py`. See the verification record below.

### OPT-006: Support CryoSPARC source metadata selection

**Status:** Completed

**Recorded:** 2026-09-14

**Scope:** Small Optimization, bounded format-support extension. Implementation authorized on 2026-09-14.

**Problem and evidence:** CryoSPARC `.cs` is supported by run and export, but `_load_run_source_metadata_for_selection` in `cryorole/select/service.py` rejects non-STAR sources. The existing CS reader can read structured fields, but currently materializes the full table, including object-backed vector columns. During README example validation, selection by `alignments3D/class` failed with the STAR-only error. README now states this limitation; CLI help and the CLI reference do not make it explicit.

**Recommended direction:** Extend the existing metadata mode to CS without adding a separate command or converting CS to STAR. Preserve STAR behavior and share selection evaluation and standard artifact writing. Treat this as format parity for categorical source fields, not a general metadata query language.

**Proposed changes:**

- Keep `--metadata-domain ref|mov`, `--metadata-column`, `--metadata-value`, and `--split-by-value`. Resolve the original source from the run and use recorded ref/mov source-row IDs; do not rematch by UID or infer the domain.
- Add a small typed source-column adapter at the shared service/I/O boundary. For CS, validate a structured array and read the requested field with compact arrays (memory mapping where appropriate), avoiding conversion of every metadata field to a DataFrame or per-row objects.
- Initially support scalar categorical fields: integer, boolean, and text/byte-string values. Preserve integer precision, including uint64 UIDs above 2^53; use explicit, lossless text/byte handling. Reject unsupported vector, nested, object, and floating-point fields with actionable errors until their comparison policies are defined. Do not flatten pose/shift vectors or silently convert values through float.
- Use actual source values without renumbering classes. A first example is `alignments3D/class` if present in the run-time CS file; do not assume every CS file contains it.
- Support multiple values and per-value splitting over the full parent landscape's recorded source rows. Reuse standard child-selection naming and overwrite protection; ensure distinct typed values cannot silently collide after name sanitization.
- Validate source identity through the existing shared identity policy before consuming metadata, and check field availability, supported dtype, and source-row bounds before writing selections. Define missing/invalid-value handling explicitly and report the requested/resolved values and counts.
- Preserve source dtype and vector fields during subsequent export. Record domain, source identity/format, field, value interpretation, source-row policy, candidate count, and selected count in the existing policy/report structure.
- Update README, CLI help/reference, the CryoSPARC workflow, and the relevant architecture contract together with implementation. Until support ships, retain the README limitation.

**CS command example:**

```bash
cryorole select --run-dir RUN --selection-id ref_classes_01 --mode metadata --metadata-domain ref --metadata-column alignments3D/class --metadata-value 0,1
```

Requires a run-time reference CS source containing the named scalar field and values.

**Acceptance criteria:**

- CS value selection, multi-value union, and per-value splitting create standard Selection artifacts that visualize and export successfully for either explicit domain.
- Tests cover reordered and partially matched inputs, source-row backtracking, missing columns, invalid row IDs, unsupported vector fields, byte/text values, boolean values, and large uint64 values without precision loss.
- Split tests cover deterministic child names, sanitized-name collisions, existing outputs, and pre-write validation. Invalid source content or requests must not leave misleading completed selections.
- End-to-end CS exports contain exactly the selected source rows with original dtype and vector-valued fields intact; source files, parent RO/SLD, and raw/canonical artifacts remain unchanged.
- Existing STAR metadata selection and export regression tests continue to pass. Include a representative wide CS-table memory check to confirm that unrelated vector fields are not materialized as per-row Python objects.

**Implementation decisions (2026-09-14):** CS supports scalar integer, boolean, and UTF-8 text. Empty strings are missing and excluded with counts; invalid UTF-8, unsupported dtypes, and invalid source rows fail. CS split is capped at 100 non-missing groups, with explicit value selection suggested above that limit; no new override option is added. STAR normalization and split cardinality remain unchanged. Float tolerance/ranges, vector-component expressions, escaped comma-containing values, external annotation joins, and new field-discovery CLI options remain deferred.

**Implementation record (2026-09-14):** Implemented mapped CS scalar-column loading, typed matching, source SHA-256 validation, and standard selection/export reuse. Removed three duplicated STAR normalization/column helpers from the service. Metadata evaluation, derived-landscape count requirements, split names, and existing targets are checked before replacing artifacts. New CS requests require verified identities; legacy STAR compatibility is retained and recorded explicitly.

**Changed files:** `cryorole/io/readers/cs_reader.py`, `cryorole/select/metadata.py` (new), `cryorole/select/service.py`, `cryorole/select/selectors.py`, `cryorole/models/policies.py`, `cryorole/provenance/source_identity.py`, `cryorole/cli/parsers.py`, and `tests/select/test_cs_metadata.py` (new). Updated `README.md`, `docs/architecture.md`, `docs/production_run_plan.md`, `docs/roadmap.md`, `docs/cli_reference.md`, `docs/cryosparc_workflow.md`, and this backlog.

**Verification:** New CS and documentation tests: 35 passed. Ruff passed; links in six documents and 32 complete command examples checked. A 100,000-row wide-CS column-preparation smoke comparison produced identical selected counts: full-table loading 0.0239 s / 139,190,272 bytes peak RSS; mapped scalar preparation 0.0129 s / 124,878,848 bytes peak RSS. This is one local column-preparation measurement, not an end-to-end benchmark; details are in `.codex_tmp/opt006_memory_check.json`. Full regression: `python -u -m pytest -q -ra --tb=short` — 718 passed, 2 skipped in 324.66 s. The skips require external ChimeraX and FFmpeg/FFprobe installations; log: `.codex_tmp/opt006_full_tests.log`.

**Completion record (2026-09-14):** Implemented and verified. Native CS value selection and splitting now support the agreed initial scalar types through the existing CLI. Deferred query extensions remain outside this completed scope.

## Verification Record: 2026-09-14

- Full regression: `python -u -m pytest -q --disable-warnings --tb=short` — 680 passed, 2 skipped. After the final small validation/quoting adjustments, the focused CLI, diagnostic, workflow, and generated-viewer checks passed: 58 passed, 159 deselected. Final persistence and documentation checks: 8 passed. Logs are under `.codex_tmp/optimization_20260913/` (`full_tests.log`, `final_targeted_tests.log`, and `persistence_docs_tests.log`).
- Static analysis: `python -m ruff check cryorole tests` passed.
- Additional test file: `tests/visualize/test_offline_viewer.py` exercises the actual generated JavaScript with Node and a canvas stub; this is separate from browser acceptance.
- Updated authoritative/user documents: `docs/architecture.md`, `docs/production_run_plan.md`, `docs/roadmap.md`, `docs/cli_reference.md`, `docs/faq.md`, and `docs/workflow_ux.md`, plus this backlog.
- A local 100,000-row diagnostic comparison measured approximately 0.063 s / 114.4 MiB peak RSS before and 0.070 s / 112.1 MiB after, with identical coincidence counts. This single measurement is a smoke check, not a performance guarantee; details are in `.codex_tmp/optimization_20260913/diagnostic_benchmark.json`.
- Remaining acceptance work: OPT-002 requires actual offline Chrome rendering and interaction checks. Computer Use stopped because it could not confidently determine the local browser URL for policy enforcement; no browser acceptance claim is made.
