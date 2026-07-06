# AGENTS.md Compatibility Guardrails

The authoritative contributor guardrails live at the repository root in
`AGENTS.md`. This docs-side copy exists so documentation contract tests and
older contributor prompts can continue to read `docs/AGENTS.md`.

Current production-run persistence guardrail:

- `raw_landscape.npz` is the production machine-readable raw landscape.
- `raw_landscape.csv` is the user-facing raw table.
- Full object-record landscape JSON is debug-only and opt-in.
- Do not restore `landscape.json` as the default production persistence contract.

For complete scientific and engineering rules, read root `AGENTS.md` first,
then `docs/architecture.md`, then `docs/production_run_plan.md`.
