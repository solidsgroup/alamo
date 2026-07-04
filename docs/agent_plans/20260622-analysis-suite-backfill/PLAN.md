# Plan: Backfill task records for the analysis suite rewrite + Phase-3 results bundles

## User goal
The `analysis/` directory was rewritten (six-phase profiling suite, `01_wallclock.py`
... `06_report.py`, `lib/`, `run_all.sh`) and several result bundles were generated
(`results_fine_grid/`, `results_wide_grid/`, `results_phase22/`,
`results_phase22_nsub1/`, `phase3_nova_bundle_20260620_212331/`) across 2026-06-19
through 2026-06-22, all of it still untracked in git and with no task record under
`docs/agent_plans/`. This plan backfills that record per `docs/llm/CONVENTIONS.md`
("if it isn't written to one of these files, it didn't happen") and writes the
standalone what/why doc for the `input_3d_flame_{128,256,512}` test cases that
`docs/llm/CURRENT.md` flagged as missing.

## Note on backfilling
This is written *after* the work, not before — task files below describe what was
actually built (verified by reading the code/output), not a prospective plan. Test
commands are real and were re-run during this backfill to confirm they still pass.

## Scope
- `tasks/001-analysis-suite-rewrite.md` — the six-phase `analysis/` suite itself.
- `tasks/002-results-bundles.md` — the five untracked results bundles, what each one
  is and what it proved (or didn't).
- `docs/agent_plans/20260620-gpu-phase3-regime-scaling/TEST_CASES.md` — what/why for
  the `input_3d_flame_{128,256,512}` inputs (owned by the existing Phase-3 folder,
  not duplicated here).

## Out of scope
No code changes. No re-running of NOVA jobs. No re-running the local analysis suite
end-to-end (it needs production CPU/GPU logs that aren't present on this box right now);
verification is limited to syntax/structure checks that don't require those logs.
