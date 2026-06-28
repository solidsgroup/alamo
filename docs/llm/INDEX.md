# docs/llm — chamber-gpu map (read this first)

Absolute-truth root for LLM planning on the `chamber-gpu` branch. Read this
file plus `CONVENTIONS.md` at session start; everything else is opt-in via
the links below. Nothing here duplicates content that lives elsewhere — it
points at it.

## Start here

- **`CURRENT.md`** — what's in flight right now. Read this second, always.
- **`CONVENTIONS.md`** — how to read/work/test/write against this map. Read once per session.

## Strategy / planning

- **`ROADMAP.md`** — the master plan: phases P0-P5, decision trees D1-D4, gates G0-G5, guiding principles. Moved from `~/Desktop/GPU-OPT-ROADMAP.txt` 2026-06-22 so it's versioned with the branch.
- `benchmark/GPU_BRANCH_GUIDE.md` — top-level build/run/branch-policy map (left in place; pre-existing).
- `docs/agent_plans/` — per-phase working folders (`PLAN.md` + `tasks/*.md` + `results/*.md`), one subfolder per dated effort. Left in place; this is where implementation detail and strategy-from-prior-work actually live:
  - `20260620-gpu-phase1-and-lever25/` — Phase 1 elastic disposition + lever 2.5 field packing
  - `20260620-gpu-phase3-regime-scaling/` — Phase 3 regime scaling / crossover hunt (currently active, see `CURRENT.md`); includes `TEST_CASES.md`, the what/why for the `input_3d_flame_{128,256,512}` sweep sizes
  - `20260621-gpu-phase4-5/` — Phase 4 dispatch decision + Phase 5 hardening
  - `20260622-analysis-suite-backfill/` — backfilled record for the `analysis/` six-phase suite rewrite and its five untracked results bundles (one of which, `phase3_nova_bundle_20260620_212331/`, is a **failed** NOVA attempt, not a success — see `tasks/002-results-bundles.md`)
- `docs/gpu_elastic_device_port_plan.md`, `docs/gpu_device_capture_conventions.md`, `docs/gpu_safe_ic_bc_matrix.md` — standing technical references (left in place).

## Past progress / completed phases (changelog)

- `changelog/2026-06-20-session-handoff.md` — Phase 2.2 commit + Phase 1/2.5 handoff. Moved from Desktop 2026-06-22.
- `changelog/2026-06-26-3d-elastic-gpu-fix.md` — **3D elastic GPU crash root-caused + fixed**: `F.inverse().transpose()` (chained Eigen expr) faults on device; 2 sites in `NeoHookean.H`. Same bug class as the elixir UAF. First time the 3D combined flame+elastic device path was exercised.
- `changelog/2026-06-26-gpu-test-suite-fixes.md` — **GPU test suite 5P/4F → 9P/0F**: 4 failures = 4 different bugs. 2 real source defects fixed in `Integrator.cpp` (`Restart()` node-fab OOB segfault; headerless restart `thermo.dat`); 3 stale-input decks rewritten; C3 MLMG stiffness-contrast divergence. Full writeup `benchmark/GPU_TEST_SUITE_FIXES.md`.
- `benchmark/PHASE1_ELASTIC_DISPOSITION.md` — D1 verdict (left in place; ~7 historical docs cite this path).
- `benchmark/PHASE3_3D_READINESS.md`, `benchmark/PHASE3_NOVA_TESTING_PROGRESS.md`, `benchmark/PHASE3_R3_crossover.md` — Phase 3 reports (left in place).
- `benchmark/PHASE4_R4_dispatch.md` — D4 dispatch decision (left in place).
- `benchmark/PHASE5_BRANCH_DONE.md` — branch definition-of-done checklist (left in place).
- `benchmark/G0_BASELINE_OF_RECORD.md` — Phase 0 baseline-of-record report (left in place).

## Performance / CPU vs GPU data

- **`perf/2026-06-20-gpu-port-report.md`** — full engineering + perf report: de-virtualization strategy, 3-campaign CPU-vs-GPU characterization, Nsight Systems device analysis. Moved from Desktop 2026-06-22.
- `benchmark/PERF_TRACKING.md` + `benchmark/perf_regression.csv` — per-commit perf-regression tracking tool (left in place; append-only CSV).
- `benchmark/GPU_TEST_PERF_TRACKING.md` — per-test wall-time/throughput baseline for the `tests/GPU/` suite, one iteration block per fix-set (started 2026-06-26 post-fix baseline). Use to catch test-suite perf regressions across fixes.
- `analysis/README.md` + `analysis/results*/REPORT.md` — the wallclock/IO/perf-stat/flamegraph/nsys analysis suite and its run bundles (left in place).

## Versioning

- **`VERSIONS.md`** — semver ledger for `chamber-gpu` eras (`gpu-vX.Y.Z` -> commit -> perf report -> changelog entry). Standing D1-D4 decisions also indexed here.

## Out of scope

`~/Desktop/{alamo_sim_status,lessons_learned,campaign_results_tracker,multifin_campaign,parameter_exploration_campaign,AUTONOMOUS_CAMPAIGN_GUIDE,stress_campaign_report_*}.md` and the `LAUNCH_PHASES_2-5.txt` / `AUTONOMOUS_CAMPAIGN_READY.txt` / `CAMPAIGN_QUICK_START.txt` files are a **separate, unrelated** propellant parameter-sweep campaign (sims 030-085) — not `chamber-gpu`/GPU-port work. Deliberately left untouched.
