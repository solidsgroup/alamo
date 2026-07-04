# Task 002 (backfilled): the five untracked Phase-2/3 results bundles

## Task
Document what each of the five untracked `analysis/results_*` / bundle directories
is, what campaign produced it, and whether it succeeded — so `docs/llm/CURRENT.md`
no longer lists them as unexplained loose state.

## Why
These are the raw artifacts backing claims already written up in
`docs/llm/perf/2026-06-20-gpu-port-report.md` and `benchmark/archive/PHASE3_R3_crossover.md`.
Without this record, a future agent can't tell which bundle is "the" source for a
cited number, or that one of them is a *failed* run, not a missing-but-fine one.

## What each bundle is

- **`analysis/results_fine_grid/`** (88K, 2026-06-20) — Campaign 2 ("deep fine-grid",
  `max_level=5`, `n_cell=64`) from the coarse/deep/wide 2D characterization. Backs
  `perf/2026-06-20-gpu-port-report.md` §6.2. Contains `REPORT.md`, `wallclock.{md,csv,json}`,
  `chart_wallclock.svg`. Headline: GPU 2.87x slower than CPU; 0 elastic solves reached.

- **`analysis/results_wide_grid/`** (5.6G, 2026-06-20) — Campaign 3 ("wide-shallow",
  `max_level=2`, `n_cell=512`, same 2048^2 finest cell as the deep run) plus the nsys
  device-level capture for *both* deep and wide (`nsys_{wide,deep}.nsys-rep/.sqlite`,
  `stats_{wide,deep}_cuda_{api,gpu_kern}_sum.csv`). Backs report §6.3 and §7 (the
  6.7x launch-count / 9.8x wall-per-step finding). Large because it carries full nsys
  traces for two 30-step captures, not just summaries — candidate for gitignore if kept
  around (see Issues below).

- **`analysis/results_phase22/`** (145M, 2026-06-20) — nsys "parity config" capture
  (30 steps, phase-field only, `wide_512_bf32_mgs128` overrides) backing Phase 2.2's
  D2 verdict in `changelog/2026-06-20-session-handoff.md` (11,592 launches/step, 4,170
  syncs/step, 18.5% sync fraction). Contains the plotfile (`plt/`), `run.log`,
  `nsys_parity.{nsys-rep,sqlite}`, `stats_cuda_{api,gpu_kern}_sum.csv`.

- **`analysis/results_phase22_nsub1/`** (62M, 2026-06-20) — same parity config with
  `amr.nsubsteps=1` instead of 2, isolating subcycling's effect on launch count
  (11,592 -> 2,026/step, 5.72x). Backs the same D2 verdict's "subcycle-multiplied, not
  ghost-exchange-bound" conclusion.

- **`analysis/phase3_nova_bundle_20260620_212331/`** (4.3M, captured 2026-06-20 21:23,
  repo commit `23c924721`) — **a failed early NOVA 3D attempt**, not a success. All
  captured 1-GPU and 2-GPU A100 jobs aborted with `CUDA error 700` (illegal memory
  access) before completing a single timestep, for all three sizes (128/256/512); the
  512 CPU-node baseline was OOM-killed. See its own `RESULTS_ANALYSIS.md` for the
  per-job table. **This is not the source of the "39x/70x crossover" numbers** in
  `benchmark/archive/PHASE3_R3_crossover.md` — those came from later NOVA jobs
  (11160767-11160774, 2026-06-21), which ran after whatever fixed this CUDA-700
  failure, and no local bundle for that later, successful run exists under
  `analysis/` — the crossover report's numbers are sourced from the NOVA job logs
  directly, not from a locally-extracted bundle.

## Test command
```bash
du -sh analysis/results_fine_grid analysis/results_wide_grid analysis/results_phase22 \
       analysis/results_phase22_nsub1 analysis/phase3_nova_bundle_20260620_212331
git status --short analysis/results/   # confirm analysis/results/ (the coarse-campaign
                                        # bundle) is gitignored via .gitignore:5 "results",
                                        # which is why it doesn't show as untracked above
```

## Acceptance criteria
Each bundle's owning campaign/report is identified; the failed NOVA bundle is
explicitly distinguished from the successful 2026-06-21 crossover run it is often
adjacent to in conversation. Re-verified during backfill: `RESULTS_ANALYSIS.md`'s own
text confirms "No usable GPU/CPU crossover timing can be computed from this bundle."

## Write results to: results/002-RESULT.md
