# Task 004: Phase 5.1 — correctness CI gate for chamber-gpu

## Goal
Add a CI workflow that runs on pushes to `chamber-gpu` and re-arms the correctness
tripwires: (a) a golden bit/tolerance compare (no-fast-math) and (b) a NaN-flag
assertion build. The branch is never merged, so this guards the branch itself.

## Context
- Existing CI patterns: `.github/workflows/linux.yml`, `performance.yml` (self-
  hosted `scooter` runner; `./configure --dim --comp=g++ && make -j8`;
  `scripts/runtests.py`).
- Golden harness: `benchmark/baseline_suite.py` (subcommand `check`) compares
  `thermo.dat` against `benchmark/baseline_references/` (canonical_step1/2,
  eta_expression_step1).
- The full golden compare is GPU-vs-CPU; a GPU CI runner may not exist. Ship the
  CPU-side correctness that runs on a standard runner now, and clearly document the
  GPU leg as runner-gated (guarded behind a job that only runs if a CUDA runner
  label is available).

## Files to read first
- `.github/workflows/linux.yml`, `.github/workflows/performance.yml`
- `benchmark/baseline_suite.py` (its CLI + what `check` needs)
- `benchmark/baseline_references/` (layout)

## Files allowed to modify
- `.github/workflows/chamber-gpu-correctness.yml` (new)
- `benchmark/ci_golden_compare.sh` (new)

## Files NOT allowed to modify
- Any existing workflow; any source; `Makefile`/`configure`.

## Implementation steps
1. Write `benchmark/ci_golden_compare.sh`: build the no-fast-math correctness build
   (CPU path is fine for the bit-compare leg), run `baseline_suite.py check`, exit
   nonzero on any mismatch; also do a NaN-flag assertion smoke (a short run that
   would abort if the device error flag / NaN tripwire fires). Make it work from a
   clean checkout; parameterize the GPU vs CPU leg with an env flag.
2. Write `.github/workflows/chamber-gpu-correctness.yml`: trigger `on: push:
   branches: [chamber-gpu]` (and workflow_dispatch). Jobs:
   - `golden-cpu` (standard runner): configure+build CPU, run ci_golden_compare.sh
     in CPU mode → fails the build on mismatch.
   - `golden-gpu` (gated): `if:` a CUDA runner is available (e.g. runs-on a
     `cuda`/`scooter` label); document that it is skipped where no GPU runner exists.
3. Keep YAML consistent with the repo's existing style (continue-on-error usage,
   step naming).

## Invariants
Must not weaken existing CI; new workflow only. Correctness gate must FAIL on a
real golden mismatch (no silent continue-on-error on the compare step).

## Build and test commands
Validate locally only by `bash -n benchmark/ci_golden_compare.sh` and a YAML
lint mentally; do not run a full build here.

## Expected result
A working correctness workflow + helper script that re-arm the golden + NaN gates
for chamber-gpu, with the GPU leg cleanly runner-gated.

## Non-goals
Setting up runners; perf tracking (task 005 owns that).

## Stop conditions
If `baseline_suite.py check` CLI differs from assumptions, adapt the script to the
real CLI (read the file) and note it in RESULT.

## Final report: write results/004-RESULT.md
