# Result: Task 002 — results bundles backfill

## Summary
Identified and recorded the source/status of all five untracked `analysis/results_*`
bundles. Four are valid artifacts backing existing reports; one
(`phase3_nova_bundle_20260620_212331/`) is a **failed** early NOVA attempt and was
at risk of being confused with the later, successful 2026-06-21 crossover run.

## Tests run and results
- `du -sh` on all five bundles -> sizes recorded in task 002 (largest:
  `results_wide_grid/` at 5.6G, carrying full nsys traces for two configs).
- Read `phase3_nova_bundle_20260620_212331/RESULTS_ANALYSIS.md` in full -> confirms,
  in its own words, "No usable GPU/CPU crossover timing can be computed from this
  bundle" — all GPU jobs hit CUDA error 700 before step 1.
- `git status --short analysis/results/` -> empty; confirmed via `.gitignore:5`
  (`results`) that the coarse-campaign bundle `analysis/results/` is intentionally
  gitignored, which is why it doesn't appear alongside the other four as untracked.

## Issues found
- `analysis/results_wide_grid/` is 5.6G of nsys traces sitting untracked in the repo
  working tree. Not a correctness issue, but worth a deliberate decision (gitignore
  it like `analysis/results/`, or move it off-repo) rather than leaving it as
  accidental untracked state. Not actioned here — flagging only.
- No local bundle exists for the successful 2026-06-21 NOVA crossover run (jobs
  11160767-11160774) that `benchmark/archive/PHASE3_R3_crossover.md` cites — those numbers
  trace back to NOVA job logs, not a locally-extracted `analysis/` bundle. If those
  logs are still on NOVA, pulling them into a `phase3_nova_bundle_20260621_*/` would
  close that gap; not done here (out of scope, no NOVA access from this task).

## Deviations from task
None.

## Follow-up needed
- Decide on `results_wide_grid/`'s size (gitignore vs. relocate vs. prune the
  `.nsys-rep` binaries once their CSV summaries are confirmed sufficient).
- Optionally retrieve/bundle the 2026-06-21 successful NOVA crossover logs for the
  same archival treatment as the failed 2026-06-20 attempt.
