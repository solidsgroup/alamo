# Result: Task 005 — Phase 5.2 performance-regression tracking

## Status: DONE

## Files created
- `benchmark/perf_regression_track.py` (new, executable via `python3`)
- `benchmark/PERF_TRACKING.md` (new)

No other files modified. `.github/workflows/**`, source, `Makefile`, `configure`
untouched (those pre-existing uncommitted diffs in the working tree were not
created or touched by this task).

## What was built

`benchmark/perf_regression_track.py`:
- Resolves `nsys` via the same precedence as the rest of the analysis tooling:
  `$NSYS` -> `.local/nsight/opt/nvidia/nsight-systems/*/target-linux-x64/nsys`
  -> `PATH`.
- `--capture`: runs the canonical parity config (64^3 base grid, max_level=3,
  blocking_factor=32, max_grid_size=128, grid_eff=0.9, nsubsteps=2,
  `elastic.tstart=1e9`) for `--steps` (default 30) under
  `nsys profile --trace=cuda,nvtx --sample=none --cpuctxsw=none`, exactly
  matching the command in `~/Desktop/SESSION_HANDOFF_2026-06-20.md` section 5,
  then `nsys stats --report cuda_api_sum --report cuda_gpu_kern_sum --format csv`.
- `--api-csv`/`--kern-csv`: parses an already-captured stats CSV pair instead
  (no GPU/nsys needed) — used for backfilling history or validating the parser.
- Computes 4 metrics: `launches_per_step` (cudaLaunchKernel Num Calls / steps),
  `kernel_avg_us` (sum of kernel Total-Time-ns / sum of Instances, from
  `cuda_gpu_kern_sum`), `sync_frac` (cudaStreamSynchronize Total-Time-ns /
  total cuda_api_sum Total-Time-ns), and accepts optional caller-supplied
  `--wall-gpu-ms`/`--wall-cpu-ms`.
- Appends one row keyed by `date, git_sha, config, steps, launches_per_step,
  kernel_avg_us, sync_frac, wall_per_step_gpu_ms, wall_per_step_cpu_ms, notes`
  to `benchmark/perf_regression.csv` (header written on first use; file is not
  pre-created or committed by this task).
- `--compare`: looks up the most recent prior row for the same `config`
  (excluding the current SHA) and flags any of the 5 numeric metrics that
  regressed (got worse) by more than `--threshold-pct` (default 10%); exits 1
  if a regression is flagged, 0 otherwise. `--no-record` allows a dry-run.
- **CI-safety**: if no usable GPU (`nvidia-smi` absent/fails) or no resolvable
  `nsys`, prints `skipped: no GPU` and exits 0 without touching the CSV. If a
  real capture attempt fails for another reason, prints
  `skipped: GPU capture unavailable/failed (...)`, cleans up the partial
  capture directory, and still exits 0 — it never hard-fails a CPU-only run.

`benchmark/PERF_TRACKING.md`: documents the 4 metrics with their exact
definitions, how to run fresh captures vs. parse existing CSVs, the
`--compare` regression workflow, the CSV schema, the CI-safety contract, and a
documented (not wired) example of how to add this as a future gated GPU job in
`.github/workflows` without touching that directory now.

## Validation performed

1. `python3 benchmark/perf_regression_track.py --help` — works, full usage text
   renders.
2. **Parser correctness against real nsys data**: ran the tool with
   `--api-csv analysis/results_phase22/stats_cuda_api_sum.csv --kern-csv
   analysis/results_phase22/stats_cuda_gpu_kern_sum.csv --steps 30` (a real
   prior capture of the exact parity config, referenced in the handoff). Result:
   `launches_per_step=11591.83`, `sync_frac=18.49%` — matching the handoff's
   quoted **11,592 cudaLaunchKernel/step** and **18.5% sync fraction** to
   rounding. This is strong evidence the CSV parsing and metric formulas are
   correct.
3. **`--compare` regression flagging**: seeded a synthetic "previous" row with
   better (lower) metric values than the real phase22 capture, confirmed
   `--compare` correctly computed percentage deltas, flagged `launches_per_step`
   and `sync_frac` as `REGRESSION` (>10% threshold) while leaving
   `kernel_avg_us` (9.2% < 10%) as `ok`, and returned exit code 1.
4. **Live end-to-end `--capture` run**: this box has a real GPU + `.local/nsight`
   nsys + a built `bin/alamo_gpu-2d-cuda86-g++` binary. `--capture --steps 5`
   actually ran nsys against the live GPU, captured, parsed, and recorded a row
   — full pipeline exercised for real, not just the parser. (This artifact was
   cleaned up afterward, see below — `--no-record` not used in this run so a
   real `perf_regression.csv`/`perf_regression_runs/` were produced and then
   deleted as test cleanup, since the task says not to leave stray generated
   artifacts.)
5. **No-GPU skip path**: ran the script with `nvidia-smi` excluded from `PATH`
   (`env -i PATH=<minimal>`) — printed exactly `skipped: no GPU`, exit code 0,
   no CSV or capture directory written.
6. **Capture-failure (not absence) path**: with `nvidia-smi` present but `nsys`
   forced to a real-but-broken invocation (stripped MPI/SSH environment caused
   the GPU binary's `mpiexec` launch to fail inside the nsys-wrapped run),
   confirmed the tool catches this distinctly, prints
   `skipped: GPU capture unavailable/failed (...)`, removes the partial capture
   directory, and exits 0 — does not leak artifacts or hard-fail.

All test artifacts (`/tmp/*.csv`, `/tmp/fakebin*`, the one real
`benchmark/perf_regression.csv` + `benchmark/perf_regression_runs/` produced
during live validation) were deleted afterward. Final `git status` for
`benchmark/` shows only `PERF_TRACKING.md` and `perf_regression_track.py` as
new files attributable to this task (plus pre-existing untracked artifacts
from other parallel tasks/sessions, unrelated to this one).

## Notes / deviations from the task spec
- The CSV's `notes` field doubles as the place to record "is the working tree
  dirty" context since `git_sha` alone doesn't capture uncommitted diffs —
  documented in PERF_TRACKING.md rather than adding a `dirty` boolean column,
  to keep the schema exactly as specified in the task (`{date, git_sha, config,
  launches_per_step, kernel_avg_us, sync_frac, wall_per_step_gpu, wall_per_step_cpu}`
  plus `steps` and `notes` added for normalization/annotation — both additions
  are backward-compatible, additive columns).
- `wall_per_step_gpu_ms`/`wall_per_step_cpu_ms` are optional and caller-supplied
  rather than self-measured, because this tool's job is the nsys-derived
  metrics; baseline_suite.py/benchmark_gpu_cpu.sh already measure wall-clock
  and can feed it in via `--wall-gpu-ms`/`--wall-cpu-ms`. Documented in
  PERF_TRACKING.md.
- GPU was available on this box, so validation went beyond "parser against a
  saved CSV" (the stop-condition fallback) to a full live capture — both paths
  are validated.

## Suggested follow-ups (not in scope here)
- Task 004 (or a future task) should decide whether/when to wire
  `perf_regression_track.py --capture --compare` into a gated GPU-runner CI
  job, per the documented (non-binding) example in PERF_TRACKING.md.
- Consider committing an initial `benchmark/perf_regression.csv` with a few
  backfilled historical rows (e.g. from `analysis/results_phase22/`,
  `analysis/results_phase22_nsub1/`) if the user wants tracked history from
  day one — not done here since this task does not commit anything.
