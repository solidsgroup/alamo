# Task 005: Phase 5.2 — performance-regression tracking

## Goal
A tool that records the roadmap's standing perf metrics per commit on the canonical
case, so regressions are visible over time: launches/step, avg kernel duration,
`cudaStreamSynchronize` fraction, and wall/step (GPU vs CPU node). Append one row
per run keyed by git SHA + date to a CSV; document usage.

## Context
- nsys tooling exists: `source benchmark/local_cuda_env.sh`; nsys at
  `.local/nsight/opt/nvidia/nsight-systems/2026.1.3/target-linux-x64/nsys`;
  `nsys stats --report cuda_api_sum --report cuda_gpu_kern_sum`. The handoff
  (`~/Desktop/SESSION_HANDOFF_2026-06-20.md`, section 5) has a full nsys capture
  command for the parity config.
- Metrics definitions: launches/step = cudaLaunchKernel count / steps;
  sync fraction = cudaStreamSynchronize time / total CUDA-API time; kernel avg from
  cuda_gpu_kern_sum.
- This box can run nsys but a CUDA runner in CI may not exist; make the script
  standalone (run by hand or wired to CI later) and document wiring.

## Files to read first
- `~/Desktop/SESSION_HANDOFF_2026-06-20.md` (section 5: nsys command + metric values)
- `benchmark/local_cuda_env.sh`, any existing `benchmark/*tinyprofiler*`/nsys parser
- `benchmark/baseline_suite.py` (for how the canonical case is launched)

## Files allowed to modify
- `benchmark/perf_regression_track.py` (new)
- `benchmark/PERF_TRACKING.md` (new)

## Files NOT allowed to modify
- `.github/workflows/**` (task 004 owns CI); any source.

## Implementation steps
1. Write `perf_regression_track.py`: run (or accept an existing nsys report for) the
   canonical case, parse cuda_api_sum + cuda_gpu_kern_sum, compute the 4 metrics,
   and append a row `{date, git_sha, config, launches_per_step, kernel_avg_us,
   sync_frac, wall_per_step_gpu, wall_per_step_cpu}` to
   `benchmark/perf_regression.csv` (create header if absent). Include a `--compare`
   mode that flags >X% regression vs the last recorded row for the same config.
   Resolve nsys via the same precedence as the analysis suite ($NSYS → .local → PATH).
2. Write `PERF_TRACKING.md`: what it measures, how to run, how to read the CSV, how
   to wire it into the correctness workflow later (a gated GPU job).
3. Make it not hard-fail when nsys/GPU is absent — print a clear "skipped: no GPU"
   and exit 0 in that mode, so it is CI-safe on CPU-only runners.

## Invariants
Standalone + CI-safe (graceful skip without GPU). Does not touch CI files.

## Build and test commands
`python3 benchmark/perf_regression_track.py --help` should work; a `--dry-run` or
no-GPU path should exit cleanly.

## Expected result
A working perf-tracking script + doc, CSV schema defined, CI-safe.

## Non-goals
Editing CI workflows; correctness gating (task 004).

## Stop conditions
If GPU is unavailable to validate a real capture, validate the parser against a
saved nsys stats CSV or the metric values quoted in the handoff, and note it.

## Final report: write results/005-RESULT.md
