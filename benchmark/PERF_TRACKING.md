# Performance-regression tracking (Phase 5.2)

`benchmark/perf_regression_track.py` records the GPU-OPT-ROADMAP's standing
perf metrics for the Flame canonical case, per git commit, so regressions are
visible over time. It appends one row per run to `benchmark/perf_regression.csv`
(created on first use; not checked in by this tool — see "Where rows live"
below).

This tool is **standalone**: it is not wired into CI by this task. Task 004
(`benchmark/ci_golden_compare.sh` + `.github/workflows/chamber-gpu-correctness.yml`)
owns correctness gating; wiring this tool into a future GPU-runner perf job is
a follow-up (see "Future CI wiring").

## What it measures

Four metrics, computed from an `nsys` capture of the canonical parity config
(64^3 base grid, max_level=3, blocking_factor=32, max_grid_size=128, grid_eff=0.9,
nsubsteps=2 — see `SESSION_HANDOFF_2026-06-20.md` section 5 and
`analysis/results_phase22/`):

| Metric                  | Definition                                                                                       | Source report      |
|--------------------------|---------------------------------------------------------------------------------------------------|---------------------|
| `launches_per_step`      | `cudaLaunchKernel` "Num Calls" / steps                                                            | `cuda_api_sum`      |
| `kernel_avg_us`          | sum(kernel "Total Time (ns)") / sum(kernel "Instances"), converted to microseconds                | `cuda_gpu_kern_sum` |
| `sync_frac`              | `cudaStreamSynchronize` "Total Time (ns)" / sum of all `cuda_api_sum` row "Total Time (ns)"       | `cuda_api_sum`      |
| `wall_per_step_gpu_ms`   | optional — supplied by the caller via `--wall-gpu-ms` (this tool does not itself time a GPU run)  | n/a                 |
| `wall_per_step_cpu_ms`   | optional — supplied by the caller via `--wall-cpu-ms`, for the matching CPU run                   | n/a                 |

`launches_per_step` and `sync_frac` are the two metrics the roadmap's Phase 2
analysis (D2) tracked explicitly: parity-config capture measured **11,592
cudaLaunchKernel/step** and **18.5% cudaStreamSynchronize fraction of CUDA-API
time**. Running this tool against that exact capture
(`analysis/results_phase22/stats_cuda_{api,gpu_kern}_sum.csv`, 30 steps)
reproduces `launches_per_step=11591.83` and `sync_frac=18.49%` — matching the
quoted values up to rounding.

## How to run

### Fresh capture (needs a GPU + nsys + a GPU binary)

```bash
source benchmark/local_cuda_env.sh   # only needed if nsys/CUDA aren't already on PATH
python3 benchmark/perf_regression_track.py --capture
```

This resolves `nsys` (precedence: `$NSYS` -> `.local/nsight/.../nsys` -> `PATH`),
finds the newest `bin/alamo_gpu-2d-cuda*-g++` binary (override with
`--gpu-bin`), profiles 30 steps of `input` with the canonical parity overrides
baked in, runs `nsys stats --report cuda_api_sum --report cuda_gpu_kern_sum`,
parses the two CSVs, prints the 4 metrics, and appends a row to
`benchmark/perf_regression.csv` keyed by the current `git rev-parse HEAD`.

Useful flags:
- `--config NAME` — label for the row (default `wide_512_bf32_mgs128`, the
  current best box-sweep config). Use a different label if you capture a
  different grid/AMR config — rows are compared per-config, never across.
- `--steps N` — must match the `max_step` actually captured (default 30).
- `--wall-gpu-ms X --wall-cpu-ms Y` — attach wall-clock numbers from
  `baseline_suite.py` / `benchmark_gpu_cpu.sh` if you have them for the same run.
- `--notes "..."` — free-text annotation (e.g. which roadmap phase/commit motivated the capture).
- `--out-dir DIR` — where capture artifacts (`nsys.nsys-rep`, `stats_*.csv`,
  `capture.log`) land (default `benchmark/perf_regression_runs/<short-sha>_<timestamp>/`).

### Parse an existing nsys stats CSV pair (no GPU run needed)

```bash
python3 benchmark/perf_regression_track.py \
  --api-csv  analysis/results_phase22/stats_cuda_api_sum.csv \
  --kern-csv analysis/results_phase22/stats_cuda_gpu_kern_sum.csv \
  --steps 30 --config wide_512_bf32_mgs128 --notes "phase2.2 fused reduction"
```

Use this to backfill historical rows from prior captures
(`analysis/results_phase22/`, `analysis/results_phase22_nsub1/`,
`analysis/results_wide_grid/`) or to record a row on a machine without a GPU,
using a CSV captured elsewhere.

### Compare against the last recorded row for the same config

```bash
python3 benchmark/perf_regression_track.py --capture --compare
# or, dry-run without writing a new row:
python3 benchmark/perf_regression_track.py --capture --compare --no-record
```

`--compare` looks up the most recent *previous* row (excluding the current
git SHA) with a matching `config`, and flags any of `launches_per_step`,
`kernel_avg_us`, `sync_frac`, `wall_per_step_gpu_ms`, `wall_per_step_cpu_ms`
that got worse by more than `--threshold-pct` (default 10%). All of these
metrics are "higher is worse." If a regression is flagged the process exits 1
(otherwise 0) — suitable for a CI gate once a GPU runner exists.

## CSV schema (`benchmark/perf_regression.csv`)

```
date, git_sha, config, steps, launches_per_step, kernel_avg_us, sync_frac, wall_per_step_gpu_ms, wall_per_step_cpu_ms, notes
```

One row per capture. `git_sha` is the full `git rev-parse HEAD` at capture
time (the working tree may have uncommitted changes — `notes` is the place to
record that). `config` partitions history: never compare rows with different
`config` values to each other.

## Where rows live

`benchmark/perf_regression.csv` is generated by running this tool — it is not
created by this task and is not checked into the repository by default. Treat
it like `benchmark/baseline_runs/`: a local, growing log. If you want history
preserved across machines/sessions, commit it deliberately (it's plain CSV,
diffs cleanly) or copy it out before cleaning a worktree.

## CI-safety (no GPU / no nsys)

If `nvidia-smi` is not on `PATH` (no usable CUDA device), or `nsys` cannot be
resolved via `$NSYS` -> `.local/nsight/opt/nvidia/nsight-systems/*/target-linux-x64/nsys`
-> `PATH`, or no `bin/alamo_gpu-2d-cuda*-g++` binary exists, the tool prints

```
skipped: no GPU
```

and exits **0** — it never hard-fails a CPU-only run. If nsys/the GPU binary
exist but the actual profiling run fails (e.g. a real hardware/driver issue),
it prints `skipped: GPU capture unavailable/failed (...)`, cleans up the
partial capture directory, and still exits 0. This makes it safe to invoke
unconditionally from a CPU-only CI job; it only becomes a real regression gate
once it is invoked on a runner that actually has the GPU + binary + nsys.

Validated on this box: with `nvidia-smi` removed from `PATH` (no GPU
detected), `--capture` prints `skipped: no GPU` and writes nothing. With the
real CUDA device + nsys present, `--capture --compare` does a real ~30-step
GPU run, captures, parses, records a row, and compares it.

## Future CI wiring

This task intentionally does not touch `.github/workflows/**` (task 004 owns
CI). When a CUDA-capable self-hosted runner is available for `chamber-gpu`,
wire it as a *separate, gated* job (not blocking the CPU correctness gate):

```yaml
  perf-regression:
    runs-on: [self-hosted, cuda]   # placeholder label; must exist
    if: ...                        # gate as appropriate (e.g. manual dispatch, label)
    steps:
      - uses: actions/checkout@v4
      - name: Build GPU binary
        run: COMP=g++ CUDA_FP=fast SMOKE=0 ./benchmark/build_alamo_local_gpu.sh
      - name: Track perf
        run: python3 benchmark/perf_regression_track.py --capture --compare --threshold-pct 10
```

On a CPU-only runner the same step is harmless (prints `skipped: no GPU`,
exit 0), so it is also safe to add unconditionally to the existing
correctness workflow now and let it become load-bearing once a GPU runner is
attached.
