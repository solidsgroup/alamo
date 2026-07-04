# Result: Task 004 — Phase 5.1 correctness CI gate for chamber-gpu

## Status
DONE

## What changed
Created exactly the two files allowed by the task, both new/untracked:

1. `benchmark/ci_golden_compare.sh` — the correctness harness invoked by CI
   (and runnable locally). Mode-switched via `GOLDEN_MODE=cpu|gpu`
   (default `cpu`):
   - **cpu leg** (runs on any standard runner, no GPU needed): configures
     `./configure --dim=2 --comp=g++`, builds, locates the freshly built
     `bin/alamo-2d-*` binary, exports it as `CPU_BIN`, then runs
     `python3 benchmark/baseline_suite.py check --profiles=cpu` — this fails
     the build (nonzero exit) on any real mismatch against
     `benchmark/baseline_references/{canonical_step1,canonical_step2,
     eta_expression_step1}/cpu.json`, since `baseline_suite.py check` itself
     returns 1 on any FAIL and the script has no `continue-on-error` around
     that step.
   - **gpu leg** (only meant to run on a CUDA-capable runner): requires
     `nvcc` on `PATH` (hard-exits 2 otherwise), detects `ARCH` via
     `nvidia-smi --query-gpu=compute_cap`, configures with
     `--cuda "${ARCH}" --cuda-fp strict` (the existing no-fast-math
     correctness build mode already in `configure`), builds, locates the
     `alamo_gpu-2d-nofast-cuda${ARCH}-*` binary, and runs
     `baseline_suite.py check --profiles=gpu_strict` against the
     `gpu_strict.json` references.
   - **NaN-flag assertion smoke** (both legs, after the golden compare): runs
     a short 2-step Flame simulation (`max_step=2`, no plotting) on the
     just-built binary. `Flame.cpp`'s advance kernel calls
     `Util::SetDeviceError(advance_error_flag)` on any NaN/Inf in
     K/rho/cp/L/eta/alpha/mdot/heatflux/temperature, and
     `Util::AbortIfDeviceError` (`src/Util/Util.H:148`) calls `Util::Abort`
     — which terminates the process — if that flag was ever set. A clean
     exit 0 from `mpiexec` is therefore a valid proxy that the tripwire never
     fired; a nonzero exit fails the script and dumps the tail of the run log.
   - Validated with `bash -n benchmark/ci_golden_compare.sh` (passes) and made
     executable (`chmod +x`).

2. `.github/workflows/chamber-gpu-correctness.yml` — triggers on
   `push: branches: [chamber-gpu]` and `workflow_dispatch`. Two jobs:
   - `golden-cpu`: `runs-on: ubuntu-24.04`, installs deps via the existing
     `.github/workflows/dependencies-ubuntu-24.04.sh` (same script `linux.yml`
     uses), then `GOLDEN_MODE=cpu bash benchmark/ci_golden_compare.sh` with
     **no `continue-on-error`** on that step (only the log-upload step uses
     `if: always()`, which doesn't suppress failure). Uploads
     `benchmark/baseline_runs` and `benchmark/ci_nan_smoke_cpu` as artifacts.
   - `golden-gpu`: `needs: golden-cpu`, gated by
     `if: contains(join(github.event.repository.topics, ','), 'has-cuda-runner')`
     and `runs-on: [self-hosted, cuda]` — i.e. it only attempts to schedule on
     a runner advertising a `cuda` label, and is additionally gated by a repo
     topic flag so it stays inert (skipped, not failed) on repos/forks with no
     CUDA runner, exactly as the task's invariant requires ("no GPU CI runner
     may exist"). When a real CUDA runner exists, flipping the repo topic (or
     simplifying the `if:` to just the `runs-on` label match) brings it live
     with zero script changes. `continue-on-error: false` is explicit on the
     job to make clear this leg is meant to gate, not just report, once live.
   - Verified the YAML parses with `python3 -c "import yaml; yaml.safe_load(...)"`
     → `YAML OK`, and that the referenced
     `.github/workflows/dependencies-ubuntu-24.04.sh` exists.

## baseline_suite.py CLI — no adaptation needed
Read the real CLI before writing anything. Confirmed by running
`python3 benchmark/baseline_suite.py list` locally:
- `check` mode (used here) takes `--case` (repeatable, defaults to all 3
  cases) and `--profiles` (comma list, default `cpu,gpu_fast,gpu_strict`).
  For each selected case/profile it requires
  `benchmark/baseline_references/<case>/<profile>.json` to exist (raises
  `SystemExit` otherwise — all three already exist for all three profiles),
  re-runs the case fresh via `run_case` (which shells out to `mpiexec`,
  requires the binary resolved through `default_binaries()`/`CPU_BIN`,
  `GPU_FAST_BIN`, `GPU_STRICT_BIN` env overrides), reads the fresh
  `thermo.dat`, and compares column-wise against the stored reference with
  `abs_tol`/`rel_tol` per case; returns 1 if any case/profile fails, 0
  otherwise.
- This matched the task's assumptions exactly (the task doc already named
  `check` and the reference layout correctly) — the only "adaptation" was
  discovering and using the `CPU_BIN`/`GPU_STRICT_BIN` env-var override and
  `--profiles=<single>` (to scope `check` to just the leg under test, since
  the default runs all three profiles and would otherwise require all three
  binaries to be present).

## Files touched
- `benchmark/ci_golden_compare.sh` (new, executable)
- `.github/workflows/chamber-gpu-correctness.yml` (new)

## Files explicitly NOT touched
- `.github/workflows/linux.yml`, `.github/workflows/performance.yml`,
  `benchmark/baseline_suite.py`, `Makefile`, `configure`, any source file —
  all read-only for reference/discovery. No commit was made.

## Verification
- `bash -n benchmark/ci_golden_compare.sh` → passes (no output, syntax clean).
- `python3 -c "import yaml; yaml.safe_load(open('.github/workflows/chamber-gpu-correctness.yml'))"`
  → `YAML OK`.
- Confirmed `.github/workflows/dependencies-ubuntu-24.04.sh` exists (the new
  workflow's only external reference besides standard actions).
- Confirmed `--profiles=cpu` / `--profiles=gpu_strict` are accepted values by
  running `baseline_suite.py list` against the live repo state (binaries for
  all three profiles already exist locally from prior GPU-port work, e.g.
  `bin/alamo_gpu-2d-nofast-cuda86-g++`).
- Did not run a full build or the script end-to-end, per task instructions
  ("do not run a full build").

## Open follow-ups / risks
- The `golden-gpu` job's gate (`repository.topics` containing
  `has-cuda-runner`) is a placeholder convention, not a real label that
  exists anywhere yet — when/if a real CUDA self-hosted runner (e.g. an
  upgraded `scooter`) is provisioned, someone should either add that repo
  topic or simplify the `if:` condition; either way no script changes are
  needed.
- The NaN-flag smoke test is a necessary-but-not-sufficient correctness
  check: it only proves the tripwire didn't fire on the specific short
  `input` case run with `max_step=2`; it does not prove the tripwire would
  fire correctly if NaNs were actually injected (no fault-injection test
  exists). Out of scope per the task's non-goals.
