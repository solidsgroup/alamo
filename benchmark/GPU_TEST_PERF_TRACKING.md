# GPU Test Suite — Performance Tracking

Per-test wall-time / throughput for `tests/GPU/` (`python3 tests/GPU/run_gpu_tests.py`),
recorded **one iteration block per fix-set** so suite-level perf regressions are
visible across changes. This complements `benchmark/PERF_TRACKING.md` (which
tracks per-commit *kernel* counters for the canonical Flame case); this file
tracks the *test suite's* end-to-end timings.

## What's tracked & why

- **P1/P2/P3 (`ms/step`)** are the real throughput metrics — these run to a fixed
  `max_step=10000` (or the 300 s per-test cap) and report `elapsed*1000/(rows-1)`.
  Compare `ms/step` across iterations; a jump is a throughput regression.
- **F*/C* tests** are pass/fail correctness/smoke tests; their wall time is a
  coarse health signal (e.g. an elastic-solver conditioning change showing up as
  a big C1/C3 slowdown), not a precision metric — they're short and launch-bound.
- P2/P3 normally hit the **300 s cap** (counts as PASS); track `steps completed`
  in the cap window as their throughput proxy (more steps in 300 s = faster).

## How to run (reproduce a row)

```bash
cd /home/jackplum/Projects/alamo
export LD_LIBRARY_PATH=.local/cuda-12.6.3-redist/lib:$LD_LIBRARY_PATH
export ALAMO_GPU_STACK=8192
export ALAMO_GPU_STRICT_BIN=$PWD/bin/alamo_gpu-2d-nofast-cuda86-g++
export ALAMO_GPU_BIN=$PWD/bin/alamo_gpu-2d-cuda86-g++
export ALAMO_GPU_3D_BIN=$PWD/bin/alamo_gpu-3d-cuda86-g++
export ALAMO_CPU_BIN=$PWD/bin/alamo-2d-g++
python3 tests/GPU/run_gpu_tests.py
```

Binary per test: F1/P1/P2 → `alamo_gpu-2d-cuda86-g++` (fast); P3 →
`alamo_gpu-3d-cuda86-g++` (fast); F2/C1/C2/C3 → `alamo_gpu-2d-nofast-cuda86-g++`
(strict); C1/C4 also run `alamo-2d-g++` (CPU) for parity.

---

## Iteration 1 — 2026-06-26 — post-fix baseline (FIRST all-green run)

- **Result: 9 passed / 0 failed / 0 skipped.** First time the suite is green
  (was 5P/4F; see `benchmark/archive/GPU_TEST_SUITE_FIXES.md`).
- **HW:** NVIDIA RTX A1000 (8188 MiB, ~50 W shared desktop GPU — absolute times
  are not exclusive-GPU numbers; use them only for relative cross-iteration
  comparison on this machine).
- **Code:** `chamber-gpu` @ `7f5095c3e` + uncommitted working-tree fixes
  (Integrator.cpp restart OOB + thermo header; F2/C1/C3 deck rewrites;
  C1/C2 test reworks). 2D CPU + GPU-strict binaries rebuilt with the source fixes.

| Test | Binary | Config | Steps | Wall (s) | ms/step | Kind |
|------|--------|--------|------:|---------:|--------:|------|
| F1_smoke_flame_only        | gpu fast   | 2D flame, ml=1, 8³ base, max_step=20      |   19 |   1.0 |    —    | smoke |
| F2_smoke_elastic           | gpu strict | 2D flame+elastic, ml=1, mgs=64, 50 steps  |   49 |   1.6 |  ~32    | smoke |
| C1_correctness_elastic     | cpu+gpu str| 2D flame+elastic, ml=2, 30 steps, CPU↔GPU |   30 |   2.6 |    —    | correctness |
| C2_restart_roundtrip       | gpu strict | 2D, ml=0, 2×20 steps (write+restart)      |   20 |   1.9 |    —    | correctness |
| C3_multibox_elastic_stress | gpu strict | 2D flame+elastic, ml=1, mgs=32, 200 steps |  199 |   3.2 |  ~16    | stress (39 elastic solves) |
| C4_amr_correctness         | cpu+gpu str| 2D AMR, 50 steps, CPU↔GPU                  |   50 |   3.0 |    —    | correctness |
| P1_perf_2d_hiRes_noAMR     | gpu fast   | 2D 256², ml=0, mgs=128, max_step=10000    | 9999 |  15.2 |  **1.52** | throughput (completed) |
| P2_perf_2d_hiRes_AMR3      | gpu fast   | 2D 64², ml=3, mgs=128, 300 s cap          | 4733 | 300.1 | **63.40** | throughput (capped) |
| P3_perf_3d_256             | gpu fast   | 3D 256²×128, ml=1, mgs=128, 300 s cap     |  236 | 300.1 | **1271.54** | throughput (capped) |

**Headline throughput numbers to watch (RTX A1000):**
- P1 (256² flat, no AMR): **1.52 ms/step** — the cleanest single-grid throughput metric.
- P2 (64² + 3 AMR levels): **63.40 ms/step** — AMR/regrid + subcycling overhead dominated; 4733 steps in the 300 s window.
- P3 (256²×128 3D): **1271.54 ms/step** — 236 steps in the 300 s window; bounded by the 8 GB card.

Notes:
- F2/C3 `ms/step` are approximate (`wall / steps`, includes init + elastic solves
  every 5 steps); they're smoke/stress, not throughput benchmarks.
- C3 does 39 elastic MLMG solves in 3.2 s on the small grid — the casing-stiffness
  softening (Bug 3) is what keeps these converging.

---

## Iteration 2 — 2026-06-28 — local phase-field structural speedup checkpoint

- **Scope:** phase-field / Flame path only; elastic solver intentionally untouched.
- **Branch/worktree:** `codex/gpu-pf-structural-speedups` in
  `/home/jackplum/Projects/alamo-pf-gpu-opt`.
- **HW:** local NVIDIA RTX A1000, sm_86, 8188 MiB. Use for A/B directionality and
  launch/API deltas; leave NOVA for major version-scale A100 benchmarks.
- **Harness:** `benchmark/local_pf_gpu_ab.sh` compares explicit `BASE_BIN` and
  `OPT_BIN`, runs local wall-clock repeats, and can collect Nsight Systems
  `cuda_api_sum` / `cuda_gpu_kern_sum` reports.

### Local A/B signal

Input decks are the existing phase-field-only GPU perf cases with
`elastic.type=disable`.

| Case | Steps | Baseline median | Optimized median | Delta |
|------|------:|----------------:|-----------------:|------:|
| P1 no-AMR thermal-on | 60 | 0.970 s | 0.940 s | 1.03x / 3.1% faster |
| P2 AMR3 thermal-on | 30 | 2.010 s | 1.470 s | **1.37x / 26.9% faster** |

Nsight Systems on P2 AMR3, 30 steps:

| Metric | Baseline | Optimized | Delta |
|--------|---------:|----------:|------:|
| `cudaLaunchKernel` calls | 108,035 | 38,555 | **64.3% fewer** |
| CUDA API total time | 1607.5 ms | 832.2 ms | **48.2% lower** |
| memcpy/memset API calls | 29,644 | 13,804 | **53.4% fewer** |
| GPU kernel aggregate time | 483.8 ms | 522.1 ms | noise / not the win source |

TinyProfiler attribution on P2 AMR3:

| Region | Baseline | Optimized | Read |
|--------|---------:|----------:|------|
| `Integrator::FillPatch` calls | 6300 | 1350 | static/non-evolving field registration removed repeated generic state fill work |
| `FillPatchTwoLevels` inclusive | 0.916 s | 0.208 s | main AMR structural win |
| `FillPatchSingleLevel` calls | 9180 | 2070 | launch count follows field-count reduction |
| `FabArray::FillBoundary()` calls | 15710 | 3650 | BC/fill launch churn reduced |
| `Integrador::Flame::Advance` inclusive | 0.299 s | 0.222 s | phase/thermal kernel path modestly faster |

Validation:

| Test | Result | Notes |
|------|--------|-------|
| F1 smoke flame-only | PASS | profile CUDA binary, 19 steps |
| C4 AMR correctness | PASS | CPU vs GPU contour compare, `eta=0.0000`, `temp=0.0000`; used the available profile CUDA binary as the strict override in this worktree |
| Harness smoke | PASS | `REPEATS=1 RUN_NSYS=0 P2_STEPS=5 P1_STEPS=5` |

Notes:
- The speedup mechanism is launch/API/fill reduction, not a raw stencil kernel
  throughput improvement.
- Forcing `thermal.on=0` onto the thermal perf input aborts in both the clean
  baseline and optimized binaries, so it is not used as an A/B signal.
- Raw local Nsight reports are intentionally ignored by git under
  `benchmark/local_ab_*/`.

---

## Iteration 3 — 2026-06-28 — skip non-evolving AMR average-down

- **Scope:** generic AMR synchronization change motivated by the phase-field
  Flame deck. `Integrator::TimeStep` now averages down registered cell/node
  fields only when their existing `evolving` flag is true, matching the behavior
  that base fields already used.
- **Why:** after fine-level subcycling, AMReX was still averaging down static
  and diagnostic fields such as `phi`, `eta_old`, `eta_0`, `L`, `mdot`,
  `alpha`, `heatflux`, and `laser`. Those fields were already removed from
  per-step FillPatch in Iteration 2, but they still paid AMR average-down launch
  cost.

### Local A/B signal

Fresh run with the same baseline binary and optimized worktree after the
average-down change:

| Case | Steps | Baseline median | Optimized median | Delta |
|------|------:|----------------:|-----------------:|------:|
| P1 no-AMR thermal-on | 60 | 0.840 s | 0.810 s | 1.04x / 3.6% faster |
| P2 AMR3 thermal-on | 30 | 1.650 s | 1.140 s | **1.45x / 30.9% faster** |

Nsight Systems on P2 AMR3, 30 steps:

| Metric | Baseline | Optimized | Delta |
|--------|---------:|----------:|------:|
| `cudaLaunchKernel` calls | 108,035 | 31,835 | **70.5% fewer** |
| CUDA API total time | 1205.3 ms | 743.4 ms | **38.3% lower** |
| memcpy/memset API calls | 29,644 | 8,042 | **72.9% fewer** |
| `amrex::average_down_w_geom` calls | 2769 | 669 | **75.8% fewer** |
| `amrex::average_down_w_geom` inclusive | 0.139 s | 0.033 s | **76.1% lower** |

Validation:

| Test | Result | Notes |
|------|--------|-------|
| CUDA profile rebuild | PASS | `make -j8 bin/alamo_gpu` |
| F1 smoke flame-only | PASS | profile CUDA binary, 19 steps |
| C4 AMR correctness | PASS | CPU vs GPU contour compare, `eta=0.0000`, `temp=0.0000`; same profile-CUDA strict override caveat |

---

## Local phase-field GPU perf gate

Use the workstation gate for structural phase-field changes before spending NOVA
time. It is intentionally local and phase-field-only; NOVA remains the venue for
major version-scale A100 benchmarks and large 3D confirmation.

```bash
cd /home/jackplum/Projects/alamo-pf-gpu-opt
BASE_BIN=/home/jackplum/Projects/alamo-pf-gpu-base/bin/alamo_gpu-2d-profile-cuda86-g++ \
OPT_BIN=/home/jackplum/Projects/alamo-pf-gpu-opt/bin/alamo_gpu-2d-profile-cuda86-g++ \
REPEATS=3 RUN_NSYS=1 P2_STEPS=30 P1_STEPS=60 \
PERF_GATE=1 GATE_MIN_P2_SPEEDUP=1.10 \
benchmark/local_pf_gpu_ab.sh
```

Artifacts:
- `wall_repeats.csv` has raw timings.
- `summary.csv` has median/mean/stdev by case and label.
- `speedups.csv` has median baseline/optimized speedups.
- `nsys_p2_{base,opt}/` contains CUDA API/kernel summaries when `RUN_NSYS=1`.

Gate behavior:
- Default mode (`PERF_GATE=0`) is exploratory and never fails on performance.
- Gate mode (`PERF_GATE=1`) fails if `p2_amr_thermal_on` speedup is below
  `GATE_MIN_P2_SPEEDUP`.
- Smoke-verified locally on 2026-06-28 with `REPEATS=1 RUN_NSYS=0 P2_STEPS=5
  P1_STEPS=5 PERF_GATE=1 GATE_MIN_P2_SPEEDUP=1.05`; observed P2 speedup 1.289x.

---

## Iteration 4 — 2026-06-28 — quiet chamber logs and skip disabled thermo diagnostics

- **Scope:** Flame/chamber host-side wall-clock cleanup plus disabled-output
  diagnostic gating.
- **Changes:**
  - Added `chamber.verbose` (default `0`) and gated the four per-step chamber
    pressure/mass/volume/dpdt messages behind it.
  - `Flame::TimeStepComplete` now computes max/min thermo diagnostics only when
    a thermo row will be written on the next step. Profiling/perf runs that pass
    `amr.thermo.plot_int=-1` no longer pay those reductions.
- **Why:** long phase-field perf runs were writing four chamber messages every
  step, and no-output A/B runs still computed non-extensive thermo diagnostics
  that had no consumer.

### Local A/B signal

Fresh local run after Iteration 4:

| Case | Steps | Baseline median | Optimized median | Delta |
|------|------:|----------------:|-----------------:|------:|
| P1 no-AMR thermal-on | 60 | 0.980 s | 0.940 s | 1.04x / 4.1% faster |
| P2 AMR3 thermal-on | 30 | 2.030 s | 1.330 s | **1.53x / 34.5% faster** |

Nsight Systems on P2 AMR3, 30 steps:

| Metric | Baseline | Optimized | Delta |
|--------|---------:|----------:|------:|
| `cudaLaunchKernel` calls | 108,035 | 31,505 | **70.8% fewer** |
| CUDA API total time | 1244.4 ms | 736.6 ms | **40.8% lower** |
| memcpy/memset API calls | 29,644 | 8,042 | **72.9% fewer** |
| chamber log lines | 120 | 0 | removed by default |

Validation:

| Test | Result | Notes |
|------|--------|-------|
| CUDA profile rebuild | PASS | `make -j8 bin/alamo_gpu` |
| F1 smoke flame-only | PASS | normal thermo output still enabled |
| C4 AMR correctness | PASS | CPU vs GPU contour compare, `eta=0.0000`, `temp=0.0000` |

Notes:
- Set `chamber.verbose=1` to restore the old per-step chamber messages.
- The thermo-diagnostic gate affects only the non-extensive diagnostic reductions
  in `TimeStepComplete`; chamber integrated variables are still computed before
  advance as required by the variable-pressure model.

---

## Iteration 5 — 2026-06-28 — local AMR-shape sweep and rejected code probes

- **Scope:** local phase-field-only profiling after Iteration 4. No solver-code
  checkpoint was kept in this iteration.
- **Why:** the remaining optimized P2 profile is still launch/sync dominated:
  `cudaLaunchKernel` and `cudaStreamSynchronize` account for most CUDA API time,
  while AMR FillPatch/FillPatchInterp remains the largest non-Advance region.

### Local AMR-shape signal

Same optimized binary as Iteration 4, P2 AMR thermal-on input, 30 steps,
3 wall-clock repeats on the local RTX A1000:

| Variant | Median wall | Speedup vs current P2 default |
|---------|------------:|------------------------------:|
| Current P2 default (`max_level=3`, default subcycling) | 1.320 s | 1.00x |
| `amr.nsubsteps=1` | 0.880 s | **1.50x** |
| `amr.max_level=1 amr.n_cell="256 256 4" amr.nsubsteps=1` | 0.830 s | **1.59x** |
| Uniform fine, `amr.max_level=0 amr.n_cell="512 512 4"` | 0.860 s | **1.54x** |

Nsight Systems on the default P2 shape vs non-subcycling AMR:

| Metric | Default | `amr.nsubsteps=1` | Delta |
|--------|--------:|------------------:|------:|
| `cudaLaunchKernel` calls | 31,505 | 8,585 | **72.8% fewer** |
| CUDA API total time | 752.2 ms | 399.3 ms | **46.9% lower** |
| memcpy/memset API calls | 8,042 | 2,042 | **74.6% fewer** |
| `cudaStreamSynchronize` calls | 38,010 | 9,840 | **74.1% fewer** |
| GPU kernel aggregate time | 495.2 ms | 125.0 ms | **74.8% lower** |

Read: the v2 roadmap AMR hypothesis holds locally. Deep subcycling AMR is still
the dominant remaining structural wall-clock lever for the phase-field GPU path.
This should become an accuracy-checked B2 recommendation before changing
production input defaults.

### Rejected probes

| Probe | Result | Decision |
|-------|--------|----------|
| Replace Flame `MFIter(..., true)` with `amrex::TilingIfNotGPU()` | Direct previous-vs-current A/B: P2 1.360 s → 1.350 s, P1 flat; launch/mem/sync counts unchanged | Reverted; no structural wall-clock gain |
| Mark `temps_mf` non-evolving | C4 failed (`eta` max rel error 0.1937, `temp` 0.2879), direct A/B was neutral/slower | Reverted; `temps_mf` participates in AMR state consistency despite being zero-ghost/write-disabled |

Validation after reverting rejected probes:

| Test | Result | Notes |
|------|--------|-------|
| CUDA profile rebuild | PASS | `make -j8` |
| C4 AMR correctness | PASS | CPU vs GPU contour compare, `eta=0.0000`, `temp=0.0000` |

Artifacts are local-only and ignored under `benchmark/local_amr_shape_*/`.

---

## Iteration template (copy for the next fix-set)

```
## Iteration N — YYYY-MM-DD — <what changed>

- Result: <P>/<F>/<S>. Code: <branch>@<hash> (+ uncommitted? y/n).
- HW: <gpu>.
- Notable deltas vs previous iteration: <P1 ms/step X→Y, etc.>

| Test | Binary | Config | Steps | Wall (s) | ms/step | Kind |
| ... |
```
