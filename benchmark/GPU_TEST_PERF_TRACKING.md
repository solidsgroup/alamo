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
  (was 5P/4F; see `benchmark/GPU_TEST_SUITE_FIXES.md`).
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

## Iteration template (copy for the next fix-set)

```
## Iteration N — YYYY-MM-DD — <what changed>

- Result: <P>/<F>/<S>. Code: <branch>@<hash> (+ uncommitted? y/n).
- HW: <gpu>.
- Notable deltas vs previous iteration: <P1 ms/step X→Y, etc.>

| Test | Binary | Config | Steps | Wall (s) | ms/step | Kind |
| ... |
```
