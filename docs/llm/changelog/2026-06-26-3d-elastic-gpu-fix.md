# 3D Elastic GPU Fix — `F.inverse().transpose()` device fault — 2026-06-26

Branch `chamber-gpu`. Repo `/home/jackplum/Projects/alamo`. Context: exercising the
**3D combined flame+elastic** path on the local A1000 (toward roadmap-v2 task A2)
surfaced a hard GPU crash. Root-caused, fixed, verified.

---

## TL;DR

3D elastic aborts on the **first elastic solve** with `CUDA error 719: unspecified
launch failure`. Root cause is a **single chained Eigen expression**:

```cpp
Eigen::Matrix3d FinvT = F.inverse().transpose();   // nested Transpose<Inverse<Matrix3>>
```

This nested expression template **faults on the GPU** (not on CPU, not in 2D).
**Fix:** materialize the inverse before transposing:

```cpp
Eigen::Matrix3d Finv = F.inverse();
Eigen::Matrix3d FinvT = Finv.transpose();
```

Two sites, both in `src/Model/Solid/Finite/NeoHookean.H`, both inside
`#elif AMREX_SPACEDIM==3` blocks: **line 75 (`DW`)** and **line 127 (`DDW`)**.
Committed as comments-in-place explaining why.

---

## Symptom

- `bin/alamo_gpu-3d-*-cuda86-g++` on any elastic input (e.g. `input_3d_flame_128_elastic`,
  even shrunk to `n_cell=8 8 8 max_level=0`) aborts immediately on `TimeStepBegin`'s
  elastic solve: `ELASTIC SOLVE: t=0 step=0` → `CUDA error 719 ... unspecified launch
  failure` at `AMReX_GpuDevice.H:444`.
- compute-sanitizer memcheck attributes the launch to `Newton.H:196`
  (`Solver::Nonlocal::Newton<NeoHookeanPredeformed>::prepareForSolve`) but reports
  **no illegal-access / no assert** — only the launch failure surfacing at the next
  `cudaLaunchKernel`. (That "719 with no memcheck annotation" signature is what a
  device fault inside an inlined expression-template evaluation looks like.)

## What it is NOT (ruled out, with evidence)

- **Not a device stack overflow.** `ALAMO_GPU_STACK=8192` and `65536` both fail
  identically; the kernel's own frame is `STACK:0` (cuobjdump `-res-usage`). The
  diagnostic probe + comment in `src/alamo_gpu.cc` (which hypothesized a stack
  overflow in `prepareForSolve`) is **wrong** and can be removed.
- **Not launch resources / local memory.** `prepareForSolve` is `REG:72 LOCAL:0`.
- **Not an out-of-bounds access.** compute-sanitizer memcheck finds none. The
  `Matrix4<3,Major>::operator()` index map covers all 81 `uid` values (verified).
- **Not Eigen-on-device in general.** Eigen 3.4's size-3 `compute_inverse` IS
  `EIGEN_DEVICE_FUNC`; `F.inverse()`, `.determinant()`, `(F*Fᵀ).trace()`, `std::pow`
  each run clean on device in isolation. The build uses `--Werror
  cross-execution-space-call`, so a host-only-on-device call would have been a
  *compile* error, not a runtime fault.
- **Not config / physics.** A 3D **CPU** build (`bin/alamo-3d-g++`, built this
  session) runs the identical input to completion.

## How it was isolated

1. Reproduced on the real binary; localized to `prepareForSolve` (`Newton.H:197`),
   model `NeoHookeanPredeformed`, which calls `DW`/`DDW`.
2. Proved GPU-specific via the 3D CPU build (same input → clean).
3. Replicated `NeoHookean::DDW` 3D math in a standalone `.cu` compiled with the exact
   nvcc flags (`sm_86`, `--expt-relaxed-constexpr --expt-extended-lambda --fmad=false
   -maxrregcount=255 -DALAMO_GPU`) → reproduced `rc=719`.
4. Bisected the math in isolated processes (719 is sticky, so each op in a fresh
   process):
   - `F.inverse()` → rc=0
   - `F.inverse().transpose()` → **rc=719**
   - `Fi = F.inverse(); Ft = Fi.transpose()` → rc=0
   - `determinant`, `(F*Fᵀ).trace()`, `pow` → rc=0
5. Verified the split-statement fix on the full DDW replica (rc=0) **and** end-to-end:
   rebuilt `bin/alamo_gpu-3d-cuda86-g++`, ran 3 steps of `input_3d_flame_128_elastic`
   — three elastic solves complete, no fault.

## Why 2D never hit it

The 2D branches of `NeoHookean::DW`/`DDW` and `NeoHookeanPredeformed::RelativeF`
**hand-roll every inverse** (explicit `FinvT(0,0)=...` etc.) and never form an
`.inverse().<chained-op>` expression. Only the 3D branches use Eigen `.inverse()`,
and only `.inverse().transpose()` (chained) is the faulting form. The 3D GPU
crossover runs to date all used `elastic.type=disable`, so this device path had
**never been exercised** before — consistent with roadmap-v2's "combined 3D
flame+elastic is unmeasured."

## Blast radius

`grep -rn "\.inverse()\." src --include=*.H --include=*.cpp` → exactly the two
`NeoHookean.H` sites. No other chained-inverse-expression in a device path today.
Recommend keeping this grep as a standing check (see roadmap-v2 A4).

## Build / artifact notes

- Fix is **inside `#if AMREX_SPACEDIM==3` only** → 2D binaries are behavior-identical
  and were NOT rebuilt.
- Rebuilt `bin/alamo_gpu-3d-cuda86-g++` (fast) **with the fix**. Built a new
  `bin/alamo-3d-g++` (CPU) for the cross-check. The 3D-nofast/strict GPU binary was
  **not** yet rebuilt with the fix.
- The freshly-built `alamo_gpu-3d-cuda86-g++` has **no embedded RUNPATH** (unlike the
  older binaries) → it needs `LD_LIBRARY_PATH=.local/cuda-12.6.3-redist/lib` at
  runtime. nvcc must be on `PATH` to build: `export PATH=$PWD/.local/cuda/bin:$PATH`.
- This is the **same bug class** as the elixir UAF (`Operator.cpp`): a GPU-only
  failure that "works on CPU." See [[gpu-elastic-fixed]].

## Related sibling bug spotted (not fixed)

`Set::Field<NeoHookeanPredeformed>::NComp()` returns `2 + DIM*DIM` (= **11 in 3D**) but
`Copy()`/`Name()` only handle 6 components in both 2D and 3D (the 3D arm still has
`//Util::Abort(INFO, "Not implemented")`). This is **plotfile-output** only (does not
affect the solve — CPU 3D runs fine) but means 3D model-field plot components 6–10 are
unwritten/garbage. Deferred.

## GPU test suite snapshot — ✅ RESOLVED 2026-06-26

> **UPDATE:** the deferred analysis below is **done** →
> `changelog/2026-06-26-gpu-test-suite-fixes.md` +
> `benchmark/GPU_TEST_SUITE_FIXES.md`. All 4 failures root-caused and fixed
> (**9 passed / 0 failed**). They were four distinct bugs (3 stale-input decks +
> the C2 checkpoint deck); two **real source defects** surfaced en route — an
> `Integrator::Restart` node-fab out-of-bounds segfault and a headerless
> restart-`thermo.dat` bug (both `src/Integrator/Integrator.cpp`). The F2
> "stale input" prediction was correct, as was the "C1/C2/C3 deck issues are
> pre-existing/harness" caveat. Original deferred notes kept below for provenance.

## GPU test suite snapshot — ANALYSIS DEFERRED (original notes, superseded)

Ran `python3 tests/GPU/run_gpu_tests.py` after the fix (env:
`LD_LIBRARY_PATH=.local/cuda-12.6.3-redist/lib`, `ALAMO_GPU_STACK=8192`).
**Result: 5 passed, 4 failed, 0 skipped.**

| Test | Result | Note |
| --- | --- | --- |
| C1_correctness_elastic | FAIL (2.2s) | analyze later |
| C2_restart_roundtrip | FAIL (11.3s) | analyze later |
| C3_multibox_elastic_stress | FAIL (10.0s) | analyze later |
| C4_amr_correctness | PASS (2.6s) | |
| F1_smoke_flame_only | PASS (0.9s) | |
| F2_smoke_elastic | FAIL (9.9s) | **known: stale input** — deck uses old `propellant.homogenize.m_ap/E_ap/...` params; current parser requires `rho_prop` etc. → `query_required ... rho_prop missing`. Test-authoring bug, unrelated to the elastic fix. |
| P1_perf_2d_hiRes_noAMR | PASS (14.9s) | |
| P2_perf_2d_hiRes_AMR3 | PASS (300.1s) | hit perf-test 300s cap (counts as pass) |
| P3_perf_3d_256 | PASS (300.1s) | **3D, ran on the rebuilt fixed binary**; hit 300s cap |

**These 4 failures are NOT yet analyzed — deliberately deferred per user.** Caveats for
whoever picks this up:
- The failures are almost certainly **pre-existing / test-harness issues, not the
  elastic fix** — they're 2D tests using 2D binaries that the fix does not touch
  (fix is `#if AMREX_SPACEDIM==3` only), and P3 (the only test exercising the rebuilt
  3D binary) passed.
- F2 = confirmed stale input (above).
- C1/C2/C3 reasons were **not captured** (the driver buffered per-test stdout behind a
  `tail`, and a direct re-run reported SKIP because `find_binary`/env differed outside
  the driver). First step next time: run each with the driver's env and capture full
  stdout, e.g.
  `ALAMO_GPU_STRICT_BIN=$PWD/bin/alamo_gpu-2d-nofast-cuda86-g++ ALAMO_CPU_BIN=$PWD/bin/alamo-2d-g++ python3 tests/GPU/<T>/test.py <outdir>`.
- The whole `tests/GPU/` suite is itself recent (`docs/agent_plans/20260625-gpu-tests/`)
  and may have its own input-staleness across C1–C3 like F2 did.
