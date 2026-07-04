# Phase 3.1 — 3D GPU Readiness Audit

**Date:** 2026-06-20  **Branch:** `chamber-gpu`  **GPU:** RTX A1000 (sm_86, 8 GB)
**Roadmap:** `GPU-OPT-ROADMAP.txt` Phase 3, step 3.1.

## Verdict: 3D GPU phase-field path is READY ✅ (elastic remains CPU-resident per D1)

The 3D CUDA build compiles and a 3D Flame step runs to completion on the device,
exercising the GPU-resident phase-field + thermal path. The device elastic solve is
NOT ready (CUDA-719, consistent with D1) and is disabled for GPU runs.

## What was done

1. **Build.** `DIM=3 PROFILE=1 SMOKE=0 ARCH=86 ./benchmark/build_alamo_local_gpu.sh`
   → `bin/alamo_gpu-3d-profile-cuda86-g++`. Clean compile (exit 0), 1772 device entry
   functions generated for sm_86. No 3D-specific source changes were required — the
   linear-algebra layer (`Set::Matrix4<3>`, `Matrix4_Isotropic` full 3×3 stress,
   quaternion ops) is already dimension-guarded and 3D-clean.
2. **Smoke run.** `mpiexec -np 1 ./bin/alamo_gpu-3d-profile-cuda86-g++ input_3d_smoke
   max_step=1 stop_time=1e-12_s` → **exit 0**, `STEP 1 ends. TIME = 0.0001`,
   `chamber.mdot = 0_kg/s`. The GPU phase-field + thermal + chamber path runs in 3D.

## Findings

### F1 — IC blocker (resolved). `src/IC/BMP.H` is 2D-only.
The production `input` initializes `pf.eta`/`phi` from `.bmp` images, which silently
collapse z in 3D. Resolution: 3D inputs use `IC::Expression` analytic ICs (vars x,y,z,t).
See `input_3d_smoke` for the proven template.

### F2 — Elastic trap (resolved, important). `elastic.on = 0` does NOT disable the solve.
`Mechanics::TimeStepBegin` (`src/Integrator/Base/Mechanics.H:186-207`) runs the Static
elastic solve regardless of `elastic.on`, and on device it aborts with
**CUDA error 719 (unspecified launch failure)** at the first solve — the same device
elastic failure recorded in D1. The correct switch is **`elastic.type = disable`**
(`Mechanics::Type::Disable`, Mechanics.H:28,58), which early-returns from every elastic
path. With `disable`, the `model_*` and `elastic.bc.*` blocks must be OMITTED or ALAMO's
strict ParmParse aborts on the unused entries.

### F3 — 3D BCs. eta and temp BCs must declare all six faces (xlo,xhi,ylo,yhi,**zlo,zhi**).
The 2D inputs only had four; the smoke input adds zlo/zhi.

### F4 — Register pressure / occupancy (preliminary, from ptxas).
The hottest kernels use up to **~255 registers/thread** (several at 230–255). At 255
regs, occupancy is capped well below peak on sm_86 — a likely scaling limiter to confirm
with `ncu` achieved-occupancy on the saturating (NOVA) regime. This is the quantitative
hook for roadmap step 3.1's "occupancy/register cost" and the `Matrix4` penalty question;
the heaviest kernels are AMReX MLTensorOp (elastic operator) — but those are not launched
in the elastic-disabled GPU configuration, so the phase-field kernels' register profile
should be re-measured in isolation on NOVA.

## Implication for D1 (elastic disposition)
D1 holds in 3D: the device elastic solve still fails (now confirmed in 3D, CUDA-719). The
GPU-resident workload for Phase 3 scaling is **phase-field + thermal with elastic
disabled**. A hybrid CPU-elastic / GPU-phase-field architecture remains the open path if
elastic is required in the scaled regime; it is out of scope for "get to testable state".

## Artifacts
- `bin/alamo_gpu-3d-profile-cuda86-g++` — verified 3D GPU binary.
- `input_3d_smoke` — minimal proven 3D GPU input (Expression ICs, 6-face BCs, elastic disabled).
- Build log: `/tmp/build_3d_cuda.log`; smoke log: `/tmp/smoke_3d.log`.

## Next (NOVA)
Production inputs (task 001), NOVA 3D build + multi-GPU launch (task 002), memory budget
+ crossover harness (task 003) → run the size×GPU sweep on A100/H200 → R3 crossover report.
