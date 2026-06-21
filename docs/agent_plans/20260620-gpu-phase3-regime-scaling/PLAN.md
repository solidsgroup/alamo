# Plan: GPU Phase 3 — Regime Scaling (bring to testable state)

## User goal
Execute Phase 3 of `GPU-OPT-ROADMAP.txt` ("Regime scaling: 2D→3D, single→multi-GPU")
to a **locally testable / NOVA-ready state**. We do NOT run the real scaling study
here — the A1000's 8 GB cannot reach the saturating regime (roadmap 3.2). The
deliverable is: a working 3D GPU binary, 3D wide-shallow inputs, NOVA 3D build +
multi-GPU launch scripts, a memory-budget calculator, and a scaling/crossover
harness — everything staged so that on NOVA (A100/H200) we just launch.

## Current architecture (from recon, 2026-06-20)
- Build: `./configure --comp=g++ --dim 3 --cuda 86 --profile && make -j` →
  `bin/alamo_gpu-3d-profile-cuda86-g++` (the `profile` infix comes from `--profile`).
  **VERIFIED 2026-06-20: this binary builds and runs a 3D GPU step successfully.**
  CUDA configure auto-restricts compilation to the
  GPU-safe source subset (`src/alamo_gpu.cc` + deps), which includes `Integrator::Flame`.
  Local CUDA env: `benchmark/local_cuda_env.sh` (CUDA 12.6.3 at `.local/cuda`, sm_86).
  Helper: `benchmark/build_alamo_local_gpu.sh` honors `DIM=3 ARCH=86 PROFILE=1 SMOKE=0`.
- 3D code readiness: linear algebra is fully 3D-ready and dimension-guarded —
  `Set::Matrix4<3>` (`src/Set/Matrix4_*.H`, `#if AMREX_SPACEDIM==3` paths),
  `Matrix4_Isotropic.H:148-153` computes the full 3×3 stress, quaternion ops are
  dimension-agnostic, `Flame.cpp` uses `AMREX_D_TERM`/`AMREX_SPACEDIM` loops throughout.
- **Blocker 1 for a 3D RUN:** `src/IC/BMP.H` is 2D-only (sets z=0, ignores k).
  The current `input` initializes `pf.eta` and `phi` from `.bmp` files, and its BC
  blocks (`pf.eta.bc`, `thermal.temp.bc`, `elastic.bc`) define only xlo/xhi/ylo/yhi.
  A 3D run therefore needs **analytic ICs** (e.g. `IC::Expression`) and **zlo/zhi BCs**.
- **Blocker 2 (the elastic trap), VERIFIED:** `elastic.on = 0` does NOT skip the elastic
  solve — `Mechanics::TimeStepBegin` runs the Static solve regardless and hits **CUDA-719**
  on device (the documented device-elastic crash). The correct switch is
  **`elastic.type = disable`** (`Mechanics::Type::Disable`, `src/Integrator/Base/Mechanics.H:28,58`),
  which early-returns from every elastic path. With `disable` the elastic `model_*` /
  `elastic.bc.*` blocks **must be OMITTED** — ALAMO's strict ParmParse aborts on the now-unused
  entries. The proven 3D template `input_3d_smoke` encodes exactly this.

## Desired architecture (end state of this plan)
1. `bin/alamo_gpu-3d-cuda86-g++` compiles locally and runs a 1-step smoke test. (3.1)
2. A canonical 3D-safe Flame input `input_3d_flame` (Expression ICs, full 3D BCs,
   wide-shallow grid + shallow AMR + large blocking_factor) plus grid-size variants
   for the crossover sweep. (3.3)
3. NOVA 3D build (`--dim 3`, sm_80 A100 + sm_90 H200) and single- + multi-GPU launch
   scripts. (3.4)
4. A memory-budget calculator that, from bytes/node, prints the largest grid that fits
   per device (A1000 8 GB, A100 40/80 GB, H100 80 GB, H200 141 GB). (3.2)
5. A scaling/crossover harness (problem-size × GPU-count driver) + R3 report skeleton. (3.5)

## Invariants and constraints
- **Branch:** all work is on `chamber-gpu` (NOT merged to master). Workers run in the
  shared working tree and create **only new, disjoint files** — no worktree isolation,
  because a worktree branched from the default branch would lack the GPU port.
- **No worker edits shared/existing files**; new files only, in the namespace assigned
  to that task. Two tasks never write the same path.
- **No worker builds or runs CUDA** — the lead owns the single 3D build + smoke + audit.
  Workers are pure artifact authors (inputs / scripts / docs) and must be runnable on a
  machine without a GPU.
- Canonical 3D input name is **`input_3d_flame`** (contract; workers reference it by
  name without needing the file to exist yet).
- Bytes/node figure for budgeting (from recon): elastic model fab `Matrix4<3>` ≈ 360 B/node;
  add eta + phi + temp (≈3 doubles = 24 B) + displacement vector (3 doubles = 24 B) +
  RHS/residual scratch. Budget conservatively ~450–512 B/node incl. ghost cells; the
  calculator must take bytes/node as a parameter, not hardcode.
- Wide-shallow box strategy (roadmap 3.3): large base grid, shallow AMR (max_level ≤ 1),
  `blocking_factor = 32`, big `max_grid_size`. Do NOT carry deep-AMR/tiny-box into 3D.
- Don't break the existing 2D workflow (don't touch `input`, `build_alamo_nova.sh`, or
  the existing `.slurm` files — create 3D-suffixed siblings).

## Files involved
- Lead (foundational): `bin/alamo_gpu-3d-cuda86-g++` (build artifact),
  `input_3d_smoke` (minimal smoke input), `benchmark/PHASE3_3D_READINESS.md` (R: 3.1).
- Task 001 (inputs): repo-root `input_3d_flame`, `input_3d_flame_*` grid variants.
- Task 002 (NOVA): `benchmark/build_alamo_nova_3d.sh`, `benchmark/nova_flame_gpu_3d.slurm`,
  `benchmark/nova_flame_gpu_3d_multi.slurm`, `benchmark/nova_flame_cpu_3d.slurm`.
- Task 003 (budget+harness): `benchmark/phase3_memory_budget.py`,
  `benchmark/phase3_scaling_sweep.sh`, `benchmark/PHASE3_R3_crossover.md` (skeleton).

## Build and test commands
- 3D build: `DIM=3 PROFILE=1 SMOKE=0 ARCH=86 ./benchmark/build_alamo_local_gpu.sh`
- 3D smoke: `source benchmark/local_cuda_env.sh && mpiexec -np 1 ./bin/alamo_gpu-3d-cuda86-g++ input_3d_smoke max_step=1 stop_time=1e-12_s amr.plot_int=-1 amr.thermo.plot_int=-1`
- Worker self-check (no GPU): bash `-n` syntax check on scripts; `python3 benchmark/phase3_memory_budget.py --help`; input files validated by structure review against the lead's `input_3d_smoke` template.

## Task graph
- **Lead-FOUND** (serial, critical path): 3D build + `input_3d_smoke` + smoke run + 3.1 audit.
  Produces the working binary, a known-good 3D input template, and the bytes/node ground truth.
- **001 inputs** (parallel-safe) — depends only on the PLAN contract + lead's smoke template.
- **002 NOVA scripts** (parallel-safe) — depends only on PLAN build command + input name.
- **003 budget+harness** (parallel-safe) — depends only on PLAN bytes/node + input name.
- **INTEGRATION** (lead): review RESULTs + diffs, consistency, recommend merge order.

## Parallelization strategy
Run Lead-FOUND first (it de-risks the build and yields the template + bytes/node).
Then dispatch 001, 002, 003 concurrently — file-disjoint, no build, no GPU.

## Integration strategy
Lead reads each `results/NNN-RESULT.md` + diff, checks each artifact against the 3D
template and the PLAN contract (input name, bytes/node, box strategy), writes
`INTEGRATION.md`. Do not auto-merge; recommend order. All artifacts are additive
(new files) so merge risk is low.

## Risks
- 3D CUDA compile failure (low — recon says LA is 3D-ready; risk is the GPU source
  subset or an untested `#if` path). Mitigation: lead owns build, fixes in `src/` if needed.
- BMP-IC trap: any 3D input that keeps `bmp` ICs will fail at runtime. Mitigation:
  Expression ICs mandated; lead's smoke input is the proven template.
- Workers diverging on input schema. Mitigation: lead publishes `input_3d_smoke` as the
  canonical template before dispatch.

## Rollback plan
All deliverables are new files on `chamber-gpu`; rollback = `git clean`/delete the new
files and the `obj-3d-cuda86-g++` build dir. No existing 2D file is modified.
