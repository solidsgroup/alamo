# Plan: chamber-gpu — Lever 2.5 (field packing) + Phase 1 (elastic disposition), parallel

## User goal
Resume the GPU-optimization roadmap (`~/Desktop/GPU-OPT-ROADMAP.txt`) from the
2026-06-20 handoff (`~/Desktop/SESSION_HANDOFF_2026-06-20.md`). User's standing
instruction: "Spin up an agent to implement 2.5 now, spin up a second agent to
move to phase 1 elastic." Run both **in parallel**. 2.5 is wanted for diligence
even though it is low-ROL — a well-documented negative result is acceptable.

## Current architecture
- Branch `chamber-gpu` @ `fe3844f31` (checkpoint of validated Phase 0 + 2.3 work).
- Tracked tree is CLEAN. Untracked benchmark/analysis tooling + input files are
  present in the main tree and are used in place (not committed).
- GPU: single RTX A1000, **8 GB**, sm_86. Binaries (built from the checkpoint):
  - `bin/alamo-2d-clang++`              — CPU, run np8
  - `bin/alamo_gpu-2d-cuda86-g++`       — fast (--use_fast_math)
  - `bin/alamo_gpu-2d-nofast-cuda86-g++`— strict (--fmad=false, no fast-math)
- Roadmap status: G0 passed (gate 9/9, box sweep 5/5). Phase 2.2 + 2.3 done.
  D2 verdict: phase-field path is launch-bound but the residual is a **regime
  problem (Phase 3)**, not fixable by 2D launch-code refactors. 2.5 is the last
  launch-code lever (~25% of launches best case, high refactor risk).

## Desired architecture
- **Track A (lever 2.5):** co-evolving cell-centered ghost-bearing fields packed
  into multi-component MultiFabs so one FillBoundary/FillPatch covers them,
  cutting the 56%-of-launches FillBoundary population — IF it stays numerically
  transparent (golden compare preserved) and yields a real wall/step win.
- **Track B (Phase 1):** elastic-path disposition decided on data (D1 verdict:
  GPU | CPU | HYBRID), with the no-fast-math high-res convergence result, ncu
  occupancy/register findings, and a fair GPU-vs-multicore-CPU benchmark.

## Invariants and constraints
- **Numerical transparency (principle 5):** every optimization re-runs the golden
  compare. no-fast-math (strict) gate must stay **bit-for-bit** vs references.
- **Device safety:** no host-loop writes into device-arena fabs reachable from
  the canonical input. Keep device-callable model evaluators free of host-only
  `Util::*` helpers.
- **Single shared GPU (8 GB):** both tracks contend for one A1000. Be robust to
  transient `cudaErrorMemoryAllocation`/OOM — retry, reduce resolution/concurrency,
  and never assume exclusive GPU access.

## Files involved
- Track A: `src/Integrator/Flame.cpp`, `src/Integrator/Flame.H`, BC construction
  for eta/psi/temp; reference only: `src/Integrator/Integrator.cpp`
  (`RegisterNewFab` ~340, `physbc_array` FillPatch ~394/1218).
- Track B: NO source edits. Input copies prefixed `input_phase1_*`; outputs under
  `benchmark/phase1_elastic_2048/` and `analysis/results_phase1*/`.

## Build and test commands
```bash
source benchmark/local_cuda_env.sh                       # project-local CUDA on PATH
# Rebuild (Track A only):
./configure --comp=clang++ --dim 2 && make -j"$(nproc)"   # CPU
COMP=g++ CUDA_FP=fast   SMOKE=0 ./benchmark/build_alamo_local_gpu.sh   # gpu fast
COMP=g++ CUDA_FP=strict SMOKE=0 ./benchmark/build_alamo_local_gpu.sh   # gpu strict
# Gates:
CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py check  # expect 9/9 ok
python3 benchmark/phase2_box_sweep.py                                             # expect 5/5
# nsys launch/sync capture: NSYS=.local/nsight/opt/nvidia/nsight-systems/2026.1.3/target-linux-x64/nsys
# ncu occupancy:           NCU=.local/nsight-compute/usr/bin/ncu  (benchmark/g0_ncu_capture.sh)
```

## Task graph
- **A — Lever 2.5 field packing** — `tasks/A-lever25-field-packing.md` — `parallel-safe` (owns src + builds)
- **B — Phase 1 elastic disposition** — `tasks/B-phase1-elastic-disposition.md` — `parallel-safe` (no src, no builds)

## Parallelization strategy (main-tree, no worktrees)
Both agents run in the main tree concurrently. Isolation is by **ownership**, not
filesystem:
- **Track A owns:** all of `src/`, all `./configure`/`make`/`build_alamo_local_gpu.sh`
  runs, all writes to `bin/`, `tmp_build_dir/`, `.alamo`, `Makefile`,
  `benchmark/phase2_box_sweep/`, and any `benchmark/golden_*` it creates.
- **Track B owns:** a **private copy** of the nofast binary (copy it first, e.g.
  `cp bin/alamo_gpu-2d-nofast-cuda86-g++ bin/alamo_gpu_phase1-nofast`), plus
  `benchmark/phase1_elastic_2048/`, `analysis/results_phase1*/`, and input copies
  prefixed `input_phase1_*`. **Track B must NOT build, configure, or edit `src/`,
  and must NEVER modify the shared `input` file** (Track A's gate reads it).
- Neither agent commits. Each leaves its changes for lead review and writes its
  RESULT file.

## Integration strategy
Lead reads `results/A-RESULT.md` + `results/B-RESULT.md` and the Track A working-tree
diff, writes `INTEGRATION.md`, and recommends a merge/commit order. Track A's 2.5
change is committed only if golden compare stayed green AND wall/step improved ≥5%;
otherwise it is reverted and the negative result is documented. Track B produces
report R1 (elastic disposition) — committed as documentation regardless of verdict.

## Risks
- 2.5 refactor breaks numerics (golden compare) — mitigate: implement the safest
  same-BC pack first, re-gate after each change, revert on any golden break.
- 2048² elastic doesn't fit in 8 GB — mitigate: Track B runs the largest fitting
  resolution and documents the memory ceiling; the convergence question is
  answerable at 1024² too.
- GPU contention between tracks — mitigate: Track A is compile-bound for its first
  ~10–15 min (no GPU) while Track B runs its GPU jobs; both retry on OOM.

## Rollback plan
- Track A: `git checkout -- src/` (Track B never touches src, so the src diff is
  purely A's). Discard A's binaries by rebuilding from `fe3844f31`.
- Track B: delete `input_phase1_*`, the private binary copy, and new output dirs.
- The checkpoint `fe3844f31` is the known-good restore point for both.
