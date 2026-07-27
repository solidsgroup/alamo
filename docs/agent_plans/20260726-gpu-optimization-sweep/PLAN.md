# TASK: gpu-optimization-sweep
# Folder: docs/agent_plans/20260726-gpu-optimization-sweep/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 3 (solver: Operator::Fsmooth, Elastic::Fapply)               |
| Model        | opus                                                         |
| Verification | partial-oracle (golden compare + device lint + sanitizer)    |
| Est. scope   | src/Operator/Operator.cpp, src/Operator/Elastic.cpp; ~120 ln |
| Parallel-safe| no (shares Elastic/Operator with any other perf task)        |

## Context budget

Read first: `benchmark/status.sh` output, this PLAN.md
Read: `src/Operator/Operator.cpp:344-420`, `src/Operator/Elastic.cpp:174-380`,
`src/Set/Matrix4_Major.H`, `src/Numeric/Stencil.H`,
`docs/agent_plans/20260721-fapply-runtime-optimization/results/RESULT.md`
Reference only if step names it: `src/Solver/Nonlocal/Newton.H`,
`benchmark/ci_golden_compare.sh`, `benchmark/local_a100_gate.sh`
Forbidden: `docs/archive/*`, unrelated task folders

## Objective

Fapply source-level work (tasks 3.1/3.2b/3.3) has taken the elastic kernel from
299 s to 255 s exclusive and the retained smoothing configuration cut MLMG solve
another 22%. What has never been attempted is the *rest* of the GPU picture:
the 2D-conservative nsys trace shows 43,212 kernel launches in a run with only
438 ms of GPU-busy time and 1.50 s of wall, and 51.9 ms (11.8% of GPU time,
~30k launches) is spent in `amrex::Copy`/`Multiply`/`Subtract` helper kernels
that exist only to stage two temporary MultiFabs inside `Operator::Fsmooth`.
After this task the smoother must compute the same values with one kernel
instead of six per Jacobi half-sweep, and the Fapply candidates listed under
Step 3 must be measured (accepted or rejected on evidence, not merged blind).

## Oracle

Command(s):
- `benchmark/lint_device_patterns.sh` (exit 0)
- `GOLDEN_MODE=cpu benchmark/ci_golden_compare.sh` (exit 0)
- `GOLDEN_MODE=gpu benchmark/ci_golden_compare.sh` (exit 0, strict build)
- `TIERS=1 benchmark/local_a100_gate.sh` (exit 0)
Covers: bit/tolerance-strict agreement with recorded CPU and GPU-strict
baselines across the golden case set; device anti-pattern regressions; runtime
strict smoke on the local sm_86 GPU.
Does NOT cover: A100/sm_80 behaviour, long-horizon stability, multi-rank,
performance itself (measured separately, medians of >= 5 reps).

## Steps

### Step 1 - Baseline capture (no source change)
DO: build `bin/alamo_gpu-{2,3}d-profile-cuda86-g++` at HEAD; record medians of
5 timed reps for the 2D-conservative and a memory-resident 3D case, plus an
nsys kernel trace and an ncu Fapply capture (LaunchStats/Occupancy/
SpeedOfLight/MemoryWorkloadAnalysis/WarpStateStats).
CHECK: `artifacts/baseline/` holds the trace, the ncu report, and a timing
summary with per-arm MAD.

### Step 2 - Fsmooth elementwise fusion
DO: in `src/Operator/Operator.cpp:356-402` delete the `Dx` and `Rx` MultiFabs
and the four `Copy`/`Multiply`/`Copy`/`Subtract` calls; compute
`Rx = Ax - x*diag` inside the existing update `ParallelFor` with the same two
IEEE operations in the same order (no FMA contraction).
CHECK: full oracle; kernel count from a fresh nsys trace drops by ~4 per
Jacobi half-sweep; timing medians vs Step 1.

### Step 3 - Fapply candidates (measure, then accept or reject)
3a. bind `DDW(i,j,k)` by const reference instead of copying a 45-double
    Matrix4 into registers (`Elastic.cpp:242`).
3b. any candidate the Step-1 ncu capture actually justifies (tiling, launch
    shape, L1/L2 policy). No speculative rewrites.
CHECK: per candidate - registers/thread and achieved occupancy from ncu, wall
medians, and the full oracle. Anything under a 3% median gain is reverted.

### Step 4 - Closeout
CHECK: `benchmark/status.sh` green; `results/RESULT.md` records accepted and
rejected candidates with evidence paths.

## Checkpoints (tier 3)

- [ ] Baseline recorded before any edit
- [ ] Step 2 diff reviewed line-by-line against the arithmetic it replaces
- [ ] Oracle output pasted before each commit
- [ ] Rejected candidates reverted, not left behind a flag

## Adversarial review

Fresh session, after implementation:
"Review the Fsmooth fusion on chamber-gpu. Assume it contains a defect. Check
the ghost-cell region the fused kernel reads versus what Copy/Multiply/
Subtract wrote, FMA contraction changing the rounding of `Ax - x*diag`,
the `relax_ghost_rows` branches, and whether any gate was weakened."

## Closeout

- [ ] Oracle passes; status.sh green
- [ ] results/RESULT.md
- [ ] changelog/ entry (append-only)
- [ ] touch results/DONE
- [ ] SESSION_LOG.tsv line appended
