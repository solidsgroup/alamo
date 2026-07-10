# TASK: fsmooth-launch-fusion
# Folder: docs/agent_plans/20260709-fsmooth-launch-fusion/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 3 (solver: MLMG smoother)                                    |
| Model        | sonnet (user-directed; orchestrator reviews + verifier gate) |
| Verification | full-oracle (bit-exact: golden compare + lint + sanitizer)   |
| Est. scope   | 1 file (src/Operator/Operator.cpp), ~15 lines                |
| Parallel-safe| yes vs 20260709-fapply-kernel-surgery (disjoint file), stacked in same worktree |

## Operating rules

Same as 20260709-fapply-kernel-surgery/PLAN.md. Work in worktree
/home/jackplum/Projects/alamo-fapply-322b (branch fapply-322b, includes
commits 9470889b1/760d7f1cd). Main tree is off-limits for code (concurrent
session). No commit — orchestrator commits.

## Context budget

Read: src/Operator/Operator.cpp:340-420 (Fsmooth), this PLAN.md
Reference only if needed: benchmark/build_alamo_local_gpu.sh header
Forbidden: docs/archive/*, src/Operator/Elastic.*, src/Solver/*

## Objective

AMReX GPU guidance: on GPU the component loop should be fused into the
kernel (ParallelFor 4D/ncomp form) instead of launching one kernel per
component. `Operator<Grid::Node>::Fsmooth` (Operator.cpp:389-410) launches
the Jacobi update once per component per box per sweep (2 sweeps/call).
Fuse into a single `amrex::ParallelFor(bx, ncomp, ...)` launch. Bit-exact by
construction: the update at each (i,j,k,n) reads only its own component of
x/b/Rx/diag and writes only x(i,j,k,n) — no cross-component or cross-cell
dependence, so launch partitioning cannot change results.

## Oracle

Same gate set as 20260709-fapply-kernel-surgery (run in the worktree):
lint, GOLDEN_MODE=cpu ci_golden_compare.sh, supplemental full-solve 2D
memcheck (`compute-sanitizer --tool memcheck ./bin/alamo_gpu-2d-profile-cuda86-g++
input max_step=55 amr.plot_dt=1e9 plot_file=<scratch>`; must finish all steps,
0 errors). Timing: same deck/protocol as prior task (input, max_step=251,
3 runs, TinyProfiler; compare Fsmooth row and MLMG/total).

## Steps

### Step 1 - Fuse
DO: Replace the `for (int n...)` + 3D ParallelFor at Operator.cpp:389-410
with one `amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j,
int k, int n) {...})`; hoist the `auto m_omega = this->m_omega;` copy above
the MFIter loop body (it is loop-invariant). Do not alter the update
expression, the domain/bx branch structure, or anything else.
CHECK: 2D + 3D GPU builds compile.

### Step 2 - Gates + timing
DO: run oracle; save logs + wall tables to this folder's results/.
CHECK: all PASS; report.

## Closeout (orchestrator)

Line-by-line diff review, verifier if warranted, commit, RESULT.md, DONE.
