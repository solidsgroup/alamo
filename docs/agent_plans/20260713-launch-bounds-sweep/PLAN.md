# TASK: launch-bounds-sweep
# Folder: docs/agent_plans/20260713-launch-bounds-sweep/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 2 (src/Operator launch plumbing; no numerics change)         |
| Model        | sonnet (opus review if wrapper design gets hairy)            |
| Verification | partial-oracle (CPU golden compare + gates scripted; A100 perf verdict is judgment) |
| Est. scope   | src/Operator/Elastic.cpp, src/Operator/Operator.cpp, 1 new small header; ~150 lines |
| Parallel-safe| no (touches same files as any elastic work)                  |

## Operating rules

1. Read ONLY the files listed in Context budget. docs/archive/ is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: STOP,
   report discrepancy, wait.
3. One commit per step unless stated. Message: `<area>: <what> (<task-folder>)`.
4. No scope expansion. New ideas go to NOTES.md here, not into code.
5. Tier 2: checkpoint at plan restatement + before each commit.

## Context budget

Read first: benchmark/status.sh output, this PLAN.md
Read: src/Operator/Elastic.cpp (Fapply ~:502, Diagonal ~:836, ALAMO_ELASTIC_OP_FOR
      macro ~:158), src/Operator/Operator.cpp (Fsmooth ~:389),
      ext/AMReX-Codes/amrex/Src/Base/AMReX_GpuLaunchGlobal.H:14-26,
      ext/AMReX-Codes/amrex/Src/Base/AMReX_GpuLaunch.H:36-44,
      ext/AMReX-Codes/amrex/Src/Base/AMReX_MFParallelForG.H:58-107
Reference only if step names it: docs/agent_plans/20260707-a100-ab-c1/results/RESULT.md
Forbidden: docs/archive/*, global AMReX header edits

## Objective

PLAN.md task 3.2 / v3 3.F. C1 edits cut registers 255->244 but occupancy is
flat (~12%, 1 block/SM). Sweep __launch_bounds__(256, {1,2,3,4}) on the
elastic kernels (Fapply, Diagonal, Fsmooth) on A100 to find whether forcing
2+ blocks/SM (register cap 128/85/64) wins wall time despite spills, or
whether spill traffic dominates. Judged per point: A100 wall/step + ncu
achieved occupancy + registers + spill bytes, plus the budget gate.

## Implementation constraints (from recon, 2026-07-13)

- AMReX provides `launch_global<MT, min_blocks>` (AMReX_GpuLaunchGlobal.H:23-25)
  with `__launch_bounds__(MT, min_blocks)`. `amrex::ParallelFor<MT>` selects MT
  only — there is NO ParallelFor path to min_blocks. Do NOT patch AMReX headers.
- Mechanism: small alamo-local launch helper (new header, e.g.
  src/Operator/ElasticLaunch.H) that replicates ParallelFor's 3D (and 4D for
  Fsmooth) box-to-thread mapping and launches via
  `AMREX_LAUNCH_KERNEL(MT, ...)` -> `launch_global<MT, MIN_BLOCKS>`, with
  MIN_BLOCKS a compile-time knob (e.g. -DALAMO_ELASTIC_MIN_BLOCKS=n, default 0
  = today's behavior via plain ParallelFor).
- Scope: exactly the launch sites Elastic.cpp:502 (Fapply), :836 (Diagonal),
  Operator.cpp:389 (Fsmooth). All other kernels untouched.
- CPU build must be completely unaffected (helper collapses to the existing
  macro/ParallelFor path when !ALAMO_GPU).
- min_blocks=1 arm must reproduce today's kernel (sanity anchor): identical
  registers to baseline, wall within noise.
- "Statically spill-free to at least 80 registers" gate: ptxas -v output per
  arm; record registers + spill stores/loads bytes per kernel per arm.

## Oracle

Command(s): benchmark/ci_golden_compare.sh (CPU, bit-exact — launch config
  cannot change results; any diff = mapping bug in the helper);
  benchmark/lint_device_patterns.sh; compute-sanitizer full-solve on A1000.
Covers: correctness of the index mapping, no device bug classes.
Does NOT cover: perf verdict (A100 wall + ncu, judgment); occupancy claims.

## Steps

### Step 1 - Launch helper + wire 3 sites, min_blocks default off
VERIFY: status.sh gates green; fapply-322b merged into chamber-gpu (this task
  builds on the merged surgery kernels).
DO: add helper header; switch the 3 sites; default build (no flag) must produce
  byte-identical PTX or at minimum identical registers/spills for all 3 kernels
  (check ptxas -v diff) and pass CPU golden compare bit-exact.
CHECK: ptxas -v diff clean vs pre-change build; golden compare green.

### Step 2 - Local ptxas sweep (no GPU runs)
DO: build min_blocks = 1,2,3,4 for sm_80; record ptxas registers +
  spill bytes per kernel per arm into results/ptxas_table.md. Any arm spilling
  while registers > 80... record; arms that spill are still swept on A100 (the
  question IS whether occupancy beats spill), but flag them.
CHECK: table complete, 4 arms x 3 kernels.

### Step 3 - A100 sweep on NOVA
DO: reuse the alamo-322b-ab NOVA setup (bundle + arm-dir pattern,
  docs/agent_plans/20260713-fapply-322b-a100/NOTES.md): one checkout per arm
  (build flag differs), deck input_3d_centre_bore_256_a2, wall job + ncu job
  (--nvtx-include "Operator::Elastic::Fapply()/") per arm. min_blocks=1 arm is
  the baseline anchor.
CHECK: all jobs COMPLETED; TinyProfiler + ncu tables per arm.

### Step 4 - Verdict + budget gate
DO: pick winner by Fapply exclusive wall; confirm stress parity (fcompare vs
  min_blocks=1 arm); run budget gate on winner; results/RESULT.md.
CHECK: oracle green on winner; RESULT.md tables complete.

## Checkpoints (tier 2)

- [ ] After plan restatement, before edits
- [ ] Before each commit: diff summary + oracle output

## Closeout

- [ ] Oracle passes; status.sh all green
- [ ] results/RESULT.md: ptxas table, wall table, ncu table, verdict
- [ ] touch results/DONE
- [ ] Session log line appended
