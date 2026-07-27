# TASK: phase0-baseline
# Folder: docs/agent_plans/20260727-phase0-baseline/

Phase 0 of the `chamber-gpu-mem` memory-strategy campaign
(`docs/agent_plans/20260727-gpu-memory-strategy/PLAN.md` §6). Runs after Phase
0.5 (`docs/agent_plans/20260727-phase05-two-rank-probe/`, PASS 2026-07-27).

---

## Header

| Field        | Value                                                                   |
|--------------|-------------------------------------------------------------------------|
| Risk tier    | 1 — measurement harness + docs. No `src/` change is authorized here.    |
| Model        | opus (campaign orchestration); capture legs are mechanical               |
| Verification | partial-oracle — `benchmark/status.sh` + figure-set completeness         |
| Est. scope   | this folder + NOVA batch script under `benchmark/`; ~400 lines           |
| Parallel-safe| no — local legs use `bin/alamo*` and race `status.sh` (campaign §2)      |

## Operating rules

1. Read ONLY the files listed in Context budget. `docs/archive/` is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: STOP,
   report discrepancy, wait.
3. One commit per step. Message: `<area>: <what> (20260727-phase0-baseline)`.
4. No scope expansion. New ideas go to the campaign `NOTES.md`, not into code.
   A defect found here gets its own tier-3 folder (campaign PLAN §1).
5. Missing knowledge → ask.
6. Tier 1: unattended, gates + spot-check diff. **Exception:** if Step L2 finds
   a step-loop region that needs a new `BL_PROFILE`, that edit is tier 2 and
   stops for a checkpoint first.

## Context budget

Read first: `benchmark/status.sh` output, this PLAN.md, campaign PLAN §6/§8/§12
Read: `benchmark/NOVA_SLURM_RUNBOOK.md`, `benchmark/nova_flame_gpu.slurm`,
`benchmark/build_alamo_nova{,_3d}.sh`, `benchmark/select_nova_resources.sh`,
`benchmark/phase3_memory_budget.py`, `benchmark/baseline_suite.py`,
`benchmark/ci_golden_compare.sh`, `configure` (profile/perf flags)
Reference only if a step names it: `ext/AMReX-Codes/amrex/Src/Base/AMReX_TinyProfiler.cpp`,
`src/Integrator/Flame.cpp`, `src/Integrator/Integrator.cpp`
Forbidden: `docs/archive/*`, unrelated task folders, the ~/Desktop sweep campaign

## Objective

Nothing in Phases 1-3 should be scoped by assumption. Phase 0 establishes what is
actually true of the current port: where time goes, what crosses the bus, what
the arena does, and how far the top kernels sit from the HBM roof. It produces
the gap table against T1-T7, sets the T7 threshold from a measured roofline
rather than a guess, and answers the campaign's own decision gate — whether
chamber throughput is worth optimizing at all versus `multicomponent/FMA`.

After this task: a populated gap table, a numeric T7 band, a reproducible
figure-set capture (one command to submit, one to collect), a footprint budget
against 8/40/64 GB, and a recorded go/no-go on Phase 1.

## Oracle

Command(s):
- `benchmark/status.sh` — three gates green, run **solo**
- Figure-set completeness: F1-F10 present for both 2D and 3D in
  `results/figures/`, plus one metrics-ledger row per `(phase, sha, case, ranks,
  device)`

Covers: that the baseline capture is complete and reproducible; that the gates
did not regress while the harness was built.

Does NOT cover: whether the *interpretation* of the baseline is right — the gap
table and the T7 threshold are judgment, reviewed by a human at the phase
boundary. No correctness property is proved by this task; it is measurement only.

## Findings already established (do not re-derive)

- **NVTX needs no source work.** Vendored AMReX 26.06
  `Src/Base/AMReX_TinyProfiler.cpp:20-27,134,211` pushes/pops an `nvtxRange`
  around every `BL_PROFILE` region under `AMREX_USE_CUDA`. No NVTX-specific
  build flag exists or is needed — `TINY_PROFILE=TRUE` (via `configure
  --profile`) is sufficient, and `bin/alamo_gpu-2d-profile-cuda86-g++` already
  has it. Verified locally 2026-07-27: an `nsys -t cuda,nvtx` capture yields
  named ranges for `Integrator::Evolve`, `Flame::TimeStepBegin`, `MLMG::solve()`,
  `MLMG::mgVcycle()`, `Operator::Fsmooth()`, `Operator::Elastic::Fapply()`,
  `FillBoundary`, and more. Campaign PLAN §6/§12 lists this as the highest-value
  QOL investment; it is already paid for.
- **Consequence:** the profile build is the *diagnostic* binary, not the timing
  binary. TinyProfiler adds a push/pop per region (16,304 `Fapply` ranges in a
  2-step 2D run) and `MODE=bench` in `nova_flame_gpu.slurm` additionally sets
  `tiny_profiler.device_synchronize_around_region=1`. **F10 wall-time must come
  from the non-profile binary**, or the headline number measures the profiler.
- Local Open MPI has no CUDA support; GPU-aware MPI is inactive on kermit
  (Phase 0.5 §2). The NOVA compute-node re-check is inherited by this task.
- `nova_flame_gpu.slurm` already carries `MODE=fast` with
  `amrex.the_arena_is_managed=0` — the T1 flip needs no new switch.

## Steps

Split into **local legs (L)** that cost no queue time and **NOVA legs (N)** that
are batched into one submission per campaign §2.

### Step L1 — NVTX verification

**Status: DONE 2026-07-27.** See "Findings already established". Evidence:
`nsys stats --report nvtx_sum` on a 2-step 2D profile-build capture.

### Step L2 — Step-loop region coverage audit

VERIFY:
```bash
grep -c BL_PROFILE src/Integrator/Flame.cpp src/Integrator/Integrator.cpp \
                   src/Solver/Nonlocal/Newton.H src/Operator/Elastic.cpp
```
DO: from the L1 trace's `nvtx_sum`, list which step-loop stages are named and
which are unattributed. Report gaps; do **not** add regions yet.
CHECK: `results/RESULT.md` §L2 lists named stages and gaps. If a gap is
load-bearing for F1, STOP and checkpoint — adding `BL_PROFILE` is a tier-2 edit.

### Step L3 — Footprint budget (campaign §8)

VERIFY: `python3 benchmark/phase3_memory_budget.py --help` exits 0.
DO: compute `resident = cells × components × (1+ghost) × sizeof(Real) × copies`
summed over levels and solver working copies, for the campaign 2D and 3D decks.
Evaluate against **8 GB** (which decks stay runnable locally, i.e. how much
correctness work can happen off-queue), **40 GB** (A100), **64 GB** (MI210),
each divided by ranks-per-GPU. Prefer a measured arena high-water over the
script's first-order model where one is available; label estimates as estimates.
CHECK: `results/RESULT.md` §L3 carries the table and an explicit verdict on
whether FP32 storage is structural or optional.

### Step L4 — `rod_and_tube_step2` gpu_strict reference

VERIFY:
```bash
python3 benchmark/baseline_suite.py --list
bash benchmark/ci_golden_compare.sh   # solo; note which legs run
```
DO: campaign §6 names this stale reference as a prerequisite. Determine whether
the GPU golden is stale, and either **re-record** it (if the current value is
demonstrably correct) or **quarantine** the case from the gpu_strict leg with a
recorded reason. Do not silently widen a tolerance.
CHECK: gpu_strict leg is either trustworthy or explicitly marked not-in-force,
with the decision written down. `benchmark/status.sh` still green.

### Step L5 — Compose the NOVA batch

VERIFY: `ssh -o BatchMode=yes nova hostname` succeeds over the ControlMaster
tunnel; `sinfo` shows a100/h200 GRES.
DO: add `benchmark/phase0_capture.slurm` (+ a local composer/collector) that, in
**one** submission per dimension, captures the whole F1-F10 set:
nsys (`-t cuda,nvtx,osrt,mpi` + unified-memory page-fault and CUDA-memory-usage
counters), an `ncu` pass restricted to the top kernels by time, an arena
high-water trace, and a non-profile timing run for F10. Pin the toolkit version
and verify every flag spelling once against the NOVA modules.
CHECK: `bash -n` clean; a dry-run prints every `srun`/`nsys`/`ncu` line without
executing. Per campaign §2 nothing is submitted until the dry-run is read.

### Step N1 — Baseline capture, managed arena as-is

DO: submit the L5 batch for 2D and 3D, `MODE=bench`, arena in its current
managed state. Record job ids; do not block on the queue.
CHECK: F1-F10 present for both dimensions; ledger rows written.

### Step N2 — Device-arena flip inventory

DO: same decks, `MODE=fast` (`amrex.the_arena_is_managed=0`) plus
`amrex.abort_on_out_of_gpu_memory=1`. **Inventory the crash set; fix nothing.**
Classify each failure by campaign §7's triage table (timestep loop / diagnostic /
checkpoint / debug leftover).
CHECK: `results/RESULT.md` §N2 lists every failure with a file:line and a class.

### Step N3 — Nsight Compute detail

VERIFY: counter permission on NOVA (`ERR_NVGPUCTRPERM` check) **before**
budgeting a job around `ncu`.
DO: `ncu` on the top kernels by time only, a handful of invocations. Produce F5,
F6, F7.
CHECK: roofline and achieved-bandwidth numbers exist for the top kernels in both
2D and 3D.

### Step N4 — GPU-aware MPI confirmation

DO: run the Phase 0.5 §2 probe on a **compute** node inside the batch.
CHECK: verdict recorded. It decides whether campaign §9's "device-buffer
Allreduce" is literal or debt.

### Step N5 — Gap table and T7 threshold

DO: populate the campaign §6 gap table for T1-T7 from N1-N4. Set the T7 band
from the measured roofline. Write the metrics-ledger rows.
CHECK: no gap-table cell is empty; T7 has a number and a derivation.

### Step D — Decision gate (campaign §6)

DO: answer, in writing, whether chamber gates simulation throughput for the SRM
paper. If it does not, the marginal value of Phase 1-3 is low and
`multicomponent/FMA` is the better claim on the time. **Readiness is not a
reason.** This is a user decision informed by N1-N5, not an agent decision.
CHECK: go/no-go recorded in `results/RESULT.md` before any Phase 1 folder opens.

## Checkpoints

- [ ] After L2, if a `BL_PROFILE` addition is proposed (tier-2 edit)
- [ ] After L5 dry-run, before the first `sbatch` (campaign §2: compose, then submit)
- [ ] At Step D, before Phase 1 opens

## Closeout

- [ ] Oracle passes; `status.sh` all green (solo run)
- [ ] Gap table complete; T7 threshold set with derivation
- [ ] F1-F10 captured in both 2D and 3D; ledger rows written
- [ ] `results/RESULT.md`: findings, evidence, deviations from plan
- [ ] `touch results/DONE`
- [ ] Session log line appended → `docs/llm/SESSION_LOG.tsv`
