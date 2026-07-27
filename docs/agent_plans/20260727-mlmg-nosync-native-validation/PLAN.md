# TASK: mlmg-nosync-native-validation
# Folder: docs/agent_plans/20260727-mlmg-nosync-native-validation/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 3 (MLMG synchronization semantics and solver configuration)  |
| Model        | opus                                                         |
| Verification | partial-oracle                                               |
| Est. scope   | 2 solver headers, 2-4 benchmark scripts, task evidence; ~240 lines |
| Parallel-safe| no (local GPU and NOVA jobs are exclusive test resources)    |

## Operating rules

1. Read only the files listed in Context budget; `docs/archive/` is forbidden.
2. Preserve unrelated dirty work, especially `src/Integrator/Flame.{cpp,H}`.
3. Do not commit unless the user separately requests it.
4. A performance claim must exclude or amortize fixed startup cost by using
   multiple steps and reporting startup separately where the log permits it.
5. The opt-in path itself must pass sanitizer and physics validation.

## Context budget

Read first: `benchmark/status.sh` output, this `PLAN.md`
Read:
- `src/Solver/Nonlocal/Linear.H`
- `src/Solver/Nonlocal/MLMGSyncStateGuard.H`
- `ext/AMReX-Codes/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG.H:130-180,280-315,393-610`
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_GpuControl.H:140-195`
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_MFIter.{H,cpp}`
- `docs/agent_plans/20260727-amrex-gpu-guide-audit/{PLAN.md,results/RESULT.md,scripts/interleaved_ab.sh,scripts/knob_sweep.sh}`
- `benchmark/{NOVA_SLURM_RUNBOOK.md,local_a100_gate.sh,local_cuda_env.sh}`
- `docs/llm/PLAN.md`
- `benchmark/validate/{README.md,cases.manifest.yaml,physics_budget.yaml,run_validation_local.py,run_validation_local.sh,run_validation_nova.py,run_validation_nova.slurm,compare_validation.py}`
- `tests/ElasticSoftVoid/{input,test}`
- build/config files directly named by the NOVA runbook
Reference only if a validation script names it: its direct helper scripts and
the selected reference manifest under `benchmark/validate/references/`
Forbidden: `docs/archive/*`, unrelated task folders, user Flame changes

## Objective

Replace the successful custom/global-stream MLMG no-sync probe with the
solver-scoped AMReX 26.06 mechanism, preserving exception safety. Validate its
speed and correctness locally and on A100, including many-box and multi-rank
coverage. Make GPU timing harnesses use multiple steps and distinguish fixed
startup from steady simulation cost.

## Oracle

Commands:
- task-local interleaved A/B harness, at least 10 steps, 2D and 3D
- task-local many-box trace/full-output comparison
- recoverable MLMG-failure state-restoration test
- enabled-path compute-sanitizer and strict physics-budget comparison
- NOVA A100 performance and at least 2-rank/2-GPU smoke validation
- `benchmark/status.sh`

Covers: local and A100 wall behavior, solver traces and output fields, device
memory safety, exception restoration, one- and two-rank execution.
Does not cover: all production geometries or full production time horizons.

## Steps

### Step 1 - Implement native scoped control
VERIFY: confirm AMReX 26.06 provides `MLMG::setNoGpuSync` and inspect its
single-stream/no-sync restoration behavior.
DO: remove the custom global-stream guard; add an Alamo solver option routed
through `Linear::Parse` and `PrepareMLMG`; preserve state on exceptions.
CHECK: build 2D/3D GPU and CPU targets; default-off behavior unchanged.

### Step 2 - Make performance timing startup-aware
VERIFY: current audit harness uses two steps and external wall only.
DO: create/update harnesses to use at least 10 steps, record total wall and
per-step/steady-state timing, and pin the local arena without changing decks.
CHECK: shell/Python syntax checks and one smoke run per dimension.

### Step 3 - Local validation
VERIFY: local GPU is idle and binaries match the edited source.
DO: interleaved 2D/3D A/B; many-box A/B; force a recoverable failure and verify
the next solve runs with restored GPU sync state.
CHECK: identical iteration/residual traces where expected, field comparison,
sanitizer clean, and performance summary excluding startup distortion.

### Step 4 - A100 and multi-rank validation
VERIFY: follow `NOVA_SLURM_RUNBOOK.md`; record exact commit/worktree patch and
AMReX checkout used remotely.
DO: run opt-in strict physics gate, compute-sanitizer, multi-step performance,
and 2-rank/2-GPU smoke comparison.
CHECK: jobs exit zero, physics budget passes, sanitizer has zero findings,
multi-rank outputs agree within the declared oracle.

### Step 5 - Review and closeout
DO: line-by-line diff review; record all raw commands, job IDs, timings,
correctness evidence, limitations, and retain/revert verdict in `RESULT.md`.
CHECK: `benchmark/status.sh`; no unrelated files changed.

## Checkpoints

- [x] Plan restated and human confirmation received: user explicitly requested
      all six validation actions plus multi-step startup-aware timing.
- [x] Before source edit: native API and exception behavior documented.
- [x] After source edit: line-by-line diff and default-off build checks.
- [x] Before remote submission: local correctness and harness smoke pass.
- [x] Before retention verdict: A100, sanitizer, physics, and MPI evidence pass.

## Closeout

- [x] Oracle passes; `benchmark/status.sh` green
- [x] `results/RESULT.md`
- [x] `results/DONE`
- [x] Session log line appended to `docs/llm/SESSION_LOG.tsv`
