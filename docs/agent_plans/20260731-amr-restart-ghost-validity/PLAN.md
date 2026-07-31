# TASK: amr-restart-ghost-validity
# Folder: docs/agent_plans/20260731-amr-restart-ghost-validity/

---

## Header

| Field        | Value |
|--------------|-------|
| Risk tier    | 2 — generic restart lifecycle in `src/Integrator/` |
| Model        | opus |
| Verification | full-oracle — CPU/GPU AMR restart plus FULL campaign gate |
| Est. scope   | 1 source file, focused regression evidence, <40 lines |
| Parallel-safe| no — gates build and run shared CPU/GPU binaries |

## Operating rules

1. Read ONLY the files listed in Context budget. `docs/archive/` is forbidden.
2. Every step's VERIFY must pass before acting. On failure: stop and report.
3. One commit per step. Message:
   `<area>: <what> (20260731-amr-restart-ghost-validity)`.
4. No scope expansion. New ideas go to this folder's `NOTES.md`.
5. Missing knowledge means stop and ask; do not weaken the restart oracle.
6. Tier 2: stop after the plan checkpoint and before the source commit.

## Context budget

Read first: `benchmark/status.sh` output and this `PLAN.md`
Read: `src/Integrator/Integrator.cpp`, `src/Integrator/Integrator.H`,
`src/Integrator/Flame.H`, `src/Integrator/Flame.cpp`,
`tests/GPU/C2_restart_roundtrip/{input,test.py}`,
`tests/GPU/testlib_gpu.py`
Reference only if Step 1 names it: AMReX FillPatch documentation and existing
`FillPatch`, `FillBoundary`, `AverageDown`, and restart call sites under
`src/Integrator/`
Forbidden: `docs/archive/*`, unrelated task folders

## Objective

A checkpoint containing level-1 cell and nodal state loads, but the first
post-restart phase-field advance observes non-finite data on level 1 on both CPU
and GPU.  Identify which restored field or halo is invalid, restore the generic
post-restart invariants without altering fresh-run evolution, and make the
two-rank AMR C2 roundtrip pass.

## Oracle

Commands:

- one-rank CPU restart from a two-level step-10 checkpoint;
- direct two-rank strict-GPU C2 continuous/restart test;
- `FULL=1 bash benchmark/status.sh`.

Covers: valid first advance after restoring level-0/level-1 cell and nodal
checkpoint state, CPU/GPU behavior, redistribution across ranks, and the full
campaign correctness gate.

Does NOT cover: restarting every integrator type or every possible refinement
layout.

## Steps

### Step 1 — Identify the missing post-restart invariant

VERIFY: reproduce the level-1 step-11 finite-value abort on CPU and strict GPU.

DO: trace generic restart construction through field registration, valid-region
readback, physical/inter-box halo fill, and first `Flame::Advance`.  Use
launch-blocking only to localize asynchronous GPU reporting.  Record the exact
field and invariant in `results/RESULT.md`.

CHECK: a minimal diagnostic distinguishes invalid valid data from invalid ghost
data and identifies the earliest lifecycle point where the invariant is
missing.

### Step 2 — Restore restart field validity

VERIFY: Step 1 has one root cause and the proposed change does not alter
fresh-run initialization.

DO: make the smallest generic restart-lifecycle change that establishes the
same field/halo state expected before a normal advance.  Do not special-case C2
or weaken the non-finite tripwire.

CHECK: one-rank CPU and two-rank strict-GPU continuous/restart runs pass; the
step-10 checkpoint contains `Level_1` for cell and nodal data.

### Step 3 — Correctness gate and closeout

VERIFY: focused restart checks are green.

DO: run `FULL=1 bash benchmark/status.sh`, record commands and evidence, and
commit only after the pre-commit checkpoint is confirmed.

CHECK: FULL gate passes with no sanitizer or golden-comparison failure.

## Checkpoints

- [ ] Plan checkpoint: approach restated; human confirms before source edit.
- [ ] Pre-commit checkpoint: focused CPU/GPU evidence, FULL output, and diff
      summary shown; human confirms before commit.

## Closeout

- [ ] Focused CPU and two-rank strict-GPU AMR restart checks pass
- [ ] `FULL=1 benchmark/status.sh` is green
- [ ] `results/RESULT.md` records root cause, change, evidence, and limitation
- [ ] `touch results/DONE`
- [ ] Session log line appended to `docs/llm/SESSION_LOG.tsv`
