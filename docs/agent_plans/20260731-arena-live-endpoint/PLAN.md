# TASK: arena-live-endpoint
# Folder: docs/agent_plans/20260731-arena-live-endpoint/

---

## Header

| Field        | Value |
|--------------|-------|
| Risk tier    | 2 — adds opt-in instrumentation to the CUDA launcher |
| Model        | opus |
| Verification | full-oracle — focused profile run plus FULL gate |
| Est. scope   | 2 source/harness files, analyzer/tests, <40 lines |
| Parallel-safe| no — rebuilds shared CUDA profile binaries |

## Operating rules

1. Read ONLY the files listed in Context budget. `docs/archive/` is forbidden.
2. Every step's VERIFY must pass before acting. On failure: stop and report.
3. One commit per step. Message:
   `benchmark: sample live arena state (20260731-arena-live-endpoint)`.
4. No scope expansion. New ideas go to this folder's `NOTES.md`.
5. The new output is opt-in and may not alter evolution or teardown order.
6. Tier 2: stop at the plan checkpoint and again before the source commit.

## Context budget

Read first: `benchmark/status.sh` output and this `PLAN.md`
Read: `src/alamo_gpu.cc`, `benchmark/phase0_capture.slurm`,
`benchmark/phase0_analyze.py`,
`ext/AMReX-Codes/amrex/Src/Base/AMReX_TinyProfiler.{H,cpp}`
Reference only if Step 2 names it: `IO/ParmParse.H`,
`benchmark/build_alamo_nova.sh`
Forbidden: `docs/archive/*`, unrelated task folders

## Objective

T5b requires the live arena allocation count at a timestep endpoint.  The
current profile prints only during finalization, after the integrator and its
persistent fields have been destroyed; because requests are balanced, AMReX
omits `Nfree` and `CurrentMem`, making the apparent flat result a false pass.
Add an opt-in pre-teardown sample in the CUDA launcher and consume it in the
Phase-0 arena leg.

## Oracle

Commands:

- one-step and five-solve profile runs with the opt-in flag;
- parser unit test proving `Nalloc`, `Nfree`, and `CurrentMem` are retained;
- `FULL=1 bash benchmark/status.sh`.

Covers: the live request count/bytes while `Integrator::Flame` and persistent
fields still exist, endpoint comparison, and no default-output behavior change.

Does NOT cover: live counts continuously between endpoints or non-Flame
launchers.

## Steps

### Step 1 — Add an opt-in pre-teardown sample

VERIFY: without the new flag, the finalizer table omits `Nfree/CurrentMem` after
balanced teardown; vendored `TinyProfiler::PrintMemoryUsage` includes those
columns when any live request remains.

DO: add a false-by-default CUDA-launcher parameter.  After `Evolve()` returns
and before deleting the integrator, call
`amrex::TinyProfiler::PrintMemoryUsage` only when enabled.  Label the output so
the endpoint table cannot be confused with finalization.

CHECK: default smoke output is unchanged; an opt-in one-step profile contains
the pre-teardown label and `Nfree/CurrentMem`.

### Step 2 — Wire capture and analysis

VERIFY: Step 1 output distinguishes the live table from the final table.

DO: enable the flag only in the Phase-0 arena leg; update the analyzer to parse
the labeled table, report live allocation count/bytes at step 1/full, and
compute the T5b delta.

CHECK: parser unit test passes; one-step/five-solve hardware output produces
non-missing T5b values.

### Step 3 — Gate and record

VERIFY: focused hardware evidence is green.

DO: run `FULL=1 bash benchmark/status.sh`; record diff/evidence and present the
pre-commit checkpoint.

CHECK: FULL is green and no source commit occurs before human confirmation.

## Checkpoints

- [ ] Plan checkpoint: human confirms opt-in launcher instrumentation before
      editing `src/alamo_gpu.cc`.
- [ ] Pre-commit checkpoint: diff, focused outputs, parser test, and FULL output
      shown for human confirmation.

## Closeout

- [ ] Default launcher behavior/output unchanged
- [ ] Pre-teardown Nalloc/Nfree/CurrentMem captured at both endpoints
- [ ] Analyzer reports T5b rather than inferring it from high-water
- [ ] FULL gate passes
- [ ] `results/RESULT.md` records evidence and limitations
- [ ] `touch results/DONE`
- [ ] Session log line appended to `docs/llm/SESSION_LOG.tsv`
