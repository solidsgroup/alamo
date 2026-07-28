# TASK: elastic-psi-stencil-oob
# Folder: docs/agent_plans/20260728-elastic-psi-stencil-oob/

---

## Header

| Field        | Value                                                                |
|--------------|----------------------------------------------------------------------|
| Risk tier    | 3 — `src/Operator/Elastic.cpp` device kernels + `src/Numeric/Stencil.H` |
| Model        | opus                                                                 |
| Verification | partial-oracle — `FULL=1 benchmark/status.sh` green + golden compare |
| Est. scope   | 1-3 files, <40 lines                                                 |
| Parallel-safe| no — rebuilds `bin/`, races every gate                               |

## Operating rules

1. Read ONLY the files in Context budget. `docs/archive/` is forbidden.
2. Every step's VERIFY must pass before acting. On failure: STOP, report, wait.
3. One commit per step. Message: `<area>: <what> (20260728-elastic-psi-stencil-oob)`.
4. No scope expansion. New ideas go to this folder's NOTES.md.
5. Missing knowledge → ask.
6. Tier 3: stop at every checkpoint. Do not continue on self-assessment.

## Context budget

Read first: `FULL=1 bash benchmark/status.sh` output, this PLAN.md
Read: `src/Operator/Elastic.cpp:190-270,430-530,650-690`,
`src/Numeric/Stencil.H:1570-1600`, `src/Solver/Nonlocal/Newton.H:1400-1420`,
`docs/agent_plans/20260727-gpu-memory-strategy/NOTES.md` (N11 only)
Reference only if a step names it: `src/Integrator/Base/Mechanics.H:200-240`,
`src/Operator/Operator.cpp:500-520`, `benchmark/local_a100_gate.sh`
Forbidden: `docs/archive/*`, unrelated task folders, the ~/Desktop sweep campaign

## Objective

`Operator::Elastic<1>::Diagonal()` reads out of bounds on the 3D deck during the
first elastic solve: 2,497 × `Invalid __global__ read of size 8 bytes`, then
`cudaErrorLaunchFailure (error 719)`. The faulting address is 18,320 bytes
*before* the nearest allocation, i.e. below the fab. This is a live GPU
correctness defect that was invisible because `benchmark/status.sh` never ran
the sanitizer tier (fixed 2026-07-28, `d7fe4a068`).

After this task: `FULL=1 bash benchmark/status.sh` is green, the fix is
justified against the psi/stencil contract rather than chosen to silence the
sanitizer, and the other `CellToNodeAverage` call sites are adjudicated
explicitly rather than left as latent copies of the same bug.

## Oracle

Command(s):
- `FULL=1 bash benchmark/status.sh` — three legs green, run **solo**
- `GOLDEN_MODE=gpu bash benchmark/ci_golden_compare.sh` — no golden regression

Covers: that the OOB read is gone under compute-sanitizer on the 3D multi-box
layout; that the fix did not move any golden value.

Does NOT cover: whether the *numerics* at psi boundaries are now correct, only
that the reads are in bounds. A fix that reads in-bounds garbage passes this
oracle. Boundary values must be argued from the psi contract and reviewed by a
human at the Step 3 checkpoint.

## Evidence already in hand (do not re-derive)

- Path: `Mechanics.H:225` → `MLMG::prepareForSolve` → `Operator.cpp:515` →
  `Elastic.cpp:441` (ParallelFor) → device frame `Elastic.cpp:461` →
  `Stencil.H:1592`.
- `Stencil.H:1587-1592`: `ilo/jlo/klo = (stencil[d] == StencilType::Lo ? 0 : 1)`,
  then reads `f(i-ilo, j-jlo, k-klo, m)`. With `DefaultType()` (all Central)
  this is `(i-1,j-1,k-1)` unconditionally.
- `Elastic.cpp:461` calls `CellToNodeAverage(psi, i, j, k, 0)` with **no**
  stencil argument, while the same lambda built the boundary-aware `sten` at
  `:444-445` and passed it to `Gradient_Diagonal` at `:448`.
- Systemic: of 19 `CellToNodeAverage` sites in `src/`, only
  `Newton.H:1411` passes a stencil. `PhaseFieldMicrostructure.cpp:270-273` has
  `//, sten);` commented out.
- Log: `benchmark/_a100_gate_20260728_125617/tier2_memcheck.log`.

## Steps

### Step 1 — Establish provenance

VERIFY:
```bash
git status --porcelain src/          # expect only Flame.{cpp,H} dirty
```
DO: stash the dirty `src/` changes, rebuild 3D CUDA, re-run `TIERS=2`
`local_a100_gate.sh`, restore the stash. Determines whether the defect is
pre-existing on the branch or introduced by uncommitted work.
CHECK: `results/RESULT.md` §1 records the clean-tree verdict with the log path.
If the clean tree is GREEN, STOP and checkpoint — the diagnosis changes.

### Step 2 — Characterize the read, do not fix yet

VERIFY: Step 1 recorded.
DO: determine what `psi` actually is at `Elastic.cpp:461` — its ghost cell
count, its BoxArray relative to `tilebox`, and whether the low-corner node of
the lowest box can legally index `(i-1,j-1,k-1)`. Establish whether the three
Elastic sites (`:255`, `:461`, `:669`) sit on the same box/ghost geometry or
different ones.
CHECK: `results/RESULT.md` §2 states, with file:line, why the read goes out of
bounds and whether the other two sites are exposed to the same condition.
**Checkpoint: report before proposing a fix.**

### Step 3 — Choose and justify the fix

VERIFY: Step 2 checkpoint cleared by the user.
DO: choose among — pass `sten`; grow psi's ghost region; clamp the read — and
write why the other two are wrong *for this call site*. Note that passing
`sten` changes the interpolation stencil at boundaries and therefore the psi
value, which is a numerics change, not just a bounds fix.
CHECK: `results/RESULT.md` §3 carries the decision and the rejected
alternatives. **Checkpoint: user confirms before any edit.**

### Step 4 — Implement

VERIFY: Step 3 checkpoint cleared.
DO: apply the fix. Apply or explicitly decline it at `:255` and `:669`, with a
recorded reason per site.
CHECK: `FULL=1 bash benchmark/status.sh` green, solo. `bash -n` n/a; device
lint clean.

### Step 5 — Adjudicate the remaining call sites

VERIFY: Step 4 green.
DO: decide whether `PhaseFieldMicrostructure.cpp:270-273`'s commented-out
`//, sten);` and the other 14 stencil-less sites are the same defect dormant.
Log verdicts; open a follow-on folder if any is live. Do not fix them here.
CHECK: `results/RESULT.md` §5 lists each site with a verdict.

## Checkpoints

- [ ] After Step 1, if the clean tree is green (diagnosis changes)
- [ ] After Step 2, before any fix is proposed
- [ ] After Step 3, before any edit
- [ ] Before commit: diff summary + `FULL=1` output

## Adversarial review (tier 3, mandatory)

Fresh session, no context from this folder:
  "Review commit <hash> on chamber-gpu-mem. Assume it contains a defect. Find
   it. Check: device-lambda captures, elixir lifetimes, Eigen expression
   chaining, MLMG hierarchy freshness, tolerance changes, boundary stencil
   semantics, and whether the sanitizer was silenced rather than the bug
   fixed. Report findings only."
Findings → `results/REVIEW.md`. Human adjudicates.

## Closeout

- [ ] Oracle passes; `FULL=1 status.sh` green solo
- [ ] `results/RESULT.md`: what changed, evidence, deviations
- [ ] changelog/ entry (append-only)
- [ ] `touch results/DONE`
- [ ] Session log line → `docs/llm/SESSION_LOG.tsv`
