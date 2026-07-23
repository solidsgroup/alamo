# TASK: face-consistent-output-stress
# Folder: docs/agent_plans/20260722-face-consistent-output-stress/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 2 (shared mechanics stress field; solver equations unchanged)|
| Model        | Codex root session                                           |
| Verification | partial-oracle                                               |
| Est. scope   | 2-4 src/test files, focused 2-D rebuild and rod/tube run      |
| Parallel-safe| no: shared mechanics and nonlinear-solver headers            |

## Operating rules

1. Read only the files listed in the context budget.
2. Preserve all pre-existing dirty-worktree changes and do not commit them.
3. Do not alter the displacement solve, nonlinear residual, material update,
   or pressure loading.
4. Keep legacy nodal stress behavior for non-conservative elasticity paths.
5. Stop if the focused baseline displacement or face-stress oracle regresses.

## Context budget

Read first: `benchmark/status.sh` output, `docs/llm/PLAN.md`, this `PLAN.md`
Read: `src/Integrator/Base/Mechanics.H:45-280`,
      `src/Solver/Nonlocal/Newton.H:70-215,280-380,1020-1070`,
      `src/Numeric/Stencil.H:830-925`,
      `src/Integrator/Integrator.H`, `src/Integrator/Integrator.cpp`,
      `src/Test/Numeric/Stencil.H`, `src/Test/Solver/Nonlocal/Newton.H`,
      `src/test.cc`, `Makefile`
Reference: `input_rt1s_ideal`,
           `docs/agent_plans/20260722-rod-tube-resolution-phi-study/analyze_interface_stress.py`,
           `output_rt1s_ideal_ncell64_casingAl_void0.5_0.5/05000node`
Forbidden: `docs/archive/*`, unrelated task folders

## Objective

Make the stress field written by conservative elasticity runs consistent with
the face stress used by the solve. Reconstruct each nodal tensor column from
the adjacent conservative face tractions while preserving the existing nodal
stress calculation for legacy/masked elasticity. Verify that the phi-boundary
line is reduced without changing displacement, chamber loading, or the solver.

## Oracle

Commands: focused 2-D build/test; one-second rod/tube run at `np=4`; independent
post-processing of its displacement/model fields; `benchmark/status.sh`.
Covers: compilation, legacy helper tests, finite output, conservative
face/output agreement, unchanged displacement to numerical tolerance, and the
phi-interface mismatch reduction.
Does NOT cover: 3-D GPU runtime performance or every stateful mechanics model.

## Steps

### Step 1 - Establish the output boundary and regression oracle
VERIFY: identify all consumers of `stress_mf`, registration/output timing, and
the conservative-face mode query without modifying source.
DO: choose the smallest implementation that cannot feed a changed plotted
stress back into the solve; add or adapt a focused face-stress reconstruction
test if feasible.
CHECK: document exact old/new semantics and expected boundary behavior.

### Step 2 - Implement conservative face-consistent stress
VERIFY: legacy/masked path remains selectable and face stencils have sufficient
ghost support at patch and physical boundaries.
DO: reconstruct tensor column `d` from the average of low/high face stresses in
direction `d`. Retain the established nodal value at physical boundaries,
where the centered tangential face stencil extends outside the domain. Share
the existing face-model/face-gradient calculation rather than duplicating
physics.
CHECK: focused build/test passes; diff confirms no residual or solve changes.

### Step 3 - Run rod/tube output oracle
VERIFY: use a unique output path and the completed baseline as the displacement
reference.
DO: run to one second at `np=4`, generate `stress_validate` and the independent
face diagnostic.
CHECK: run exits zero, output is finite, displacement remains within numerical
tolerance, and stored nodal stress agrees with independently reconstructed
face stress through the phi boundary.

## Checkpoints

- [x] After plan restatement: human confirms the output-only approach before source edits
- [x] Before completion: diff summary and oracle output reviewed

## Closeout

- [x] Oracle passes; status.sh device lint remains green
- [x] `results/RESULT.md` records implementation, evidence, and limitations
- [x] `results/DONE` created
- [x] One outcome line appended to `docs/llm/SESSION_LOG.tsv`
