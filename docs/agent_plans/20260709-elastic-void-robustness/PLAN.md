# TASK: elastic-void-robustness
# Folder: docs/agent_plans/20260709-elastic-void-robustness/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 3                                                             |
| Model        | opus                                                          |
| Verification | partial-oracle                                                |
| Est. scope   | Elastic operator/solver plus focused 2D/3D regression inputs |
| Parallel-safe| no, shares Elastic/Newton paths                              |

## Operating rules

1. Read ONLY the files listed in Context budget. docs/archive/ is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: stop,
   record the discrepancy, and refine the diagnosis before changing physics.
3. One commit per step unless stated. Message: `<area>: <what> (<task-folder>)`.
4. No scope expansion. New ideas go to NOTES.md in this folder, not into code.
5. Do not weaken the oracle by adding an artificial mask floor or stiff void.
6. Tier 3 checkpoints apply before source edits and before commits.

## Context budget

Read first: benchmark/status.sh output, this PLAN.md
Read: input_rod_and_tube_2d, input_rod_and_tube_3d,
  input_confirm_rod_and_tube, src/Operator/Elastic.{H,cpp},
  src/Integrator/Flame.{H,cpp}, src/Integrator/ThermoElastic.H,
  src/Solver/Nonlocal/Newton.H, src/Solver/MLMG.*,
  Model construction files reached from Flame/ThermoElastic,
  existing non-archived rod-and-tube logs and input decks
Reference only if step names it: arXiv:2001.04789, AMReX MLMG documentation
Forbidden: docs/archive/*

## Objective

Make the finite-difference elastic solve robust when the phase-field mask is
zero in void and the void model is physically soft. Establish the failure
mechanism from an independently repeatable 2D and 3D rod-and-tube case, then
make the smallest discretization/solver change that removes the artificial
`psi_floor` and high void-modulus requirement without changing solid-region
elasticity or silently accepting an unconverged solve.

## Oracle

Command(s): focused 2D and 3D rod-and-tube regression commands created by this
task, run with `elastic.psi_floor = 0` (or omitted) and a soft nonzero void
modulus selected from the material model, plus existing Elastic unit/regression
tests and `benchmark/status.sh`.
Covers: no MLMG/Newton convergence failure, finite displacement/stress/energy,
residual criteria met, and 2D/3D execution through the targeted physical event.
Does NOT cover: independent experimental validation of a void constitutive law;
the numerical formulation must still be reviewed for consistency with the
masked-domain continuum problem.

## Steps

### Step 1 - establish a controlled baseline
VERIFY: the current built binaries and the two rod-and-tube inputs are present.
DO: record solver configuration, exact mask/coefficient construction, and run
the current artificial-stabilization baseline plus the zero-floor/soft-void
counterpart in 2D and 3D.
CHECK: retained logs identify the first failed solve and its MLMG/Newton
residual history, rather than only a later phase-field failure.

Status: completed. See NOTES.md.

### Step 2 - isolate the operator defect
VERIFY: baseline result distinguishes a masked linear-system failure from a
time-integration or constitutive instability.
DO: use small manufactured/geometry-preserving probes to test diagonal,
null-space treatment, coefficient coarsening, and mask behavior at zero psi.
CHECK: each candidate cause makes a falsifiable prediction confirmed or
rejected by a run or operator diagnostic.

Status: completed. The residual and Jacobian use incompatible finite-difference
stencils; see NOTES.md.

### Step 3 - implement the bounded correction
VERIFY: the identified mechanism has a correction consistent with the paper's
diffuse-interface formulation and AMReX operator contract.
DO: modify only the required Elastic/Newton/MLMG code; add a focused regression
for exact-zero mask and soft void in each supported dimension.
CHECK: compile cleanly and all focused oracles pass with no implicit floor.

### Step 4 - adversarial verification
VERIFY: source diff and focused regressions are complete.
DO: test resolution, decomposition, and void-stiffness sensitivity; compare
solid-region stress/displacement to the stabilized reference where meaningful.
CHECK: both dimensions remain converged and finite; no test relies on a hidden
coefficient floor or forced iteration count.

## Checkpoints

- [x] After plan restatement: diagnosis first, no source edit before a
      reproducible failure and mechanism.
- [x] Before source edit: failure signature, proposed invariant, and affected
      path recorded in NOTES.md; awaiting human confirmation.
- [ ] Before commit: diff summary, focused oracle output, and adversarial
      review findings.

## Adversarial review

After implementation, review the final diff specifically for accidental
regularization, changed boundary conditions, loss of elasticity symmetry,
stale MLMG coefficients, tolerance changes, and tests weakened to pass. Write
findings to results/REVIEW.md.

## Closeout

- [ ] Oracle passes; status.sh all green
- [ ] results/RESULT.md states mechanism, validation matrix, and limitations
- [ ] touch results/DONE
- [ ] Session log line appended to docs/llm/SESSION_LOG.tsv
