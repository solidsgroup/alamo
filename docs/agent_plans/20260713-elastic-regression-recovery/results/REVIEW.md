# Elastic regression recovery review

## Decision

**BLOCK.**  The current delta contains valuable, focused repairs, but it is not
release-qualified.  Do not commit, refresh references, add `DONE`, or begin
performance claims until the force-balance failure and the remaining
validation gates are closed.

## Findings

### P0 — update convergence is not equilibrium convergence

The same deterministic chamber stops after two Newton iterations in 2-D
serial, 2-D MPI, and 3-D serial because the maximum update is below `2e-5`.
At that point the nonlinear residual is still about `0.236` of its initial
value and the plotted equilibrium ratio is `0.2701–0.2708`, above the `0.05`
oracle.  A stricter positive update tolerance uses all five iterations and
still reaches only `0.139951` before correctly aborting.

The constituent model is linear.  Before adding a residual convergence mode,
verify that Newton residual assembly, `Fapply`, physical boundary closure,
coarse/fine closure, and the coefficient hierarchy describe the same linear
operator.  A residual stopping mode would currently only expose exhaustion;
it is not itself the repair.

### P1 — final regression surface is incomplete

The trusted-parent restoration passed focused PlateHole,
FracturePFCZM, EshelbyFiniteKinematics, VoronoiSimplePeriodic,
SCPSpheresElastic-short, and earlier 2-D soft-void runs.  A discovery suite
before the final halo changes reached 121 runs/96 checks with only the
contradictory RubberPlateHole coverage section failing; that input was fixed.
However, no complete 2-D or 3-D suite was run on the exact final source after
the Flame/diagonal/normalize changes.  Those earlier results are supporting
evidence, not the final gate.

### P1 — CUDA correctness and performance are absent

No strict- or fast-CUDA correctness case, golden comparison, sanitizer run, or
performance repetition exists for this final delta.  The RTX A1000 is visible,
but `nvcc` is unavailable on `PATH` and was not found under `/usr/local` or
`/opt`.  CPU/GPU speedups must not be inferred from the CPU solver timings.

### P1 — no independent final Tier-3 review

Two read-only closeout reviews were requested, but the subagent quota was
exhausted.  This file is the implementation agent's self-review and does not
satisfy the plan's fresh adversarial-review gate.

### P2 — focused oracle limitations remain explicit

The new checker is substantially stronger: it proves the actual plotted
moduli are positive and sub-MPa, checks the exact constitutive mixture on each
AMR level, distinguishes outer fine-union faces from same-level chopped-box
faces, checks serial/MPI diagnostics, and includes a force-balance gate.
The invalid requirement that one deterministic 2-D solve must backtrack was
removed; backtracking is conditional behavior, while the acceptance predicate
is directly unit tested.

Remaining coverage gaps are an end-to-end failed-line-search rollback case and
a trusted mechanics-sensitive spatial/norm witness for `ElasticSoftVoid`.
Existing mechanics regressions provide broader protection, but they do not
replace these direct cases.  The 3-D symmetry assertions occur after the
equilibrium assertion, so the current failing run was also checked manually:
it is z invariant and has negligible out-of-plane response.

### P2 — halo repairs need final device scrutiny

Building the model on the entire grown nodal box is consistent with Flame's
two-model/three-source ghost allocation and fixes mixed periodic/coarse-fine
corners.  It also relies on physical-boundary source ghosts having been filled
before `UpdateModel`; periodic `FillBoundary` alone does not create nonperiodic
physical data.  Existing CPU cases are finite, but a final reviewer should
trace that lifecycle and validate the complete path on CUDA.  Likewise,
`Elastic::Diagonal` now supplies all two-ring entries read by the smoother,
while `normalize` uses valid rows; this pairing should be retained and reviewed
as one contract.

## Guardrail audit

| Guardrail | Current state |
|---|---|
| No psi floor | Satisfied: `elastic.psi_floor=0` |
| No hidden stiff void | Satisfied: plotted minima are `0.64–0.76 MPa` |
| No relaxed linear tolerance | Satisfied: final linear relative residuals are below `1.0e-5` |
| No increased nonlinear cap | Satisfied: cap remains five |
| No forced line-search acceptance | Satisfied; failed search restores and aborts |
| No reference refresh | Satisfied |
| No false completion record | Satisfied: blocked, no `DONE`, no commit |

## Required disposition

Continue from `HANDOFF.md`.  If the residual/tangent experiment proves a
source inconsistency, repair and add a discriminating low-level test before
rerunning the chamber.  If it instead proves the plotted residual oracle uses
a different, justified quantity, document that derivation and change the
oracle only with independent evidence—not by choosing a larger threshold.
