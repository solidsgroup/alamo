# Full GCC regression review — 2026-07-13

## Decision

**BLOCKED.** The focused finite-soft-void tests pass, but the paired-stencil
work is not ready to commit or merge. Do not refresh references, relax
tolerances, increase iteration caps, or add an artificial psi/coefficient
floor to hide the failures.

The detailed executable handoff is REGRESSION_RECOVERY_HANDOFF.md. The
follow-on Tier-3 task is
../../20260713-elastic-regression-recovery/PLAN.md.

## Execution record

The worktree was /tmp/alamo-elastic-void-continue on branch
elastic-void-wip-20260713 at a06c6f15d, with the source delta uncommitted.
Both 2-D and 3-D GCC targets were built, then this complete runner was used:

~~~
scripts/runtests.py --comp=g++ --no-backspace
~~~

Canonical artifacts:

~~~
report/output_2026-07-13_15.29.01_kermit.json
report/output_2026-07-13_15.29.01_kermit.html
~~~

All 144 sections executed, with no skips:

| Result | Count |
|---|---:|
| Run pass | 131 |
| Run fail | 13 |
| Check pass | 85 |
| Check fail | 17 |
| No check | 42 |

All Unit sections passed. ElasticSoftVoid 2-D and 3-D also passed. That
focused result is evidence of a narrow viable configuration, not proof that
the revised elastic operator preserves the existing solver surface.

The aggregate make test target was attempted. Its tab check passed, but the
documentation phase stopped because local doxygen cannot load libclang-18.so.18.
This environment issue is separate from the test-runner failures.

## Runtime failures

| Failure mode | Affected sections |
|---|---|
| Iteration cap | Eshelby / 3D-serial-4levels: 20, 2.097364563e-08; Eshelby / 3D-parallel-5levels: 20, 4.256658658e-08; EshelbyFiniteKinematics / 2D-serial: 20, 1.036873858e-08; PlateHole / 2D-serial: 20, 2.023381522e-06; PlateHole / 3D-serial: 20, 2.803395211e-06; PlateHole / 3D-parallel: 20, 2.912295376e-06; Scratch / 2D-serial: 40, 1.952953903e-08; Scratch / 2D-parallel: 40, 1.952953682e-08. |
| Divergence/stall | FracturePFCZM / 2d-serial-notch: 404, 1.057125712e+20; RubberPlateHole / serial-2d: t=0.5, 150, 0.241361769; RubberPressurizedHole / 2d-serial: t=2, 21, 3.536749885e+20; RubberPressurizedHole / 2d-parallel: t=2, 21, 4.119987487e+20; TopOp / parallel: step 58, 88, 1.273287521e+20 and MPI_ABORT. |

## Numerical-check failures

These are genuine output mismatches, not artifacts from a dirty worktree:

- RubberWithInclusion / serial-2d
- SCPSpheresElastic / 2d-parallel-long
- Suture / 2D-serial-4levels
- ThermoElastic / 2d
- TwinGrowth / serial
- VoronoiElastic / 2d-serial-restart, 2d-serial, 2d-parallel,
  periodic-2d-serial
- VoronoiSimplePeriodic / 2d-amr0-mg0-np1, 2d-amr2-mg0-np1,
  2d-amr0-mga-np1, 2d-amr2-mga-np1, 2d-amr0-mg0-np4,
  2d-amr2-mg0-np4, 2d-amr0-mga-np4, 2d-amr2-mga-np4

SCPSpheresElastic is the first priority. At t=0.006 / step 6000, the linear
residual drops only to about 9.6e-03, rises to 5.77e3 by the fixed iteration
limit, then the next RHS/residual reaches O(1e23). Plot 10000 onward contains
NaN displacement, stress, and residual fields, while eta/psi/model remain
finite. The current checker accepted this false success.

The other checks are finite mechanics drift. Examples include:

- RubberWithInclusion stress error 8.33e-03 versus 1e-04.
- ThermoElastic displacement/stress about 0.038 versus 0.01.
- Suture stress about 0.095553 versus 1e-04.
- VoronoiElastic displacement about 0.867–0.873 and stress about
  0.790–0.823 versus 0.01.
- VoronoiSimplePeriodic stress about 0.0668–0.0674 versus 0.05.

Retained logs include:

~~~
tests/PlateHole/output_2026-07-13_15.29.01_kermit_2D-serial/stdout
tests/SCPSpheresElastic/output_2026-07-13_15.29.01_kermit_2d-parallel-long/stdout
~~~

Every report entry maps to its corresponding retained output directory.

## Source-review findings

1. **SetPsi halo contract.** Elastic.cpp allocates three psi ghosts for the
   paired stencil. The public per-level SetPsi overload copies the caller tile
   but does not invoke FillPsiGhosts(), while the all-level Field overload
   does. Fapply can read the third low collar through a predecessor flux
   anchor. Normal vector Newton happens to use the safe overload, but direct
   callers can get default or stale collar data. No in-tree direct caller was
   found; the public API still needs a contract or a test.
2. **Masked cross-AMR effective coefficient.** With
   elastic.solver.average_down_coeffs=1, DDW is restricted separately from
   psi/theta even though Fapply multiplies them. Masked C/F interfaces may
   therefore use inconsistent theta*DDW. FracturePFCZM is the targeted
   reproducer.
3. **Direct Apply ghost precondition.** Newton checks two displacement/model
   ghosts before paired-stencil use. Direct callers have no comparable guard,
   and PairedGradient cannot make an Array4 anchor outside its allocation
   safe. CUDA was not rerun after the last source edit.
4. **Exact-zero psi remains unsupported.** The zero-diagonal diagnostic is
   not an inactive-DOF/nullspace policy. It must not become a hidden floor.

## Oracle-review findings

ElasticSoftVoid is not a sufficient acceptance test:

- The 2-D runner overrides the void to 0.2 MPa, but the checker merely
  requires contrast at least ten; the base material already meets it.
- It lacks an automated MPI section, proof that heterogeneity crosses a C/F
  interface, all-component validation, and a mechanics-sensitive reference
  such as force balance, symmetry, traction, or a manufactured solution.
- It tests finite output/residual metadata but not the late-time NaN/fixed
  iteration path exposed by SCPSpheresElastic.
- Existing unit coverage lacks direct tests for PairedGradient, stencil bounds,
  bounded interpolation, diagonal impulse agreement, and exact-zero rejection.

## Attribution boundary

Only PlateHole 2-D is confirmed pre-existing relative to the day's unstaged
continuation: it failed before the current source timestamps and showed the
same 20-iteration-cap class in the main checkout. Do not call the remaining
failures pre-existing or caused by the latest delta without a clean same-base
comparison.

Historical June 29 reports in the main checkout passed many affected families,
including PlateHole, Eshelby, Rubber, Scratch, Suture, ThermoElastic,
TwinGrowth, and Voronoi. The WIP stack is therefore a regression relative to
that history, but individual attribution still needs bisection.

## Required next action

Start with clean same-base baselines, then reduce SCPSpheresElastic at the
first bad solve and compare the residual, Fapply, Diagonal, boundary closure,
and coefficient hierarchy. Use PlateHole as the fast deterministic
linear-solver signal, then test masked AMR/MPI separately. The full ordered
plan, reproducer commands, guardrails, and acceptance gates are in
REGRESSION_RECOVERY_HANDOFF.md.
