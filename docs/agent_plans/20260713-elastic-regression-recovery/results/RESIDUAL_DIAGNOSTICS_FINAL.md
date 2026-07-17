# Elastic soft-void residual diagnostics and recovery

## Outcome

The 0.2 MPa, unmasked (`elastic.use_psi=0`) ElasticSoftVoid problem now
converges to the unchanged force-balance oracle in both 2-D and 3-D.  The
repair makes the nonlinear residual and tangent use the same conservative
face-stress discretization.  It does not add a psi or stiffness floor, relax a
tolerance, increase an iteration cap, force Newton acceptance, or refresh a
reference.

The generic Mechanics path, including `m_psi`/masked solves, retains the exact
trusted-parent algebra and ghost-row behavior.  Flame selects the new
conservative face-flux formulation only for its unmasked (`use_psi=0`) solve.
Binding the AMR policy to the formulation, rather than merely to the presence
of `m_psi`, is required by both the legacy elasticity oracles and the unchanged
soft-void oracle.

| Focused current-source gate | Result | Measured behavior |
|---|---|---|
| 2-D unit executable | PASS | zero failed tests, including face-stencil and Newton decision tests |
| 3-D unit executable | PASS | zero failed tests, including 3-D face-stencil identities |
| ElasticSoftVoid / 2-D serial | PASS/PASS | MLMG 37/45; final nonlinear ratio `1.91078e-6` |
| ElasticSoftVoid / 2-D, 2 MPI ranks | PASS/PASS | final nonlinear ratio `1.91153e-6` |
| ElasticSoftVoid / 3-D serial | PASS/PASS | MLMG 37/35; final nonlinear ratio `1.46238e-6` |
| PlateHole / 2-D serial | PASS/PASS | trusted-parent result restored: 14 MLMG iterations, `4.756775594e-7` |
| RubberPlateHole / 2-D serial | PASS | nonlinear run completes |
| SCPSpheresElastic / 2-D short serial | PASS/PASS | finite mechanics and numerical checker pass |
| Complete GCC regression suite | PASS | 145 runs, 113 run-and-verified, zero failures |

The exact-current-source focused runs were made on 2026-07-16.  In
particular, the final 2-D pair is recorded by test ID
`output_2026-07-16_23.29.02_kermit`, the final 3-D oracle is
`tests/ElasticSoftVoid/output_2026-07-16_23.25.36_kermit_3d-serial`, and the
requested PlateHole/RubberPlateHole/SCPSpheresElastic gate is recorded by test
ID `output_2026-07-16_23.30.02_kermit`.

## Retained failure localization

The retained failing plots were read without changing them.  Both maxima lie
in the diffuse material transition on an outer row of the level-1 fine-grid
union.  They are coarse/fine-interface rows, not physical-boundary rows or
ordinary interior rows.

| Retained plot | Maximum residual | AMR level and index | Coordinate (m) | Classification | `max|res|/max|rhs|` |
|---|---:|---|---|---|---:|
| 2-D | `res_x = +7.083374394e7` | L1 `(64,24,0)` in a max-L2 hierarchy | `(0.1315500,0.0767375)` | diffuse interface; outer C/F-union row | `0.270766694` |
| 3-D | `res_y = -7.067038657e7` | L1 `(40,0,6)` in a max-L1 hierarchy | `(0.0986625,0.0438500,0.0164625)` | diffuse interface; outer C/F-union row | `0.270142249` |

The common RHS maximum is `2.616043465e8`.  The corresponding final nonlinear
ratios were `0.236769` in 2-D and `0.236027` in 3-D even though MLMG converged
the individual linear corrections.  The maxima identify where AMR amplified
the defect; they do not establish AMR as its cause.  Old single-level hard
controls still had plotted ratios near `0.23425`.

## Controlled matrix after the repair

The expanded diagnostic harness is `benchmark/elastic_void_matrix.py`.  It
does not weaken or replace `tests/ElasticSoftVoid/test`; it stages independent
inputs and products under `/tmp`, captures every Newton and MLMG diagnostic,
and classifies the plotted residual maximum.  The final complete matrix is
`/tmp/alamo-elastic-faceflux-matrix-20260716-dispatch-final/matrix.json`
(SHA-256
`f38129b547c86cc34fabe952c1b714566844acd3d27d09755a81ce505b426771`).
All 16 combinations completed and every post-repair residual maximum was an
interior node.

`Soft` means 0.2 MPa void bulk and shear modulus; `hard` means 10 MPa.  The
Newton column lists MLMG iteration counts for successive corrections and the
final composite nonlinear ratio.  `Plot ratio` is the independently measured
true plotted force imbalance.

| Dim. | Void | Thermal | AMR | MLMG iterations | Final nonlinear ratio | Plot ratio | Maximum residual field/value |
|---:|---|---|---|---|---:|---:|---:|
| 2 | soft | on | on | 40 | `3.86823e-3` | `6.45181e-4` | `res_x=1.68782e5` |
| 2 | soft | on | off | 5 | `3.09765e-3` | `7.06118e-4` | `res_x=1.52990e5` |
| 2 | soft | off | on | 37, 45 | `1.91078e-6` | `1.91078e-6` | `res=4.99868e2` |
| 2 | soft | off | off | 5, 5 | `8.76643e-7` | `8.76643e-7` | `res=1.89937e2` |
| 2 | hard | on | on | 40 | `2.95048e-3` | `5.95123e-4` | `res=1.55687e5` |
| 2 | hard | on | off | 4 | `2.87855e-3` | `8.03837e-4` | `res=1.74162e5` |
| 2 | hard | off | on | 36, 37 | `1.99500e-7` | `1.99500e-7` | `res=5.21900e1` |
| 2 | hard | off | off | 5, 5 | `2.52131e-7` | `2.52131e-7` | `res=5.46277e1` |
| 3 | soft | on | on | 36, 39 | `1.35340e-5` | `3.11282e-6` | `res=8.14327e2` |
| 3 | soft | on | off | 6, 25 | `1.14637e-5` | `2.91897e-6` | `res=6.32436e2` |
| 3 | soft | off | on | 37, 35 | `1.46238e-6` | `1.46238e-6` | `res=3.82565e2` |
| 3 | soft | off | off | 8, 24 | `9.00337e-7` | `9.00337e-7` | `res=1.95070e2` |
| 3 | hard | on | on | 39 | `3.25792e-3` | `6.22511e-4` | `res=1.62852e5` |
| 3 | hard | on | off | 5 | `3.31411e-3` | `7.47297e-4` | `res=1.61912e5` |
| 3 | hard | off | on | 37, 35 | `2.04645e-7` | `2.04645e-7` | `res=5.35360e1` |
| 3 | hard | off | off | 5, 11 | `2.63860e-7` | `2.63860e-7` | `res=5.71689e1` |

Every reported MLMG solve reached a relative linear residual between
`1.48e-6` and `9.997e-6`.  Before the repair this fact was misleading: MLMG
was accurately solving a tangent equation that was not the derivative of the
reported nonlinear residual.  After the repair, Newton uses one or two
corrections in these cases and the independent force-balance oracle passes.
The larger thermal-on nonlinear ratios are normalized by the initial
composite residual, while the plotted ratio is normalized by the plotted RHS;
the denominators are intentionally different.

Two 20-step 2-D soft runs also completed:

| Loading | Final nonlinear ratio | Final plot ratio | Final MLMG iterations |
|---|---:|---:|---:|
| controlled thermal/eigenstrain | `3.75676e-3` | `7.00655e-4` | 45 |
| no imposed eigenstrain | `1.91078e-6` | `1.91078e-6` | 45 |

Their records are
`/tmp/alamo-elastic-faceflux-long-thermal-20260716-dispatch-final/matrix.json`
(SHA-256
`535edd8f7b954dc479bf1385fdd2159bad6b50dcc8d8dea87c9e5316c917cb58`)
and
`/tmp/alamo-elastic-faceflux-long-mechanical-20260716-dispatch-final/matrix.json`
(SHA-256
`567b6381eadaca5e95ef74627d33402a72d33a71fbd396faf7d396d25fa57505`).

### What each ratio means

- The nonlinear ratio is the refluxed AMR-composite infinity norm of
  `b - div(P(u))` after an accepted update divided by its initial composite
  norm for that solve.  It measures progress toward the nonlinear discrete
  equilibrium represented by the solver.
- The plot ratio is `max|true residual field| / max|elastic RHS|` over valid
  plotted nodes and vector components.  It is the test's independent global
  force-balance measure.
- The MLMG relative residual applies only to one tangent correction equation.
  A value near `1e-5` says that linear equation was solved to its requested
  accuracy; by itself it says nothing about whether the tangent is the
  derivative of the nonlinear residual.
- A small Newton update is a step-size condition, not a force-balance proof.
  The old solver stopped on small updates with roughly 24% nonlinear and 27%
  plotted imbalance.  The unchanged plot oracle is what rejected that false
  success.

## Root cause and repair

The decisive defect existed even on a uniform single level.  The nonlinear
residual differentiated stress made from a centered displacement gradient,
whereas `Elastic::Fapply` used a one-cell Hessian plus an explicit
coefficient-gradient product-rule term.  With constant one-dimensional
coefficient `C`, these reduce to different operators:

```
dR/du: C (u[i+2] - 2 u[i] + u[i-2]) / (4 h^2)
Fapply: C (u[i+1] - 2 u[i] + u[i-1]) / h^2
```

The required invariant is now enforced on active algebraic rows:

```
J_h(u) v = d R_h(u + epsilon v) / d epsilon at epsilon = 0
```

The conservative Flame path selected by `elastic.use_psi=0` now:

1. computes displacement gradients and constitutive stress at faces;
2. forms the nonlinear residual as a conservative face-stress divergence;
3. forms the tangent from the derivative of the same face constitutive
   evaluation;
4. uses the same physical-boundary closure and AMR composite-row policy in
   residuals, smoothing, normalization, and diagnostics;
5. uses quadratic coarse-to-fine ghost interpolation and retains the coarse
   composite equation row while refluxing fine residual information; and
6. fills periodic plot ghosts and supports invariant periodic projection for
   the 3-D extrusion check.

Manufactured tests cover mixed-polynomial response and rigid/null modes.  The
Newton tests cover line-search acceptance and termination decisions.  The
diagnostic test was expanded without changing its force-balance threshold or
reference data.

Generic Mechanics, including `m_psi_set == true`, deliberately selects the
exact trusted-parent product-rule residual/tangent, linear interpolation, C/F
residual replacement, ghost-row relaxation, coefficient extent, and
one-component tangent preparation.  This restored PlateHole exactly and made
RubberPlateHole and SCPSpheresElastic complete.  Flame explicitly selects the
conservative formulation when `use_psi=0`; that formulation carries its
matching quadratic interpolation, retained shared coarse C/F row, and
valid-row relaxation policy.  The unchanged unmasked 2-D and 3-D soft-void
tests continue to pass.

The first broad compatibility run exposed one additional diagnostic-only
defect: `compResidual` used a field's maximum allocated AMR capacity rather
than its active `finest_level`.  TopOp allocates room for three levels but has
only level 0 active on its first solve, so true-residual plotting dereferenced
the intentionally null level-1 slot after an otherwise converged solve.  Both
residual overloads now allocate, reflux, and copy only active levels.  The
unchanged TopOp serial, 4-rank 100-step/checker, and coverage sections then all
passed (test ID `output_2026-07-17_00.21.20_kermit`).

## Assessment of the publication error report

### High-contrast items

| Item | Assessment against source and experiments |
|---|---|
| 1.1 C/F stencil inconsistency | **Plausible and confirmed as an amplifier, not the root cause.** Linear coarse-to-fine ghost filling is not second-derivative consistent at the interface, and the failing maxima localized there.  However, old no-AMR controls retained about 23.4% force imbalance.  Quadratic filling improved the 3-D soft/no-thermal pair from 55/53 to 37/35 MLMG iterations, while an attempted blanket `2/3` C/F row scaling worsened convergence and was reverted.  The conservative formulation therefore carries quadratic interpolation plus its consistent composite-row policy; the generic trusted-parent formulation, masked or unmasked, retains its prior AMR policy.  Applying the conservative AMR bundle globally caused broad regression failures, while forcing the trusted-parent bundle onto conservative Flame made 3-D soft void fail its unchanged time gate. |
| 1.2 non-divergence product-rule split | **Strongest confirmed cause, in a more precise form.** The old nonlinear residual and tangent were different discrete operators even for constant coefficients and no AMR.  Conservative face stress plus its exact face tangent repaired all selected cases.  This is stronger evidence than a generic loss-of-diagonal-dominance hypothesis. |
| 1.3 weak smoother | **Secondary, and the report's implementation detail is stale.** Current `Fsmooth` is weighted Jacobi with default `omega=2/3`; the soft-void input requests four pre- and four post-smoothing calls.  MLMG converged every post-repair matrix solve.  Smoother work may improve cost or extreme-contrast margin, but cannot repair a residual/Jacobian mismatch. |
| 1.4 coarsening below interface resolution | **Plausible secondary robustness issue, not causal here.** Prior maximum-coarsening caps, bottom-solver changes, and extra smoothing did not restore force balance.  AMR-off failures established a finer-level inconsistency.  A coefficient-aware coarsening sweep remains useful only as a performance/robustness study after algebraic consistency. |
| 1.5 arithmetic mixture | **Present by design; proposed harmonic replacement is not established.** Flame's phase-field material interpolation is arithmetic and is asserted by the regression.  A scalar Reuss average is not automatically correct for a full nonlinear tensor mixture or its free energy.  Harmonic face homogenization could be studied separately, but silently changing the constituent mixture would change physics and was not part of this repair. |

The report's proposed modulus and psi floors were not used.  They may improve
conditioning by changing the physical problem, but they cannot validate the
zero-floor/unmasked formulation requested here.

### Publication-only and expository items

| Item | Assessment |
|---|---|
| 2.1-2.2 ellipsoid indicator and error-function argument | The printed formulas are internally inconsistent; the current ellipse/sphere initial-condition implementations use geometric radii and a surface-centered transition.  These are paper errors, not the elastic crash mechanism. |
| 2.3 strain-energy factor and fracture normalization | The missing `1/2` is a conventional energy-definition error if the displayed quadratic form is intended literally.  The AT normalization concern affects the quantitative meaning of `G_c`; neither explains the soft-void linear-solve mismatch. |
| 2.4 derivative written with respect to time | Typographical: the phase-field variational derivative must be with respect to the order parameter. |
| 2.5 crack-field convention | The prose is reversed relative to the displayed degradation and code convention.  This is expository. |
| 2.6 undefined refinement field `c` | Typographical; the context indicates the phase/crack field. |
| 2.7 relaxation signature | Notational mismatch between declaration and use in the paper; current C++ interfaces determine the actual call signature. |
| 2.8-2.9 nonlocality/restriction notation | The set/codomain statements and restriction/interpolation terminology are malformed, but they do not specify the current implementation closely enough to diagnose this failure. |
| 2.10 2-D figure colorbar | `sigma_33` is inconsistent with the 2-D text/caption and is a figure-label typo. |
| 2.11 elastic versus total strain in energy | Valid consistency question when eigenstrain is present.  It should be checked against the precise constitutive free-energy definition, but it is separate from the now-demonstrated discrete tangent mismatch. |
| 2.12 element-equivalence assertion | Unsupported as printed; it is not an implementation defect by itself. |

## Broad qualification and remaining chamber run

The definitive complete GCC regression suite passed: 145 tests ran, 113 were
run and verified against checkers, and none failed (test ID
`output_2026-07-17_00.24.59_kermit`; report JSON
`report/output_2026-07-17_00.24.59_kermit.json`).  This includes the unchanged
2-D serial/MPI and 3-D ElasticSoftVoid cases, all requested compatibility
cases, deep AMR and periodic controls, and the corrected TopOp diagnostics.

The remaining qualification is the canonical 6.5 s chamber run from
`input_nova_centre_bore` with 0.2 MPa `model_void.kappa` and `model_void.mu`,
`elastic.use_psi=0`, and no positive psi/stiffness floor.  Its result must not
be inferred from the short regression: it contains 65,000 coarse steps,
approximately 1,300 elastic events, AMR level 3, burn evolution, and a
25,000:1 void-to-casing modulus contrast.

No commit was made, no retained reference or prior artifact was refreshed,
and no `DONE` marker was created.
