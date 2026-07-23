# Stress-ripple diagnostic results

## Bottom line

The observed pattern is not one single defect.  The tests separate it into
three effects:

1. The large line produced by recomputing stress from the nodal displacement
   gradient is a postprocessing artifact.  Reconstructing the tensor from the
   conservative face tractions reduces it substantially.
2. The original final plot has a separate AMR lifecycle inconsistency: level 2
   is regridded after the last elastic solve.  Displacement, model, strain, and
   saved output stress are interpolated independently onto the new hierarchy,
   so the plotted state is no longer the state that converged.
3. A smaller feature at each `eta` interface remains in the raw face traction.
   It is symmetric between x and y, persists several nodes away from the
   symmetry planes, and survives a solve performed on the final hierarchy.
   That part is a discrete interface/load feature, not solely an output error.

Solver convergence does **not** imply that a stress field must be visually
smooth.  It only checks the discrete equation (and this run terminates on the
Newton update, not directly on the nonlinear residual).  A sharp face stress
can balance `rhs = -P grad(eta)`, or adjacent variations can cancel under the
discrete divergence.  A post-solve regrid is also invisible to a solver that
has already returned.

## Tests, in order

| Test | Result | Interpretation |
|---|---:|---|
| Output off/on, same short run | Every common field at steps 0, 50, and 100 was bitwise identical | Stress plotting does not feed back into the solution. |
| Offline reproduction of the face stencil | Levels 0 and 1 agreed with saved stress to about `2e-5 Pa`; original level 2 differed by `96.96 kPa` overall and `48.60 kPa` in the outer-interface band | The output formula is reproducible; only the finest level had become stale/inconsistent. |
| Saved strain versus gradient of saved displacement | 5,496 level-2 values differed above `1e-12`; maximum difference was `0.0195` | Independently remapped derived fields, not constitutive physics, caused part of the final-plot inconsistency. |
| Original final hierarchy equilibrium | Interior `Linf(div(P)-rhs) = 815.7 MPa/m`, while the last solve logged `2.322 MPa/m` | The plotted hierarchy is not the hierarchy on which the residual was measured. |
| Solve on the final/non-regridded hierarchy | Saved versus reconstructed stress returned to `2.60e-5 Pa`; interior residual was `2.32284 MPa/m`, matching the logged `2.32284 MPa/m` | This is the decisive confirmation of the post-solve-regrid defect. |
| Uniform affine manufactured field | Stress variation was at roundoff (`8.24e-5 Pa`) | The face stencil exactly preserves the affine case. |
| Radial diffuse manufactured field | x/y transpose error was `4.70e-5 Pa` maximum | The stencil has no intrinsic x-versus-y bias. |
| Axis layers 0 through 4 | The `eta`-interface roughness remained about `135-140 kPa` on all layers; x/y counterparts agreed closely | It is not caused by a one-row symmetry-boundary stencil. |
| Raw face traction versus nodal reconstruction | At `eta`, face roughness was `135-140 kPa` versus `285-304 kPa` nodally.  At the casing interface, face roughness was `15-50 kPa` versus roughly `1.2-1.3 MPa` nodally | The nodal method strongly amplifies the feature, but the smaller `eta` feature exists before nodal reconstruction. |
| Resolved analytic circular `phi` | Outer-interface face roughness fell from `15-50 kPa` to `1.4-5.1 kPa`; `phi` transpose RMS fell from `3.36e-2` to roundoff.  The `eta` features stayed at `131-139 kPa` | Bitmap/circle anisotropy explains much of the casing-boundary variation, but not the `eta` feature. |
| Legacy product-rule operator | It aborted in MLMG on the same unmasked high-contrast case | This A/B comparison is inconclusive for ripple shape; the legacy operator cannot solve this case as configured. |

The apparent component localization is consistent with radial projection.  On
the horizontal centerline, radial normal stress is `stress_xx`; on the vertical
centerline it is `stress_yy`.  Seeing the pattern in those respective
components is therefore not, by itself, evidence of an x/y indexing error.

The face-reconstructed tensor also need not be exactly symmetric at a node:
each column is assembled from a different staggered face family.  The
manufactured radial test produced up to `20.2 kPa` of nodal `Pxy-Pyx` even while
its x/y transpose symmetry was at roundoff.  Native face traction, rather than
the symmetry of an assembled nodal first-Piola tensor, is the relevant quantity
for checking the conservative solve.

## Recommended next actions

1. Fix the mechanics/regrid ordering.  If a regrid occurs after a static
   solve, solve mechanics again on the new hierarchy before writing a plot.
   Merely recomputing stress is insufficient because the interpolated
   displacement does not satisfy the new hierarchy's discrete equation.
2. Treat strain and stress as derived output.  Recompute them from the final
   displacement and model instead of independently prolongating stale derived
   fields.  Add a plot-time check that saved stress matches face reconstruction
   and that the interior face residual matches the solver's reported residual.
3. Keep the symmetry reconstruction explicitly enabled for this quarter deck:
   `elastic.output_stress_symmetry.xhi = 1` and
   `elastic.output_stress_symmetry.ylo = 1`.
4. Prefer an analytic circular `phi`, or ensure its diffuse width spans enough
   finest-grid cells.  The analytic-circle control materially reduced the
   outer casing artifact.
5. Isolate the remaining `eta` feature with a fixed, uniform hierarchy.  Sweep
   `pf.eps/dx` while holding pressure and material properties fixed, and plot
   the two native face-traction columns together with `-P grad(eta)`.  Repeat
   with analytic circular `eta`.  This distinguishes diffuse-interface
   discretization from bitmap and AMR coarse/fine effects.
6. As a convergence-control experiment, require both update and residual
   tolerances (or run in residual mode) and verify whether the `eta` feature is
   unchanged.  A smaller residual may improve equilibrium accuracy, but it is
   not expected to remove a balanced discrete interface feature.

## Artifacts

- Reproducible analysis: `../check_stress_ripple.py`
- Original frozen-field metrics: `ripple_checks/metrics.json`
- Synchronized final-hierarchy metrics: `ripple_no_postsolve_regrid/metrics.json`
- Analytic-circle metrics: `ripple_exact_phi_fixedgrid/metrics.json`
- Synchronized spatial plot: `ripple_no_postsolve_regrid/spatial.png`
- Synchronized analytical comparison:
  `ripple_no_postsolve_regrid/analytical/compare_alamo.png`
