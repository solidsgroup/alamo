# Face-consistent output stress

Date: 2026-07-22

## Outcome

Conservative elasticity runs now write a stress tensor reconstructed from the
same directional face tractions used by the nonlinear solve. In the one-second
64 x 64 / max-level-2 rod-and-tube case, the maximum radial-profile mismatch
through the diffuse `phi` interface fell from 0.2108 MPa to 0.0207 MPa, a
90.2% reduction. Visual inspection shows that the localized old nodal bump and
undershoot are absent from the new curve; the new output lies nearly on top of
the independently recovered face stress.

The displacement solve and all evolving fields are unchanged. A paired
four-rank 0.02-second run with stress plotting enabled and disabled produced
exactly zero maximum absolute difference in displacement, RHS, strain,
material parameters, `phi`, `eta`, `psi`, and temperature.

## Implementation

- `Solver::Nonlocal::Detail::NodalStressFromFaces` reconstructs tensor column
  `d` by averaging the low and high face-stress columns in direction `d`. It
  reuses `FaceStress`, and therefore exactly shares the conservative face
  model interpolation, face gradient, kinematic-variable dispatch, and
  constitutive `DW` evaluation.
- `Mechanics` retains `stress_mf` as the evolving internal nodal constitutive
  stress. Material advancement and integrated diagnostics continue to consume
  this established field.
- When `elastic.plot_stress=1`, a separate zero-ghost, non-evolving output
  field is registered under the existing `stress_*` plot names. Conservative,
  unmasked solves fill it with the face reconstruction. Masked/nonconservative
  paths preserve the legacy nodal output.
- Physical-domain boundary nodes retain the legacy nodal value because the
  centered tangential part of the face-gradient stencil is not defined beyond
  the physical domain. The reported `phi` interface is internal and uses the
  new reconstruction.
- Dynamic mechanics copies its established internal stress to the output
  field, preserving its previous semantics.

## Verification

### Build and unit tests

`make -j4 bin/test-2d-g++ bin/alamo-2d-g++` completed and linked both targets.
`./bin/test-2d-g++` exited zero with `0 tests failed`. The new
`conservative nodal stress reconstruction` test passed. `git diff --check`
also passed. This environment did not provide `clang-format` or `eclint`.

### Output-only paired oracle

Two 0.02-second runs used the same executable and four MPI ranks, differing
only in `elastic.plot_stress`:

- `output_face_stress_on_short`
- `output_face_stress_off_short`

All 18 common plotted components compared exactly, with maximum absolute
difference 0.0. The complete component report is
`short_solution_comparison.json`.

### Full rod-and-tube oracle

The new four-rank run
`output_rt1s_face_consistent_ncell64_maxlevel2` reached step 5000 / `t=1.0 s`,
finalized normally, and reports `Status = Complete (99%)`. Runtime was
493.658 seconds.

| Metric | Legacy nodal output | Face-consistent output |
|---|---:|---:|
| Propellant-side max (output - independently recovered face) | +0.2108 MPa | +0.0207 MPa |
| Interface minimum (output - independently recovered face) | -0.0772 MPa | -0.00618 MPa |
| Interface maximum absolute mismatch | 0.2108 MPa | 0.0207 MPa |
| `stress_validate` radial RMS difference | 0.139 MPa | 0.1288 MPa |
| `stress_validate` hoop RMS difference | 0.184 MPa | 0.1848 MPa |

The final analytical validation retains the prior global behavior: traction
transmission is 0.9246, and the tube peak-hoop ratio is 0.464. The change is
therefore localized to output reconstruction rather than equilibrium or
loading.

`benchmark/status.sh` reports `device-lint: PASS` and `a100-sanitizer: PASS`.
Its golden comparison currently fails while compiling an unrelated dirty-tree
change in `src/Operator/Elastic.cpp` (`AMREX_D_TERM` receives too many macro
arguments); this task does not modify that file.

## Limitations

- The output field adds storage for one nodal stress tensor on each active AMR
  level whenever stress plotting is enabled, although it carries no ghosts.
- Physical-domain boundary nodes deliberately retain nodal stress.
- In this driver, plotfiles are written after the evolution step but the
  static elastic solve occurs at the beginning of a mechanics interval. Thus
  the saved temperature/material fields can be up to one mechanics interval
  newer than the stress and displacement. An independent reconstruction using
  the final saved model therefore retains a small temporal mismatch; it is not
  the removed nodal/face spatial reconstruction error.
- The focused runtime/build oracle was two-dimensional CPU. The repository's
  static device lint and existing A100 sanitizer pass, but a fresh 3-D/CUDA
  build was not run for this task.

## Artifacts

- Interface comparison: [interface_stress_composite.png](interface/interface_stress_composite.png)
- Interface metrics: [interface_metrics.json](interface/interface_metrics.json)
- Analytical validation: [compare_alamo.png](stress_validate/compare_alamo.png)
- Analytical CSV: [compare_alamo.csv](stress_validate/compare_alamo.csv)
- Paired-field comparison: [short_solution_comparison.json](../short_solution_comparison.json)
- Full run log: [full.log](../full.log)
- Unit-test log: [unit.log](../unit.log)
- Reproducible plot-component comparator: [compare_components.py](../compare_components.py)

