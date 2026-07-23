# Rod-and-tube resolution and phi-interface study

Date: 2026-07-22

## Outcome

All three new four-rank simulations reached step 5000 / `t=1.0 s` and exited
zero. The stress feature at the propellant side of the casing `phi` transition
is present at both requested resolutions. It does not converge away when the
finest spacing is reduced from 0.3426 mm to 0.2284 mm.

The feature is principally a plotted nodal-stress reconstruction error. A
spatial AMR composite of the stored nodal stress differs locally from stress
recovered using the conservative face gradients and face-averaged material
model used by the elastic solve. The recovered face stress is smooth through
the location of the visible nodal feature.

Replacing the broad bitmap mask with an analytic circular mask whose tanh
scale is `pf.eps = 0.125 mm` does not help at the baseline mesh. It narrows the
measured 1--99% transition from about 11.5 mm to 0.5 mm, only about 1.5 finest
cells, and increases the maximum absolute nodal/face mismatch from 0.211 MPa
to 0.942 MPa.

## Runs

| Case | Base cells | max_level | Finest dx | Wall time | Result |
|---|---:|---:|---:|---:|---|
| bitmap baseline (existing) | 64 x 64 | 2 | 0.3426 mm | existing | complete |
| bitmap requested | 96 x 96 | 2 | 0.2284 mm | 12m41s | complete |
| bitmap requested | 128 x 128 | 1 | 0.3426 mm | 4m22s | complete |
| analytic `phi`, tanh scale=`pf.eps` | 64 x 64 | 2 | 0.3426 mm | 9m14s | complete |

The 96 and 128 cases ran concurrently, each with `mpiexec -np 4`. Their driver
timestamps were 12:41:04--12:53:45 and 12:41:04--12:45:26 CDT,
respectively. Each new `thermo.dat` ends at 0.9998 s and each output contains a
one-second `05000node` plotfile. Metadata records the requested base grids and
AMR levels.

## Interface diagnostic

The custom diagnostic first masks coarse nodes covered by finer AMR patches,
then area-weights the surviving hierarchy in radial bins. It reconstructs the
elastic face stress with the same endpoint-averaged
`NeoHookeanPredeformed` model, `FaceGradient`, and `DW` calculation used by the
conservative solve; adjacent face-traction columns are averaged back to nodes
only for comparison with the plotted nodal tensor.

| Case | Measured phi 1--99% width | Propellant-side max (nodal - face) | Interface max absolute mismatch |
|---|---:|---:|---:|
| 64 / level 2 bitmap | 11.50 mm | +0.211 MPa | 0.211 MPa |
| 96 / level 2 bitmap | 11.25 mm | +0.237 MPa | 0.237 MPa |
| 128 / level 1 bitmap | 11.25 mm | +0.212 MPa | 0.212 MPa |
| 64 / level 2 epsilon mask | 0.50 mm | +0.912 MPa | 0.942 MPa |

The 64/level-2 and 128/level-1 cases have the same finest spacing and nearly
identical mismatch (0.211 versus 0.212 MPa), despite different base grids and
AMR hierarchies. The 96/level-2 case has 1.5 times finer spacing but a slightly
larger mismatch. This rules out a simple base-mesh convergence explanation.

For the epsilon-width arm, the composited nodal radial stress jumps to
-1.820 MPa on the propellant side and -3.674 MPa on the casing side, while the
recovered face values remain near -2.73 MPa at both locations. The alternating
sign and smooth face result identify a reconstruction overshoot rather than a
solver-resolved physical concentration.

## Existing stress_validate results at one second

| Case | radial RMS diff | hoop RMS diff | traction transmission | peak tube hoop ratio |
|---|---:|---:|---:|---:|
| 64 / level 2 bitmap | 0.139 MPa | 0.184 MPa | 0.9245 | 0.467 |
| 96 / level 2 bitmap | 0.164 MPa | 0.178 MPa | 0.9286 | 0.464 |
| 128 / level 1 bitmap | 0.136 MPa | 0.177 MPa | 0.9245 | 0.468 |
| 64 / level 2 epsilon mask | 0.150 MPa | 0.092 MPa | 0.9230 | 0.915 |

The epsilon mask makes the analytical sharp-interface comparison look better
in hoop stress while simultaneously making the plotted local stress spike much
worse. The global analytical metric therefore must not be used as evidence
that the nodal interface stress is trustworthy.

`compare_alamo.py` currently keeps the maximum AMR level found anywhere in an
entire radial bin. That can discard valid coarse-level angular coverage and
exaggerate isolated points. Correct spatial compositing reduces this plotting
bias, but does not eliminate the nodal/face mismatch documented above.

The simulation emits an existing warning that MPI `thermo.dat` traction and
displacement integrals may be incorrect. The analytical validation uses its
pressure value, so those pressure-derived ratios retain that caveat. The
nodal-versus-face interface comparison uses only the one-second plotfile and is
not affected by the thermo integral warning.

## Recommendation

1. Do not use a `phi` tanh scale equal to `pf.eps` at the present mesh. It is
   under-resolved and amplifies the displayed-stress overshoot.
2. Make plotted/output stress consistent with the conservative elasticity
   operator: compute face stress with `FaceGradient` and the face model, then
   expose face traction directly or reconstruct a nodal tensor from adjacent
   face columns for visualization.
3. Fix `stress_validate` AMR reduction to spatially mask covered coarse nodes
   and area-weight the surviving hierarchy, rather than choosing one maximum
   level per radial bin.
4. Treat `phi` as a material-interface resolution choice independent of the
   burn phase-field epsilon. If a diffuse material transition is retained, test
   a width resolved by at least 4--6 finest cells; for the baseline grid that
   means a tanh scale around 0.30--0.45 mm (1--99% width about 1.4--2.1 mm).

## Artifacts

- Final four-case composite: [interface_stress_composite.png](composite_all/interface_stress_composite.png)
- Metrics: [interface_metrics.json](composite_all/interface_metrics.json)
- Reproducible diagnostic: [analyze_interface_stress.py](../analyze_interface_stress.py)
- 96/level-2 validation: [compare_alamo.png](ncell96_maxlevel2/compare_alamo.png) and [sweep](ncell96_maxlevel2/compare_alamo_sweep.png)
- 128/level-1 validation: [compare_alamo.png](ncell128_maxlevel1/compare_alamo.png) and [sweep](ncell128_maxlevel1/compare_alamo_sweep.png)
- Epsilon-mask validation: [compare_alamo.png](ncell64_maxlevel2_phi_eps/compare_alamo.png) and [sweep](ncell64_maxlevel2_phi_eps/compare_alamo_sweep.png)

## Verification

- All three new logs contain `STEP 5000 ends. TIME = 1` followed by AMReX
  finalization.
- All three new output metadata files report `Status = Complete (99%)`.
- All final composite CSV numeric values are finite.
- All four one-second `compare_alamo.py` runs exited zero and wrote PNG/CSV
  outputs; requested new cases also wrote time-sweep PNG/CSV outputs.
- No production source file was modified.
