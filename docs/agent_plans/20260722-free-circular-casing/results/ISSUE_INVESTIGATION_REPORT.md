# Rod-and-tube stress artifact investigation

Date: 2026-07-22

## Executive conclusion

Two visually related but distinct problems were present.

1. The thin alternating stress line on the propellant side of the `phi`
   transition was primarily an output-stress reconstruction error.  Replacing
   the legacy nodal constitutive output with stress reconstructed from the
   conservative face tractions reduced the local mismatch by 90.2%.
2. The remaining loss of axisymmetry near the casing is in the displacement
   solution itself.  The strongest evidence points to `base_circle0.bmp`: its
   `phi=0.5` contour is elliptical by about 1.5 mm and its material transition
   is extremely broad.  The new circular casing-support field is still exactly
   one where the largest polar shear occurs, so it is not the source of that
   peak.

## Technical explanation of the face-consistent output change

### Old plotted stress

The original plot path evaluated one constitutive stress directly at each
node.  For the finite-strain model this was

```text
F_i       = I + grad_node(u)_i
P_old,i   = DW(model_i, F_i)
```

Here `grad_node` is a nodal finite-difference gradient and `model_i` is the
already blended nodal material model.  This tensor was useful as a local
constitutive diagnostic, but it was not the quantity used by the conservative
elastic operator to enforce equilibrium.

### Stress used by the conservative solve

The conservative operator evaluates stress on coordinate faces.  On an
`x`-face between nodes `(i,j)` and `(i+1,j)`, for example, it uses

```text
model_(i+1/2,j) = 0.5 * (model_(i,j) + model_(i+1,j))

du/dx |_(i+1/2,j) = (u_(i+1,j) - u_(i,j)) / dx

du/dy |_(i+1/2,j) =
    [(u_(i,j+1)   - u_(i,j-1))
   + (u_(i+1,j+1) - u_(i+1,j-1))] / (4 dy)

P_xface = DW(model_(i+1/2,j), I + grad_face(u))
```

The `y`-face formula is the coordinate-swapped analogue.  Equilibrium at a
node uses the difference of the traction columns on its low and high faces:

```text
R_i = [P_xface(:,x)_(i+1/2) - P_xface(:,x)_(i-1/2)] / dx
    + [P_yface(:,y)_(j+1/2) - P_yface(:,y)_(j-1/2)] / dy
    - b_i
```

Only column `d` of the stress on a face normal to coordinate direction `d`
enters that directional flux.

Near a material transition, constitutive evaluation, interpolation, and
differentiation do not commute:

```text
DW(model_i, grad_node(u)_i)
    != average[DW(model_face, grad_face(u))].
```

That mismatch was largest where `phi` blended propellant into a much stiffer
casing.  It produced the visible nodal overshoot even though the face
tractions driving the solve were smooth.

### New plotted stress

The new output reconstructs each nodal tensor column from the same two face
tractions used by the conservative operator:

```text
P_out,i(:,x) = 0.5 * [P_xface(:,x)_(i-1/2) + P_xface(:,x)_(i+1/2)]
P_out,i(:,y) = 0.5 * [P_yface(:,y)_(j-1/2) + P_yface(:,y)_(j+1/2)]
```

Thus the plotted `x` column represents the centered `x`-face traction and the
plotted `y` column represents the centered `y`-face traction.  It is a nodal
visualization of the solver's directional fluxes, not a fresh constitutive
evaluation at the node.

There are two important consequences:

1. The reconstructed columns are collocated at the node only by averaging;
   they originate on different face families.  The assembled tensor need not
   be exactly symmetric near a steep interface, even if the underlying
   continuum Cauchy stress is symmetric.
2. The output is intended for visualization and comparison with discrete
   equilibrium.  Applying a second nodal divergence operator to it is not a
   substitute for using the original face fluxes directly.

### Why it cannot feed back into the solve

The internal `stress_mf` is still computed with the legacy nodal constitutive
formula and remains available to mechanics internals.  A separate
zero-ghost, non-evolving `output_stress_mf` is registered under the plotted
`stress_*` names.  Only conservative, unmasked static solves fill that output
with the face reconstruction.  Dynamic and nonconservative/masked paths retain
their previous output semantics.

This separation was verified with paired runs: enabling versus disabling
stress plotting changed displacement, RHS, strain, material fields, `phi`,
`eta`, `psi`, and temperature by exactly zero.

### Symmetry-plane treatment

The face formula needs values on both sides of a node.  At the two genuine
quarter-domain reflection planes, the current code supplies those values by
parity:

```text
material:             even
normal displacement:  odd
tangential displacement: even
```

For this domain that means `u_x` is odd across `xhi`, while `u_y` is even;
`u_y` is odd across `ylo`, while `u_x` is even.  The identical face formula can
then be evaluated on both sides of the symmetry node.  `phi`, `eta`, and the
casing-support contour are not treated as boundaries.

On a physical boundary that is not configured as a reflection plane, the
code still uses the legacy nodal value because the centered tangential face
stencil has no defined exterior state.

## Hypotheses and tests

| Hypothesis | Why it seemed plausible | Test or change | Outcome |
|---|---|---|---|
| Insufficient base-grid resolution | Both `phi` and `eta` contours looked stair-stepped on the Cartesian AMR mesh | Compared 64/level-2, 96/level-2, and 128/level-1 at 1 s | Did not fix the line.  The 96 case had 1.5x finer finest spacing but a slightly larger nodal/face mismatch.  The two cases with identical finest spacing had nearly identical errors. |
| AMR hierarchy rather than finest spacing | 64/level-2 and 128/level-1 have different base grids and refinement layouts | Compared those two cases at the same 0.3426 mm finest spacing | Not the primary cause; interface mismatch was 0.211 versus 0.212 MPa.  Faint coarse/fine banding remains a secondary visualization/discretization effect. |
| A genuine stress concentration caused by the modulus jump | `phi` and `eta` both blend materials with large stiffness contrasts | Independently reconstructed the stress using the exact conservative face gradient and face-averaged model used by the solve | The recovered face stress was smooth through the thin visible line.  The alternating nodal overshoot/undershoot was therefore not a solver-resolved physical concentration.  A broader physical casing hoop-stress transition remains. |
| Cell-to-node or nodal stress reconstruction | The stored tensor was calculated at nodes while equilibrium is enforced through face tractions | Implemented face-consistent output stress and compared it with an independent face reconstruction | Substantially correct for the original thin line: maximum interface mismatch fell from 0.2108 to 0.0207 MPa.  It did not explain all subsequent boundary or axisymmetry behavior. |
| The broad bitmap transition was causing the feature | The original bitmap has an approximately 11.5 mm measured 1--99% `phi` transition | Replaced it in an isolation run with an analytic circle using tanh scale `pf.eps=0.125 mm` | Did not work at the existing mesh.  The approximately 0.5 mm transition was only 1.5 finest cells wide and increased the nodal/face mismatch to 0.942 MPa. |
| The phase-field epsilon should also be the material-interface width | `pf.eps` is the natural diffuse length in the burn model | Same epsilon-width isolation run | Rejected at the present resolution.  Burn-interface epsilon and elastic material-interface width cannot be equated unless the latter is adequately resolved. |
| The initial face-output change solved the entire problem | It removed the internal `phi` line in radial profiles | Recreated the case and inspected physical symmetry boundaries | Only partially successful.  Interior nodes used face stress while physical boundary nodes retained legacy nodal stress, creating a scheme discontinuity one node inside the boundary.  At the finest bottom boundary the phase-band jump grew from about 0.010 to 0.552 MPa. |
| A material/phase-aware boundary gradient was needed | The bad line followed `phi` and appeared close to a boundary | Considered a boundary-aware face gradient | Not implemented.  Treating `phi` or `eta` as geometric boundaries is inappropriate for this phase-field formulation. |
| Symmetry parity was missing at the quarter-domain boundaries | The right and bottom edges are actual reflection planes, unlike `phi` and `eta` | Extended the material evenly and displacement with odd-normal/even-tangential parity across `xhi` and `ylo` for output reconstruction | Worked for the boundary scheme mismatch.  Finest-level adjacent-node phase-band jumps are now about 8--9 kPa and symmetry-plane shear is only a few kPa. |
| The quarter-domain constraints were oriented incorrectly | An incorrect normal constraint would create tearing or rigid motion | Confirmed `xhi`: normal `u_x=0`, tangential traction free; `ylo`: normal `u_y=0`, tangential traction free | Constraints were already oriented correctly; not the cause. |
| The clamped square exterior was creating nonphysical casing stress | The analytical reference was circular, whereas the Alamo casing filled square corners and was fixed externally | Added a circular casing-support field, filled corners with soft void, and made the external rectangle traction-free | Improved the physical model and radial free-surface comparison, but did not remove the casing-axisymmetry residual. |
| The new circular casing/void support caused the axisymmetry peak | It introduced another diffuse modulus transition | Correlated polar shear with `casing_support` and radius | Ruled out as the main peak.  Maximum mean `|sigma_rtheta|` occurs near `r=0.083 m`, where `casing_support=1`; the outer support transition occurs later, near 0.086--0.088 m. |
| Face-consistent output itself causes the remaining axisymmetry error | Its tensor columns are reconstructed from directionally staggered faces | Reconstructed the legacy nodal constitutive stress from displacement and model data | Ruled out as the main remaining cause.  Mean casing `|sigma_rtheta|` was about 0.120 MPa for face output and 0.136 MPa for nodal constitutive stress.  The displacement also has a tangential component, so the asymmetry is in the solution. |
| `stress_validate` was reading the wrong fields after adding `casing_support` | The extra component shifted every subsequent hard-coded index | Changed the reader to resolve components by Header name | This was a validation-tool defect, not a simulation defect.  Fixing it prevents incorrect comparisons but does not alter Alamo stress. |
| Plotfile ghost nodes were contaminating radial bins | Rectangular plot seams initially resembled ghost-region overlays | Paired `Cell_H` valid boxes with FAB records and cropped ghosts | Robustness improvement only.  The current plotfiles have `nGrow=0`, so comparison CSV hashes and metrics were unchanged. |
| The analytical free-casing reference was wrong | The default JSON retained the old square's area-equivalent radius `a4=0.09896 m` | Added a free-circular config with `a4=0.0877 m` | The previous analytical curve was indeed misleading.  Correcting the radius reduced radial RMS error from about 0.436 to 0.195 MPa, but it did not change the simulated axisymmetry. |

## Attempts that did not solve the stress artifact

### Mesh and interface changes

- Increasing the base mesh to 96 x 96 while retaining `max_level=2`.
- Increasing the base mesh to 128 x 128 while reducing `max_level` to 1.
- Replacing the broad bitmap with an under-resolved `pf.eps`-scale analytic
  `phi` interface.
- Correcting AMR compositing and plotfile component selection in
  `stress_validate`.

These changed resolution or visualization details but did not remove the
underlying feature.  The epsilon-width interface made the local output error
worse because it was represented by only about 1.5 finest cells.

### Boundary and casing changes

- Retaining fixed exterior displacement while using the new circular support.
- Switching the external edges to zero traction.
- Replacing square casing corners with soft void.
- Updating the analytical model from fixed `u(a4)=0` to free
  `sigma_r(a4)=0`.

These were useful physical-model corrections, but none removed the remaining
angular casing variation.  That variation is already present inside the
fully supported region at the `phi` transition.

### Initial face-output boundary treatment

The first face-consistent output implementation deliberately fell back to
legacy nodal stress on every physical-domain boundary.  This removed most of
the internal line but created a new one-node discontinuity between boundary
and interior stress schemes.  The run was stopped after visual inspection.
Symmetry parity, applied only at the true reflection planes, replaced that
fallback and removed the large boundary jump.

### Solver configurations tried for the circular void case

The new `26 GPa -> 0.5 MPa` casing/exterior contrast destabilized the original
rediscretized multigrid hierarchy.  The following did not converge:

- restoring fixed outer displacement;
- disabling casing-support-specific refinement;
- changing only the bottom solver to BiCGStab while retaining coarsening;
- enabling `normalize_ddw` and `average_down_coeffs`;
- reducing smoother `omega` and increasing smoothing;
- disabling intra-level coarsening while retaining the bottom smoother.

The stable combination was `elastic.max_coarsening_level=0` with a BiCGStab
bottom solve.  This was a solver-enablement issue for the circular case, not a
fix for the stress-axisymmetry problem.

## What did work

1. Face-consistent stress output reduced the original internal `phi`-line
   reconstruction mismatch by 90.2% without changing displacement or any
   evolving field.
2. Symmetry-parity reconstruction eliminated the mixed boundary/interior
   stress scheme on `xhi` and `ylo`.
3. Circular casing support plus soft corners and free external traction
   produced the intended physical casing geometry.
4. The correct free-circular analytical configuration gives a good radial
   comparison: approximately 0.195 MPa RMS at 1 s, with radial stress tending
   to zero at the casing exterior.

## Evidence for the remaining cause

`base_circle0.bmp`, which defines `phi`, is not circular in physical space:

- horizontal `phi=0.5` radius: approximately 79.75 mm;
- vertical `phi=0.5` radius: approximately 81.25 mm;
- axis difference: approximately 1.5 mm, or 4.4 finest cells;
- approximate 10--90% transition width: 6.8 mm;
- image/transpose RMS difference: 0.0463, maximum 0.239.

For comparison, `blur5_rod_and_tube.bmp`, which initializes `eta`, has an
image/transpose RMS difference of only 0.00034.

Near the casing, mean `|sigma_rtheta|` is approximately 0.12 MPa against
28--30 MPa hoop stress, or about 0.4%.  The axisymmetry subplot visually
amplifies it by using a `10^5 Pa` scale.  It is nevertheless genuine: the mean
casing tangential displacement is about 8.4% of the mean radial displacement.

## Recommended decisive follow-up

Replace `phi.ic.bmp` with an exact radial expression centered at
`(0.0877,0.0877)` and rerun the one-second validation.  At the present finest
spacing of 0.3426 mm, use a resolved tanh scale around 0.35 mm initially.  A
literal 0.125 mm scale is under-resolved; approximately two additional AMR
levels would be needed to represent its 10--90% transition cleanly.

No production commit has been created.
