# LowMach Eulerian Solids Findings

This file summarizes the recent debugging work around the LowMach Eulerian
solids/reference-map implementation. It is intended as a handoff note for a
future agent or developer. The code was reverted to the most recent commit
before this file was written, so the experimental changes described below
should be treated as observations, not as currently active implementation.

## Current Problem

The current LowMach Eulerian-solids model evolves a phase field `eta`, a
reference map `xi`, computes deformation from `grad(xi)`, evaluates a finite
strain solid model, and couples the solid deviatoric Cauchy stress back into
the momentum equation.

The main observed failure mode is not simply elastic CFL stiffness. Reducing
the elastic wave speed makes the simulation smoother, but the same qualitative
problem remains:

- Cauchy stress accumulates in the diffuse/interfacial region, especially
  `0.5 < eta < 1`.
- The reference map also develops visible distortion in that same region.
- Once the interfacial distortion grows, stress singularities and velocity
  artifacts follow.

This suggests the diffuse boundary/reference-map treatment is the dominant
issue, not just the bulk elastic wave speed.

## Important Observations

1. Lowering solid wave speeds helps stability, but does not remove the
   interfacial stress accumulation.

2. Adding Newtonian viscosity everywhere, including in the solid, helped a lot.
   It damps the velocity/stress artifacts without changing the basic reference
   map model.

3. Higher-order advection, especially WENO5, seemed to help somewhat, but did
   not eliminate the diffuse-boundary artifacts.

4. The problematic stress growth often appears at the phase-field boundary and
   at corners of the solid body. Patch-boundary artifacts were also observed in
   parallel runs, so boundary fills and ghost-cell validity should remain on the
   checklist.

5. A prior diagnostic showed that `solid_first_piola_kirchhoff_stress_xy`
   appeared to have a sharp cutoff at `eta = 0.5`. The cause was that the raw
   solid stress was only being evaluated where the solid weight was positive.
   This is important diagnostically: output fields can show an artificial
   cutoff even if the intended physical coupling is weighted smoothly.

6. Treating a higher threshold, such as `eta > 0.7` or `eta > 0.9`, as the
   "truth" region for `xi` was considered. This would smooth/reconstruct more
   of the boundary region, but it risks erasing real deformation history near
   the material surface and making the solid artificially soft/slippery.

7. An attempted change using `reference_map.eta_core = 0.7` made the result
   worse. The likely reason is that part of the mechanically important
   interior shell stopped carrying true advected deformation history and was
   instead reconstructed.

8. Another attempted change remapped the mechanical solid weight so that
   `eta >= 0.5` was fully solid and the blurred transition lived only in
   `eta < 0.5`. Conceptually this matches a level-set blurred-boundary view,
   but the concrete implementation broke the current behavior and was reverted.
   Do not blindly reapply it.

9. Explicitly snapping/resetting `xi` or `alpha` type fields tends to introduce
   artificial gradients. Earlier work suggested that if smoothing is needed, it
   should be applied to stress evaluation or extrapolation, not by hard
   overwriting the transported field in a way that creates discontinuities.

10. Advecting `eta` as a passive phase-field marker without maintaining its
    structure may be part of the problem. However, a previous phase-field
    evolution attempt diffused the boundary too much, so any Allen-Cahn or
    Cahn-Hilliard regularization needs to be conservative and carefully tuned.

## Experiments That Were Tried And Reverted

These should be regarded as debugging attempts, not known-good fixes.

- Stress-only smoothing of the reference map:
  smoothing `q = xi - x` before stress evaluation, while not overwriting the
  transported `xi`. This is still conceptually plausible, but the exact
  implementation should be revisited carefully.

- Using `reference_map.eta_core > eta_cutoff`:
  making only a deeper solid core protect the true map, while reconstructing
  the shell. This made behavior worse in the tested case.

- Remapping `solid_weight` so `eta >= 0.5` was fully solid and only
  `eta < 0.5` was the blurred layer:
  conceptually plausible, but the implementation tested here broke the run.
  If revisited, it should be done with much more diagnostic output and probably
  separated from changes to reference-map reconstruction.

## Recommended Next Steps

1. Add diagnostics before changing the model again. Useful fields:
   - `det_grad_xi`
   - `J = det(F)`
   - `|grad_xi - I|`
   - raw elastic Cauchy stress
   - weighted/coupled elastic deviatoric stress
   - `solid_weight`
   - `eta`
   - possibly `grad(eta)` magnitude

2. Separate raw stress diagnostics from coupled stress diagnostics. It should
   be clear whether a large stress exists only in the raw reconstructed field
   or is actually entering `div(sigma)`.

3. Verify exactly where `xi` is advected, reconstructed, smoothed, or reset.
   The main suspect is the interfacial band where the map is neither cleanly
   true material history nor cleanly extrapolated boundary data.

4. If smoothing is added, prefer a stress-evaluation smoothing path:
   make a temporary `xi_for_stress`, smooth/extrapolate it, compute stress from
   that, and leave the transported `xi` untouched. Avoid hard snapping the live
   `xi` field.

5. If the blurred-boundary idea is revisited, change only one thing at a time:
   first change the stress weight, then inspect raw and weighted stresses, then
   separately modify reference-map extrapolation if still needed.

6. Keep the simple `input.lm.couette_solid` style case as the primary debug
   case. It is much cheaper and easier to interpret than the driven cavity with
   a moving solid.

7. Before reintroducing any aggressive model change, verify the no-solid or
   eta-marker-only LowMach flow still behaves correctly with AMR and parallel
   runs. This solver has previously shown sensitivity at coarse/fine and domain
   boundaries.

## Working Assumptions

- The baseline branch after revert is the best starting point.
- Use the 2D clang build unless explicitly changed.
- Do not reintroduce custom average-down, prolongation, or boundary-fill logic
  unless there is a demonstrated gap in the existing Alamo/AMReX machinery.
- The ultimate goal is still a phase-field/blurry-boundary Eulerian solid
  method, not a level-set rewrite.

