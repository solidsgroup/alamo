# LowMach Eulerian Solids Findings

Last updated: 2026-07-14

This file summarizes the current LowMach Eulerian-solid investigation. Unlike
the previous version of this note, the reference-map, AMR, and phase-field
changes described as active below are present in the current worktree. The main
implementation is in `src/Integrator/LowMach.cpp` and `LowMach.H`; the tuned
debug case is `input.lm.couette_solid`.

## Current Status

The solver evolves the solid phase field `eta` and reference map `xi`, computes
`F = inverse(grad(xi))`, evaluates the finite-strain solid model, and couples the
weighted deviatoric Cauchy stress into momentum.

Three related failure modes have been investigated:

1. The transported reference map degraded in the diffuse boundary. This first
   appeared as boundary stress concentration, then visible `xi` distortion,
   and eventually instability.
2. Enabling AMR changed the solution because the interface crossed coarse/fine
   transitions, the refinement mask lagged the rotating body, and projection
   needed to be composite across levels.
3. Conservative phase-profile compression generated directional lobes at the
   original upper-right and lower-left corners. Suppressing that compression
   protected the `eta=0.5` contour but left low-eta trails.

The current implementation addresses all three mechanisms. In the tested AMR
case through `t=20`, the long corner lobes are gone, the effective phase width
remains near its target, and the `eta>0.5` area changes by about one percent.
Longer runs are still needed before treating this as a final model.

## Reference-Map Reconstruction

### Failure mechanism

The old reconstruction divided the domain into discrete known and unknown
regions using eta thresholds. A cell could therefore switch abruptly from
transported material history to reconstructed data as the diffuse boundary
moved. The resulting mismatch in `grad(xi)` was mechanically important even
when the map looked visually reasonable.

Raising `reference_map.eta_core` alone was not a fix. A test with a deeper
truth region made the result worse because too much of the mechanically active
shell stopped carrying advected deformation history.

### Active graded algorithm

Reconstruction and smoothing are now continuous functions of eta. For a core
value `eta_core`, the reconstruction grade is

```text
g(eta) = 1 - SmootherStep(clamp(eta / eta_core, 0, 1)).
```

The reconstruction and smoothing strengths are

```text
repair = reconstruction_alpha * g^reconstruction_power
relax  = smoothing_alpha      * g^smoothing_power.
```

This keeps the interior map unchanged, applies progressively stronger repair
toward the exterior, and removes artificial switching surfaces.

Each reconstruction sweep:

- Uses only coordinate neighbors with larger eta, so information propagates
  outward from the material interior.
- Weights candidates by the eta increase toward the neighbor.
- Uses the nearest inward value as the robust fallback.
- Extends the inward slope when deeper samples are available.
- Preserves an affine map across low-eta cells when consecutive inward slopes
  have sufficiently low relative curvature.
- Blends the reconstruction into the live map using the eta grade rather than
  replacing it at a threshold.

Smoothing is also eta graded. It is strongest in the exterior extension and
goes continuously to zero at the protected core. Boundary conditions and
periodic fills are refreshed between every sweep.

The current Couette-solid settings are:

```ini
reference_map.eta_core = 0.7
reference_map.eta_extension = 1.0e-3
reference_map.extrapolation_sweeps = 4
reference_map.reconstruction_alpha = 1.0
reference_map.reconstruction_power = 1.0
reference_map.affine_tolerance = 1.0e-7
reference_map.smoothing_sweeps = 4
reference_map.smoothing_alpha = 1.0
reference_map.smoothing_power = 2.0
```

This construction is intended to support phase growth as well as advection:
newly occupied low-eta cells receive an outward continuation of the existing
reference map instead of an identity reset or an abrupt copied value.

## AMR Findings

AMR changed more than resolution in this problem:

- Coarse/fine interpolation and average-down alter transported `eta` and `xi`.
- Regridding can place the diffuse layer directly on a patch transition.
- Subcycling changes the number and timing of reconstruction and phase-source
  updates by level.
- A level-local projection does not enforce one composite divergence constraint
  across coarse/fine interfaces.

`projection.amr_enabled = 1` is now used so the correction solve is composite
across all active levels before velocity is averaged down.

The original AMR phase trigger was also too narrow. With
`eta_refinement_criterion = 0.1`, the `eta=0.1` contour sat at the edge of the
fine patch. With `amr.regrid_int = 1000` and `dt = 0.004`, the grid was rebuilt
only every four time units, allowing the rotating corners to outrun the fine
halo. This produced stair-stepped outer contours and contributed to trailing
material.

The active AMR settings are:

```ini
amr.n_cell = 64 32
amr.max_level = 2
amr.nsubsteps = 2
amr.regrid_int = 50
eta_refinement_criterion = 0.01
projection.amr_enabled = 1
```

Lowering the eta criterion made the `eta=0.01` contour substantially smoother.
Frequent regridding gave a smaller additional improvement by keeping the wider
halo centered on the moving body. These changes do not, by themselves, remove
phase trails; they remove AMR-induced stair-stepping and mesh lag.

## Phase-Field Stabilization

### Initial-condition correction

The old box expression used

```text
max(abs(x-x0)-hx, abs(y-y0)-hy)
```

as though it were a signed distance. It is an L-infinity distance, so all
diffuse contours inherit sharp corner-normal jumps. The active expression is
the Euclidean signed distance to a rectangle. The `eta=0.5` geometry remains a
sharp rectangle, while exterior diffuse contours are correctly rounded.

### Mapped-distance conservative flux

For the target profile

```text
eta = 0.5 * (1 + tanh(d / epsilon)),
```

the code recovers the distance coordinate

```text
d = 0.5 * epsilon * log(eta / (1 - eta)).
```

Using `d` makes the conservative diffusion/compression balance exact for a
resolved planar tanh profile. With `w = eta * (1-eta)`, the ungraded mapped flux
has the form

```text
q = w * grad(d) * (1 - counter_curvature / |grad(d)|).
```

This was more reliable than applying the old eta-gradient flux directly, but
full compression still distorted sharp corners because a corner has no unique
normal.

### Normal-coherence grading

The code now measures agreement among mapped-distance normals in a neighborhood
of each face. If `c` is the magnitude of their average, the compression grade is

```text
C = clamp((c - coherence_threshold) / (1 - coherence_threshold), 0, 1)
    ^ coherence_power.
```

The face flux is

```text
q = w * grad(d) * (D - C * counter_curvature / |grad(d)|)
D = C + (1-C) * incoherent_diffusion.
```

Thus smooth interface sections retain the full balanced flux. Compression is
removed continuously where neighboring normals disagree, while a bounded
fraction of diffusion remains to prevent unresolved corner noise.

### Exterior signed-distance relaxation

Coherence grading fixes the material contour but, by itself, leaves low-eta
material behind the two extensional corners. Restoring the conservative
compression there only sharpened those trails into thin filaments.

The active algorithm therefore adds a local signed-distance relaxation only for
`eta < 0.5` and only where normal coherence is poor. Its source is proportional
to

```text
mobility * exterior_reinitialization * (1-C)
* (2*w/epsilon) * smooth_sign(d) * (1-|grad(d)|).
```

This drives the exterior profile toward `|grad(d)| = 1`. The source is zero at
`eta=0.5`, so it does not directly move the material contour. It is deliberately
not globally conservative: it removes low-eta trail mass. Restricting it to the
exterior was important. A two-sided version expanded the enclosed area.

The active phase settings are:

```ini
eta.phase_field.enabled = 1
eta.phase_field.epsilon = 0.015
eta.phase_field.mobility = 0.01_m/s
eta.phase_field.counter_curvature = 1.0
eta.phase_field.band = 1.0e-3
eta.phase_field.mapped_distance = 1
eta.phase_field.normal_coherence_threshold = 0.9
eta.phase_field.normal_coherence_power = 2.0
eta.phase_field.exterior_reinitialization = 1.0
eta.phase_field.incoherent_diffusion = 0.25
```

## Parameter Screens And Rejected Variants

The following conclusions came from matched runs of the Couette-solid case:

- Standard conservative compression at mobility `0.05_m/s` recreated a corner
  lobe by about `t=4`. Mobility `0.02_m/s` was also less clean than `0.01_m/s`.
- Disabling all corner flux protected `eta=0.5` but left pointed low-eta tails.
- Retaining all corner diffusion removed some tails but rounded the material
  corners and broadened the interface.
- `incoherent_diffusion = 0.25` was the best tested compromise. A value of
  `0.5` broadened the effective epsilon more and developed outer oscillations.
- Restoring conservative compression only in the exterior narrowed trails into
  filaments rather than removing them. That experiment is not in the final code.
- Two-sided local signed-distance relaxation reduced width error but increased
  the `eta>0.5` area by about 1.7 to 2.2 percent at `t=10`.
- Exterior-only relaxation avoided that expansion. Strength `0.5` retained
  visible lobes at `t=20` and lost more enclosed area than strength `1.0`, so
  the active strength is `1.0`.
- Changing only the AMR eta threshold from `0.1` to `0.01` improved contour
  smoothness and reduced AMR phase-mass drift slightly. Regridding every 50
  steps gave a further smaller improvement.

## Quantitative Results

Metrics were evaluated on a level-2 covering grid. Effective epsilon was
computed as

```text
epsilon_eff = 2 * integral(eta*(1-eta)) / integral(|grad(eta)|).
```

The nominal epsilon is `0.015`. The reported width is
`area(0.1 < eta < 0.9) / integral(|grad(eta)|)`.

| Case | Time | Phase mass change | `eta>0.5` area change | Effective epsilon | 0.1-0.9 width |
| --- | ---: | ---: | ---: | ---: | ---: |
| Final exterior relaxation | 10 | -0.494% | +0.403% | 0.015389 | 0.033702 |
| Final exterior relaxation | 20 | -2.307% | -1.048% | 0.015326 | 0.033531 |
| Same AMR, exterior relaxation off | 10 | +0.513% | +0.645% | 0.015696 | 0.034102 |

At `t=20`, the final run retained a compact `eta=0.1` contour and had no loop
or lobe in the `eta=0.5` contour. The phase-mass loss is larger than the material
area loss because the exterior relaxation intentionally removes diffuse trails.

For comparison, a uniform-grid conservative run without exterior relaxation
changed phase mass by approximately `1.3e-4` through `t=8`. This supports the
conclusion that most earlier positive mass drift came from AMR transport and
regridding rather than the conservative face flux.

## Integrator Streamlining

The LowMach hot path previously treated derived diagnostics as integration
state. Every RK RHS and post-stage callback recomputed and ghost-filled
momentum, total energy, vorticity, solid weight, deformation gradient,
first-Piola stress, total Cauchy stress, and two nodal stress tensors. Most of
those fields were never consumed by the evolution equations.

The streamlined implementation now:

- computes only density, mole fraction, and weighted cell solid stress in the
  RHS;
- keeps the required post-stage and post-projection reference-map passes;
- removes the registered projection RHS scratch field and allocates it only in
  the single-level projection where it is used;
- removes the two unused nodal stress fields, including the unnecessary nodal
  restart field;
- updates plot diagnostics through an integrator plot-preparation hook, so they
  are current at initial, periodic, and final plot writes without being updated
  on every timestep;
- registers momentum, energy, vorticity, solid weight, `F`, first-Piola stress,
  total Cauchy stress, and mole-fraction output only when
  `diagnostics.extended_fields = 1`.

The default registry is reduced from 23 cell fields plus 2 nodal fields to 15
cell fields and no nodal fields. The compact plotfile retains velocity,
temperature, mass fraction, `eta`, `xi`, density, pressure, pressure correction,
and weighted solid deviatoric Cauchy stress. Set
`diagnostics.extended_fields = 1` to recover the larger cell-diagnostic output
set for focused analysis.

On the 12-rank, 60-step AMR Couette-solid benchmark, wall time decreased from
`6.14 s` to `5.47 s` without plot output. With the old and extended diagnostic
plot sets enabled at step 60, wall time decreased from `6.28 s` to `5.58 s`.
These short-run measurements indicate an approximately 11 percent speedup;
startup and final I/O are included. A separate six-rank Couette run was active
on the host during both measurements, so this is an indicative paired result,
not a controlled performance benchmark.

An AMReX plotfile comparison against the preserved pre-streamlining executable
at step 60 found identical temperature, composition, and density. Maximum
differences were approximately `3e-8` in `xi`, `1e-6` in `eta`, and `9e-8` in
velocity, consistent with floating-point propagation through the parallel
projection and stress kernels. A trial that removed RK post-stage or the second
post-projection reference-map pass produced materially larger differences and
was rejected.

## Validation Performed

- Clean 2D clang build passed.
- 3D clang compilation passed; no full 3D solid run has been performed.
- New mapped-distance and legacy `mapped_distance=0` paths passed one-step MPI
  smoke tests.
- The selected AMR configuration was run through `t=20`.
- The installed `bin/lowmach-2d-clang++` matches the selected tested build.
- `git diff --check` passes.

The selected `t=20` plotfile is under:

```text
/tmp/alamo-phase-outer-local-reinit-100-amr-t20-cont-20260714/05000cell
```

## Remaining Risks And Next Steps

1. Run the selected configuration through at least `t=50` and preferably the
   configured `stop_time = 100`. Track phase mass, `eta>0.5` area, effective
   epsilon, and corner positions to determine whether the `t=20` area loss
   saturates or continues.
2. Re-evaluate reference-map and stress diagnostics with the stabilized phase
   boundary. The phase fix removes a major source of boundary motion, but it
   does not prove that long-time `xi` degradation is eliminated.
3. Add plot fields for `det(grad(xi))`, `J = det(F)`, `|grad(xi)-I|`, raw solid
   stress, weighted solid stress, normal coherence, and exterior phase source.
4. Quantify AMR cost from `amr.regrid_int = 50`. A larger interval may be safe
   if `amr.n_error_buf` is increased enough to keep the full diffuse layer on
   the finest level.
5. Run a true 3D corner test. The 3D path compiles, but its 27-cell coherence
   neighborhood has not been validated dynamically.
6. Keep `eta.phase_field.exterior_reinitialization` explicit in production
   inputs. It is non-conservative by design and should not be enabled silently.

## Historical Observations Still Relevant

- Lower solid wave speeds improve explicit stability but do not remove the
  diffuse-boundary failure mechanism.
- Newtonian viscosity, including in the solid, damps velocity and stress
  artifacts but is not a reference-map repair.
- WENO5 advection helped somewhat but did not eliminate boundary artifacts.
- Raw and weighted stress outputs must remain separate. A previous raw
  first-Piola field showed an artificial cutoff because stress was evaluated
  only where the solid weight was active.
- Hard resetting `xi` or displacement-like fields creates artificial gradients.
  Reconstruction should remain graded and should propagate material history
  outward rather than snapping to identity.
- Remapping the mechanical weight so `eta>=0.5` was fully solid was tested and
  broke the current behavior. Do not reapply it without isolated diagnostics.

## Working Assumptions

- `input.lm.couette_solid` remains the primary debug case.
- Use the 2D clang build unless explicitly testing dimensional behavior.
- Preserve the current graded reference-map algorithm while evaluating the
  phase fix; changing both again would make regressions difficult to attribute.
- Do not add custom AMR prolongation, average-down, or boundary-fill logic
  without first demonstrating a gap in AMReX's existing machinery.
- The target remains a diffuse-boundary Eulerian-solid method, not a sharp
  level-set rewrite.
