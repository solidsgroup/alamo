# LowMach Eulerian Solids Findings

Last updated: 2026-07-21

This file is a technical handoff for the current `eulerian-solids` LowMach
implementation. It describes code active in the present branch, results
verified through 2026-07-20, and remaining modeling gaps. Older mapped-distance
reinitialization, nodal-velocity, custom regridding, and implicit-elastic
experiments are historical only and are identified as such below.

The principal files are:

- `src/Integrator/LowMach.H` and `src/Integrator/LowMach.cpp`
- `src/Operator/PressurePoisson.H` and `PressurePoisson.cpp`
- `src/Numeric/ReferenceMap/Reconstruction.H` and `Reconstruction.cpp`
- `src/Numeric/Advect/`
- `src/Model/PhaseField/AllenCahn.H`
- `src/Model/Mechanism/PhaseChange.H`
- `src/Model/Chemistry/`

The most useful inputs are:

- `tests/LMDrivenCavity/input`: pure-fluid AMR regression
- `tests/LowMachChemistry/input`: finite-rate chemistry regression
- `input.lm.couette_solid`: compact deformable-solid mechanics case
- `input.lm`: larger deformable-solid driven-cavity case
- `input.lm.rigid`: fixed rigid inclusion
- `input.lm.rigid_effusion`: rigid-to-gas phase change in crossflow
- `input.lm.ap_htpb`: AP/HTPB regression and gas combustion

## Executive Status

LowMach now has one species representation for fluids and solids. Conserved
partial densities are authoritative; `eta`, total density, mass fractions, and
mole fractions are derived from them. Each named species has one mechanics
classification:

```text
fluid | deformable_solid | rigid_solid
```

The current solver supports:

- any number of gas species described by one gas model;
- one deformable-solid species with one reference map and a Neo-Hookean model;
- multiple rigid-solid species sharing one prescribed rigid velocity;
- mass-conservative condensed-to-gas phase-change mechanisms;
- frozen, finite-rate, or six-species Rocfire gas chemistry;
- explicit or locally implicit chemistry;
- a hierarchy-wide variable-coefficient pressure projection;
- selectable collocated advection operators; and
- AMR refinement based on velocity, pressure, temperature, or reconstructed
  solid volume fraction.

The pure-fluid and chemistry regressions pass in serial/MPI AMR configurations.
The AP/HTPB MPI/AMR case reaches 100 microseconds without the former diffuse-
interface temperature singularity or timestep collapse. Deformable-solid
mechanics remain experimental: the reference-map treatment is substantially
better than the original thresholded reconstruction, but long-time stability
with explicit elastic feedback is not yet established.

## Authoritative State And Species

`component_density_mf` contains one Eulerian partial density for every named
species. Species names are user-supplied identifiers, and each species has its
own initial condition, for example:

```ini
species.names = AP_gas HTPB_gas AP_solid HTPB_solid

AP_gas.mechanics = fluid
HTPB_gas.mechanics = fluid
AP_solid.mechanics = rigid_solid
HTPB_solid.mechanics = rigid_solid

AP_solid.density.ic.type = expression
AP_solid.density.ic.expression.region0 = "..."
```

The first `gas.nspecies` entries must be fluid species and must match the gas
property arrays. Later entries are condensed species. `UpdateComponentState`
reconstructs

```text
rho             = sum_n rho_n
eta_deformable  = rho_deformable / rho_ref_deformable
eta_rigid       = sum_{n in rigid} rho_n / rho_ref_n.
```

Mass and mole fractions are derived only for gas species. They are diagnostics,
not evolved state. `eta_mf` and `rigid_eta_mf` are likewise derived fields; no
independent eta transport equation remains.

This design is important for phase change. A mechanism transfers equal mass
between named partial densities. The volume fraction changes because the
source and product have different equations of state or reference densities,
not because an independently evolved eta is adjusted afterward.

Current restrictions:

- only one deformable-solid species is allowed;
- all rigid species use the same target velocity and relaxation time;
- condensed species require a positive reference density; and
- one velocity and one temperature field are shared by the mixture.

## LowMach Evolution

All primary evolved fields are cell centered. Velocity was experimentally
migrated to nodes, but that path created widespread stencil, boundary, AMR, and
projection inconsistencies and was reverted. The present method is collocated,
not staggered.

The explicit momentum predictor contains

```text
-u dot grad(u) + g + (mu/rho) laplacian(u)
+ sign/rho * div(sigma_solid_dev).
```

Pressure is not added to this predictor. It is applied once by the
nonincremental pressure projection. This avoids feeding the stored pressure
correction back into the following predictor.

The Newtonian viscosity from the gas transport model is applied throughout the
mixture, including the solid region. The deformable-solid model can add a
second solid viscosity and an interface-localized deviatoric damping stress.

Partial densities use conservative advection,

```text
partial_rho_dot = -div(partial_rho * u) + mechanisms + diffusion + chemistry,
```

while velocity, temperature, and the reference map use the advective form.
SSPRK3 is used by the current principal inputs. Chemistry may additionally use
Strang splitting around the transport update.

## Advection Interface

`Numeric::Advect::Advect` is a GPU-safe pseudo-polymorphic functor. LowMach
selects an operator in `Parse` and invokes the same object for scalars and
vectors. Each implementation exposes its required ghost width and the expected
locations of phi and velocity. LowMach currently rejects anything other than
collocated cell-centered data.

Available implementations are:

- `upwind`
- `centered`
- `quick`
- `muscl`, with MC, minmod, superbee, van Leer, van Albada, Koren, and UMIST
  limiters
- `weno5`

MUSCL and WENO5 reconstruct face states but obtain face velocities by averaging
neighboring cell velocities. This is not exactly the same face flux retained by
the pressure projection. That mismatch is one source of mixture-volume drift
and is discussed below.

## Deformable-Solid Mechanics

The deformable species carries a cell-centered inverse reference map `xi`.
Inside mechanically trusted material,

```text
grad_xi = grad(xi)
F       = inverse(grad_xi)
J       = det(F)
sigma   = P(F) F^T / J.
```

The current constitutive model is `Model::Solid::Finite::NeoHookean`. Only the
deviatoric part of its Cauchy stress is returned to the explicit momentum RHS.
Volumetric elastic stress is not also applied explicitly; the projection
enforces the mixture incompressibility constraint. This avoids double counting
a pressure-like volumetric contribution.

The mechanical weight is smooth above `solid.model.eta_threshold`:

```text
w_s(eta) = SmootherStep((eta - eta_threshold)/(1 - eta_threshold))
           for eta > eta_threshold,
w_s(eta) = 0 otherwise.
```

Stress is evaluated wherever `eta > reference_map.eta_extension`, then the
deviatoric stress is multiplied by `w_s`. Interface damping is

```text
mu_interface * 4 eta (1-eta) * dev(grad(u) + grad(u)^T).
```

The retained tensor state is a typed `Set::Field<Set::Matrix>` named
`solid_deviatoric_stress`. The divergence is computed with the matrix
divergence stencil. Extended diagnostics retain `F`; the previous collection of
raw Piola, duplicate Cauchy, weighted stress, and nodal stress fields was
removed.

Elastic feedback remains explicit. The dynamic timestep includes an elastic
wave estimate and viscous restrictions, but there is no implicit elastic
operator. The former `Operator/ElasticLowMach` and
`ImplicitElasticVelocitySolve` were deliberately removed because their scalar
linearization was not a general implicit discretization of the selected finite
strain model.

## Reference-Map Reconstruction

Directly advecting `xi` is not sufficient near a diffuse boundary. Cells enter
and leave the mechanically active region, and AMR coarsening can replace fine
material history with interpolated coarse data. Hard resets, identity snaps,
and abrupt eta cutoffs produced large artificial `grad(xi)` and stress jumps.

The active reconstruction is isolated in
`Numeric::ReferenceMap::Reconstruction`. For
`reference_map.eta_core = eta_core`, define

```text
g(eta) = 1 - SmootherStep(clamp(eta/eta_core, 0, 1)).
```

Reconstruction and smoothing strengths are

```text
repair = reconstruction_alpha * g^reconstruction_power
relax  = smoothing_alpha      * g^smoothing_power.
```

The interior approaches zero repair continuously. Exterior cells extrapolate
from coordinate neighbors with larger eta. The nearest inward continuation is
the fallback; deeper inward samples permit a linear continuation when the map
is locally affine. Several sweeps propagate this information through the
extension region. Smoothing uses the same eta grading and is confined to cells
with nearby material support.

An optional second smoothing pass operates on a temporary copy used only for
stress calculation. It smooths `xi-x`, rather than absolute coordinates, so an
identity reference map remains identity. The transported map itself is not
modified by this stress-only pass.

The commonly used current settings are:

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
reference_map.stress_smoothing_sweeps = 2
reference_map.stress_smoothing_alpha = 0.05
```

Reconstruction runs for RK stage states and again after the completed update.
The repeated boundary fills inside this numerical processor are intentional;
the old integrator-wide custom state fill/regrid machinery is not.

## Rigid Solids

Rigid solids use an implicit local Brinkman relaxation embedded in the
projection rather than an explicit penalty-force field. For prescribed rigid
velocity `U_r`, rigid volume fraction `eta_r`, and relaxation time `tau`,

```text
m = 1 / (1 + dt * eta_r/tau)
u_star = m*u + (1-m)*U_r.
```

The pressure coefficient is correspondingly

```text
beta = m/rho.
```

Thus the projection accounts for the reduced velocity mobility inside the
rigid phase. The former `rigid_penalty_force_mf` was removed. This is suitable
for fixed or prescribed-motion diffuse solids, but it is not a six-degree-of-
freedom rigid-body solver and does not calculate hydrodynamic force or torque.

## Phase Field And Phase Change

The old standalone eta reinitialization and mapped-distance counter-curvature
implementation were removed. The only active phase-field model is the standard
Allen-Cahn model in `Model/PhaseField/AllenCahn.H`, parameterized either by
`lambda,kappa` or by `sigma,epsilon`.

Allen-Cahn is used through a named `Model::Mechanism::PhaseChange`. A mechanism
identifies one condensed input species and one or more gas output species. It
computes an eta rate from the input partial density, converts it to a mass
source, subtracts that mass from the condensed species, and distributes the
same mass among the products. Product mass fractions must sum to one.

The current phase-change implementation deliberately applies

```text
eta_dot = min(eta_dot_AllenCahn, 0).
```

It therefore supports regression, decomposition, sublimation, or pyrolysis,
but not condensation or phase growth. With the default `interface_only=1`, the
rate is additionally weighted by `4 eta (1-eta)`. Optional temperature cutoff,
Arrhenius activation, and a rate multiplier are available.

Mass transfer is conservative even though Allen-Cahn is not conservative in
eta. The mechanism also supplies the projection with the volume source

```text
mass_source * (1/rho_ref_condensed - R_product*T/p0).
```

The deformable-solid demonstration inputs use `interface_only=0` and zero
driving force mainly to maintain a regular profile while transferring any lost
solid mass to the gas species. The AP/HTPB input uses large artificial rate
multipliers so regression is observable over a short test. Those multipliers
are testing parameters, not calibrated propellant kinetics.

## Pressure Projection And AMR

`LowMach::ProjectVelocity` now contains the physical construction of the
projection source and mobility. AMReX setup and scratch storage are encapsulated
in `Operator::PressurePoisson`.

The operator:

1. follows the current hierarchy layout;
2. owns cell RHS, coefficient, solution, divergence, and face scratch fields;
3. fills coarse/fine coefficient ghosts;
4. averages the cell coefficient to faces;
5. constructs one `MLABecLaplacian` over all active levels;
6. computes the RHS from the divergence of averaged face velocity;
7. solves with MLMG; and
8. forms the correction from the same face coefficient and face pressure
   gradient before averaging the correction back to cell velocity.

Using the same face discretization for divergence, coefficient, and pressure
correction removed an early one-cell interface velocity spike. Pure-fluid
projection behavior is protected by `tests/LMDrivenCavity`.

LowMach no longer owns custom state prolongation, average-down, regrid, or
`FillStateBoundaries` logic. Primary state communication is left to the
Integrator/AMReX framework. `PressurePoisson` still uses a local
`FillPatchTwoLevels` for its private coefficient scratch field; that is solver
scratch preparation, not an alternate state-management path.

## Chemistry Port And LowMach Coupling

The finite-rate and Rocfire kinetics were ported from
`origin/flame-with-multicomponent` with limited model changes:

- compile-time `NSPECIES` became runtime `ngas_species`, backed by a GPU-safe
  `MAX_SPECIES=32` array;
- the pseudo-polymorphic wrapper gained `Reactive`, `Implicit`, iteration, and
  tolerance queries;
- parsing now receives the active gas-species count;
- finite-rate heat release uses species enthalpy instead of internal energy,
  because LowMach advances a constant-pressure temperature/enthalpy equation;
- `Equilibrium` chemistry was not ported; and
- the Cantera YAML parser itself was retained unchanged.

LowMach stores Eulerian partial densities, while the chemistry models expect
intrinsic gas density. At fixed thermodynamic pressure `p0`, it reconstructs

```text
alpha_g = rho_g * R(Y) * T / p0
rhoY_intrinsic = rhoY_eulerian / alpha_g.
```

Chemistry is evaluated on the intrinsic state, then species and heat sources
are multiplied by `alpha_g` when returned to mixture-volume equations.

For implicit chemistry, LowMach uses Strang splitting. Each half-step solves
gas mass fractions and temperature together with a local backward-Euler/Newton
solve at fixed `p0`. A cell first attempts the full half-step and halves its
local chemistry step only after a failed positive Newton solve. This removes
the chemical-kinetics restriction from the global flow timestep; transport,
phase field, viscosity, diffusion, or advection can still limit it.

## Diffuse-Interface Chemistry Failure And Fix

The first AP/HTPB implicit runs developed a hot spot in partially solid cells.
Temperature rose above 13,000 K near the AP/HTPB/eta interface, species
diffusivity increased rapidly, and the explicit transport restriction drove the
global timestep toward `1e-10 s`. The instability was not an acoustic CFL or a
failure of the local chemistry Newton solve.

The problem was thermodynamic weighting. A cell containing a small gas volume
and a large condensed partial density received the intrinsic gas adiabatic
temperature rise as though only the small gas mass contributed heat capacity.
Chemistry dilatation was also initially treated as a full-cell source.

The active coupling now uses

```text
T_dot_reaction = alpha_g * qdot_intrinsic / (rho_mixture * cp_gas)
```

and the implicit residual uses the equivalent expression

```text
dt * rho_g * qdot_intrinsic
-----------------------------------------
rho_g_intrinsic * rho_mixture * cp_gas
```

The integrated implicit dilatation is

```text
rho_g * (1/rho_g_intrinsic,new - 1/rho_g_intrinsic,old)
```

and explicit thermal and molar dilatation terms are multiplied by `alpha_g`.
This prevents a tiny gas pocket in the diffuse solid boundary from receiving a
full-cell reaction heat or expansion source.

Conservative scalar advection and the collocated projection still accumulate a
small split volume error. In mixed-phase cases only, the projection therefore
adds the standard one-step discrepancy correction

```text
V = alpha_g + eta_deformable + eta_rigid
S_discrepancy = (V-1) / (dt * max(V,0.1)).
```

Scoping this correction to mixed-phase cases was necessary: applying it to a
pure gas slightly changed the small transverse velocity in the driven-cavity
regression without solving a relevant problem there.

The present heat-capacity treatment is still approximate. Condensed material
contributes `rho_condensed * cp_gas` because no condensed-species heat-capacity
model exists yet. The volume and source weighting are structurally correct,
but quantitatively accurate combustion/regression will require per-species
condensed thermodynamics and phase-change enthalpy.

## Verified Results

The following results were checked on 2026-07-20 with the 2D clang build.

| Case | Configuration | Result |
| --- | --- | --- |
| `LMDrivenCavity` | Re=100, AMR levels 0-1, serial | Run and reference-profile check passed |
| `LMDrivenCavity` | Re=100, AMR levels 0-1, MPI 2 | Run and reference-profile check passed |
| `LowMachChemistry` | 29-reaction H2/O2, explicit, MPI 2, AMR | Product, positivity, pressure, and AMR checks passed |
| `LowMachChemistry` | Same case, implicit with 10x flow timestep | Product, positivity, pressure, and AMR checks passed |
| `input.lm.ap_htpb` | MPI 2, AMR levels 0-3 | Reached 100.041 microseconds in 895 steps |

Final AP/HTPB metrics were:

```text
flow timestep            84.76 ns
maximum temperature      1687.95 K
maximum speed            8.88 m/s
mixture volume range     0.942 to 1.050
eta_rigid at Tmax        5.5e-17
minimum partial density  0.0
```

The maximum temperature moved out of the diffuse solid boundary and was in
pure gas by 20 microseconds. At 20.028 microseconds, serial and two-rank fields
agreed to at most approximately `2.3e-14` relative for velocity, temperature,
eta, and the checked partial densities.

The full pseudocolor result and timestep/volume history are in
`reports/lowmach_stress/lowmach_report.html`. When its local server is active,
the report is available at `http://127.0.0.1:8765/lowmach_report.html`.

Useful regression commands are:

```bash
./scripts/runtests.py tests/LMDrivenCavity \
  --sections 2d-amr 2d-amr-parallel --comp clang++ --no-clean --no-backspace

./scripts/runtests.py tests/LowMachChemistry \
  --sections 2d 2d-implicit --comp clang++ --no-clean --no-backspace
```

## Historical Findings Still Relevant

- Do not reintroduce integrator-local state `FillPatch`, average-down, or
  regridding implementations without first demonstrating a missing framework
  operation. Several coarse/fine zeroing and asymmetric high-face failures came
  from competing state-management paths.
- Velocity is genuinely cell centered in the current branch. Treating nodal
  velocity as cell centered caused high-side stencil omissions; trying to make
  every downstream operator nodal produced a much larger inconsistent system.
- A hierarchy-wide projection is required. Independent level solves produced
  coarse/fine diffusion and vorticity errors.
- Hard resetting or snapping `xi` creates artificial gradients. Reconstruction
  must extend material history continuously into the diffuse exterior.
- A sharp stress cutoff at `eta=0.5` creates a stress jump. Stress evaluation
  and mechanical weighting are separate: evaluate on the reconstructed
  extension, then apply a smooth mechanical weight.
- Applying full elastic Cauchy stress divergence alongside the pressure
  projection double counts the pressure-like volumetric response. The current
  explicit feedback is deviatoric only.
- Lower elastic wave speeds, Newtonian viscosity, interface viscosity, WENO5,
  and stress-map smoothing can damp symptoms, but none repairs an inconsistent
  reference map or thermodynamic interface source.
- The removed mapped-distance/counter-curvature reinitialization could preserve
  a contour only by adding nonconservative exterior cleanup and still produced
  corner-specific behavior. The current branch intentionally uses the simpler
  Allen-Cahn mechanism instead.
- The removed implicit elastic solve was not a general solve for arbitrary
  constitutive models. Do not count it as implicit mechanics if revisiting that
  direction.

## Remaining Risks And Next Steps

1. Add condensed-species heat capacities and phase-change enthalpies. The
   present use of gas `cp` for the entire mixture is only a stabilizing first
   model.
2. Use one projected face-velocity field for both scalar fluxes and the
   projection. The current one-step volume discrepancy correction leaves about
   five percent pointwise volume error in the 100-microsecond AP/HTPB case.
3. Revalidate `input.lm.couette_solid` and `input.lm` for long times after the
   species, projection, and chemistry refactors. The most recent exhaustive
   validation focused on pure flow and rigid-phase combustion.
4. Quantify conservation across AMR regrids for each partial density, not only
   serial/MPI agreement at a fixed time.
5. Add a focused deformable-solid regression that checks reference-map
   deformation, stress, and mechanical response rather than relying on visual
   Couette output.
6. Add a rigid-effusion regression that checks condensed mass loss, gas product
   gain, and the prescribed rigid velocity.
7. Generalize beyond one deformable solid and one shared reference map if
   multiple mobile condensed bodies or phases are required.
8. Generalize phase change to permit positive eta rates for condensation or
   growth while preserving bounds and mass conservation.
9. Add calibrated AP and HTPB surface kinetics. The current rate multipliers
   intentionally accelerate regression for testing.
10. Revisit the explicit elastic timestep only after the spatial stress and
    reference-map discretizations have dedicated regressions. Damping cannot
    remove the elastic wave CFL in general.
11. Run and validate the complete method in 3D. Compilation alone is not a
    mechanics, chemistry, or AMR validation.
12. Keep the LowMach hot path compact. Diagnostics and temporary checks should
    not become permanent state unless the evolution equations consume them.

## Working Conventions

- Partial densities are the definitive species state.
- `eta` is reconstructed from condensed partial density; do not independently
  advect or snap it.
- Distinguish physical terms from numerical processing. Reference-map repair
  belongs in `Numeric::ReferenceMap`; constitutive behavior belongs in
  `Model::Solid`; phase kinetics belongs in `Model::PhaseField` and
  `Model::Mechanism`; projection setup belongs in `Operator::PressurePoisson`.
- Prefer the existing Integrator and AMReX hierarchy communication machinery.
- Keep pure-fluid behavior covered by `LMDrivenCavity` while changing mixed
  mechanics or chemistry.
- Keep explicit and implicit chemistry covered by `LowMachChemistry`.
- The target remains a diffuse-boundary LowMach method, not a level-set rewrite.
