# LowMach Eulerian Solids Findings

Last updated: 2026-07-22

This file is a technical handoff for the current `eulerian-solids` LowMach
implementation. It describes code active in the present branch, results
verified through 2026-07-22, and remaining modeling gaps. Older mapped-distance
reinitialization, nodal-velocity, custom regridding, and implicit-elastic
experiments are historical only and are identified as such below.

The principal files are:

- `src/Integrator/LowMach.H` and `src/Integrator/LowMach.cpp`
- `src/Operator/Diffusion.H` and `Diffusion.cpp`
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
- composite AMR implicit species, thermal, and optional momentum diffusion;
- a hierarchy-wide variable-coefficient pressure projection;
- selectable collocated advection operators; and
- AMR refinement based on velocity, pressure, temperature, or reconstructed
  solid volume fraction.

The pure-fluid and chemistry regressions pass in serial/MPI AMR configurations.
The AP/HTPB MPI/AMR case reaches 500 microseconds without the former diffuse-
interface temperature singularity or diffusion timestep collapse. Deformable-solid
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

With explicit viscosity, the momentum predictor contains

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
When `implicit_viscosity=1`, the Newtonian Laplacian is removed from this RHS
and advanced by the composite diffusion solve described below. Elastic stress
and solid/interface damping remain explicit.

Partial densities use conservative advection,

```text
partial_rho_dot = -div(partial_rho * u) + mechanisms + diffusion + chemistry,
```

while velocity, temperature, and the reference map use the advective form.
SSPRK3 is used by the current principal inputs. Chemistry may additionally use
Strang splitting around the transport update.

## Implicit Transport Diffusion

`Operator::Diffusion` owns the AMReX face coefficients, hierarchy fill, MLMG
operator, and synchronization needed for a composite variable-coefficient
backward-Euler solve. It advances a cell-centered field `q` according to

```text
a * (q_new - q_old) - dt * div(b * grad(q_new)) = 0.
```

LowMach uses this same operator with different local coefficients:

```text
species mass fractions:  a = rho_g,          b = rho_g * D
temperature:             a = rho_mixture*cp, b = conductivity
velocity (optional):     a = rho_mixture,    b = dynamic viscosity.
```

The species solve currently requires a transport model with one common
diffusivity. It solves all gas mass fractions together and then reconstructs
partial densities using the pre-diffusion gas density, preserving the local
sum of gas partial densities. Thermal and species diffusion contribute their
integrated EOS volume change to the following pressure projection.

Thermal conduction and common-coefficient Rocfire species diffusion are
implicit automatically when active. Newtonian momentum diffusion is opt-in
through `implicit_viscosity=1`; this preserves the established explicit cavity
discretization for existing inputs. The AP/HTPB input enables all three.

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
identifies one condensed input species and one or more gas output species. For
a rigid input, Allen-Cahn acts on the total rigid volume fraction so that only
the exterior rigid/gas interface evolves. The resulting eta rate is allocated
to the input species according to its local share of rigid volume. This avoids
both spurious reaction at internal condensed-species interfaces and pinning of
the dilute exterior tail. For a deformable input, Allen-Cahn acts directly on
that species' reconstructed volume fraction.

The eta rate is converted to a mass source, that mass is subtracted from the
condensed species, and the same mass is distributed among the gas products.
Product mass fractions must sum to one.

The current phase-change implementation deliberately applies

```text
eta_dot = min(eta_dot_AllenCahn, 0).
```

It therefore supports regression, decomposition, sublimation, or pyrolysis,
but not condensation or phase growth. No additional interface mask is applied:
the Allen-Cahn operator already contains the required interfacial degeneracy.
Optional temperature cutoff, Arrhenius activation, and a rate multiplier are
available.

Mass transfer is conservative even though Allen-Cahn is not conservative in
eta. The mechanism also supplies the projection with the volume source

```text
mass_source * (1/rho_ref_condensed - R_product*T/p0).
```

The deformable-solid demonstration inputs use zero driving force mainly to
maintain a regular profile while transferring any lost solid mass to the gas
species. The AP/HTPB input uses large artificial rate multipliers so regression
is observable over a short test. Those multipliers are testing parameters, not
calibrated propellant kinetics.

An earlier implementation multiplied the Allen-Cahn rate by an additional
`4*eta*(1-eta)` interface mask. Since Allen-Cahn's chemical potential already
degenerates at the pure phases, the extra factor made the dilute gas-side tail
move as eta squared. The `eta=0.5` surface regressed while that tail, and the
associated heating front, remained near its initial location. The extra mask
and its `interface_only` input were removed. The AP/HTPB initial condition now
uses the stationary `0.5*(1-tanh(2*y/epsilon))` profile corresponding to the
configured `sigma,epsilon` coefficients instead of an unmatched error-function
profile.

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

Condensed species may define `specific_heat` and `thermal_conductivity`. When
these are present, every condensed species must define both and the implicit
temperature solve uses

```text
a = rho_g * cp_gas + sum_s(rho_s * cp_s)
b = alpha_g * k_gas(T) + sum_s(eta_s * k_s).
```

Here `rho_s` is condensed partial density and `eta_s=rho_s/rho_ref,s`. Inputs
that define no condensed thermal properties retain the previous gas-property
closure exactly. There is no eta cutoff in either path, and the composite MLMG
operator conducts across the diffuse interface.

`input.lm.ap_htpb` uses the established Flame propellant values: AP has
`k=0.4186 W/(m K)` and `cp=1297.90 J/(kg K)`; HTPB has
`k=0.1300 W/(m K)` and `cp=2418.29 J/(kg K)`. Their corresponding pure-solid
diffusivities are `1.65e-7` and `5.84e-8 m^2/s`. The former gas-property
closure gave `2.64e-8` and `5.60e-8 m^2/s`, underdiffusing AP by about 6.3
times and reversing the AP/HTPB ordering.

In the physical scalar-conduction run at 6 ms, cells with raw condensed volume
fraction between 0.45 and 0.55 span approximately 996--1445 K. The configured
solid coefficients are therefore active and the gas flame supplies substantial
conductive feedback without an artificial bridge.

The checked-in AP and HTPB mechanisms now have distinct Arrhenius parameters.
They remain provisional because the AP fit used the removed tensor bridge.
LowMach still has no phase-change latent-enthalpy source, so quantitative
grooving and burn-rate claims require the conservative enthalpy work described
below.

### AP monopropellant pressure calibration

Figure 5a of `main.pdf` gives approximate pure-AP experimental rates of zero
through 1.5 MPa and 2.65, 3.8, 5.5, and 7.9 mm/s at 2, 3, 4.5, and 6 MPa. The
paper uses a temperature-only Arrhenius mobility with a 11000 K activation
temperature and no direct pressure factor.

A previous calibration sweep used the normal tensor thermal bridge and appeared
to recover a useful pressure-dependent interface temperature. The curved
AP/HTPB run later showed that the bridge does not preserve temperature bounds,
so that calibration is no longer accepted. The phase-change coefficients must
be recalibrated against the physical scalar conduction model.

The provisional phase-change law retained in the inputs is

```text
rate_multiplier = 2450
activation_temperature = 3145 K.
```

Pressure still does not enter the mobility explicitly. It changes the Rocfire
reaction structure, normal heat flux, interface temperature, and hence the
Arrhenius factor. The multiplier is specific to the current Allen-Cahn
normalization (`mobility=0.01 1/Pa/s`, `sigma=0.001 J/m2`, `epsilon=20 um`, and
`driving_force=200 Pa`) and is not directly comparable to the paper's
pre-exponential factor.

The historical bridge-based sweep used an 800 micrometer column with one
eta-tracked AMR level,
giving isotropic 3.125 micrometer finest cells through the interface, and runs
for 20 ms. The 10--20 ms average rates at 2, 3, 4, and 6 MPa are 2.425, 3.392,
4.727, and 6.997 mm/s. Rates over the final two-millisecond interval are 2.530,
3.853, 5.590, and 7.914 mm/s. The latter compare with approximate experimental
values of 2.6, 3.7, 5.0, and 7.9 mm/s in Figure 5a. The corresponding final
`eta=0.5` temperatures are 778.3, 878.1, 990.2, and 1119.6 K.

The former full-tensor heat-flux values are bridge-dependent and are not valid
physical calibration observables. A replacement comparison must use the
physical scalar conductive flux and a conservative energy balance across the
diffuse layer.

The 1.5 MPa extinction point is outside the paper's stated 2--6 MPa calibration
range and is not reproduced by this continuously active law. Reproducing it
requires an ignition/extinction construction; it should not be forced by
adding pressure directly to the mobility.

### HTPB sandwich calibration

The original HTPB law (`A=2.333e5`, `Ta=7500 K`, cutoff `500 K`) produced
essentially no HTPB regression or binder gas in the 3 MPa sandwich. This
starved the gas chemistry of binder reactant and contributed to the observed
intermittent burning. With the same Allen-Cahn normalization used for AP, a
short 3 MPa sweep selected

```text
rate_multiplier       = 1.85e5
activation_temperature = 3145 K
temperature_cutoff     = 360 K
```

The resulting HTPB `eta=0.5` front regressed at 3.95 mm/s over 0.75--1.5 ms
and 4.07 mm/s over 0.9--1.5 ms. The AP front moved at 5.71 mm/s over the latter
window. This is an initial 3 MPa, short-time calibration for continuous binder
supply; it is not yet a validation of long-time sandwich morphology or of the
HTPB pressure response.

The sweep also exposed a numerical issue in the implicit phase-field solve.
Inactive cells had been assigned a very large artificial mass coefficient.
That coefficient dominated MLMG's composite residual normalization and could
make an active fine level return exactly zero phase update. The solve now
scales each equation by the local gradient coefficient with a floor based on
the global active maximum, and applies the phase change only where the physical
coefficient exceeds that floor. Thermal and species diffusion are applied
before phase change so the Arrhenius law sees the current diffused temperature.

That historical sweep was a 20 ms front-rate fit, not an asymptotic steady-state
fit, and must not be treated as validation of the current scalar model.
The 3 and 4 MPa fronts continue to accelerate after 20 ms. The NaNs formerly
seen during extension runs came from mole-fraction normalization, not AMR patch
motion. Species diffusion creates positive gas-density tails that eventually
become subnormal in the condensed region. `Model::Gas::MoleFraction` previously
stored the reciprocal of total molar density; that reciprocal overflowed, and
zero species densities multiplied by infinity produced NaNs in `cp_mass` and
then in the implicit thermal solve. Dividing each molar density directly by the
total avoids the overflow. A strict 6 MPa extension now completes through 30 ms
with all 30 plot fields finite; its 28--30 ms front rate is 8.017 mm/s, compared
with 7.914 mm/s over 18--20 ms. The two-process implicit LowMach chemistry
regression also passes. The plotfile diagnostic now advances field names by
component count, so the former `velocityy` message is correctly identified as
temperature. No NaN check or solver tolerance was relaxed. The transverse
domain was widened so the AMR-refined periodic and normal spacings are both
3.125 micrometers; the earlier anisotropic calibration layout is discarded.

The Gross-Beckstead comparison requires more care than the raw values suggest.
Their reported heat flux is calculated on a geometrically sharp propellant
surface by a steady, detailed multicomponent gas-phase model. The gas-side
temperature gradient and CHEMKIN conductivity supply the local Fourier heat
feedback. Species diffusion at the surface is included in the gas solution and
coupled boundary conditions, while the Dufour effect is explicitly neglected.
The papers do not document a separate species-enthalpy contribution to the
plotted flux. The resulting conductive heat flux is fed into
precomputed one-dimensional AP or AP/HTPB gas/condensed correlations. Those
correlations update surface temperature, injected species, and mass flux; the
two-dimensional gas solution and surface state are iterated to convergence.
At 20 atm they report 12.45 MW/m2 for pure AP and 12.46 MW/m2 at the AP-particle
centerline. Their smallest gas-side surface cell is 0.008 micrometers.

By contrast, the former bridge values were full-tensor fluxes at the center of a
roughly 20-micrometer diffuse interface on a 3.125-micrometer mesh. They are not
the same observable as the sharp gas-side Gross-Beckstead flux. Equation (15a)
and Table 2 of our paper fit the Gross-Beckstead calculation as
`q = 1e7*(0.46*p + 0.42) W/m2`, giving 13.4 MW/m2 at 2 MPa and 31.8 MW/m2 at 6
MPa, but the remaining difference from the current `eta=0.5` values cannot be
interpreted as a validated heat-transfer deficit.
A like-for-like diagnostic must either extrapolate the pure-gas conductive
flux to the representative surface or perform a conservative energy balance
across a pillbox containing the diffuse layer and phase change. The latter must
also account for species energy transport even though it is not part of the
documented Gross-Beckstead Fourier-flux postprocessing.
These values and the corrected interpretation are recorded in
`reports/lowmach_stress/lowmach_report.html`.

The outward gas flow does not invalidate conduction as the surface-feedback
mechanism. In a steady one-dimensional gas layer, reaction supplies energy,
outward enthalpy advection removes part of it, and the remainder can conduct
down the gas-side temperature gradient toward the propellant. The relevant
measure is the thermal Peclet number `Pe=u*L/alpha`, or equivalently the upstream
thermal penetration length `L_th=alpha/u`.

The common-temperature advection closure transports only the gas heat capacity
at the gas velocity. For a stationary rigid solid and moving gas under local
thermal equilibrium, the equation is

```text
C_mix * dT/dt + rho_g*cp_g*u_g*grad(T) = div(k_mix*grad(T)) + Q.
```

At the previously sampled midpoint,
`rho_g*cp_g/C_mix` was approximately 0.005, making this weighting essential.
The current implementation applies this ratio to the advective temperature
term.

### Tensor thermal bridge failure and removal

The `Pe=0.01` tensor bridge passed planar tests but failed on the curved HTPB
edges in `output.lm.ap_htpb`. At 4.9703 ms the solution still had
`Tmin=700.01 K`. It then fell to 594.8 K at 4.9796 ms and 44.85 K at
4.9898 ms. The two cold cells were symmetric, near
`x=+/-0.255 mm, y=-0.013 mm`, inside the diffuse HTPB interface where its
normal is oblique to the grid. The chemistry model is exothermic, so this was a
numerical undershoot rather than chemical cooling. The local chemistry solve
failed shortly afterward.

The bridge conductivity there was hundreds of W/(m K), compared with physical
conductivities of order 0.1 W/(m K). Although the continuum tensor was positive
definite, its strong off-diagonal cross derivatives did not give a
maximum-principle-preserving discrete operator. Positive definiteness alone was
therefore insufficient at anisotropy ratios of thousands.

A fresh run with the bridge disabled completed through 6 ms. It retained
`Tmin=700.04 K`, reached `Tmax=2941 K`, and had temperatures of approximately
996--1445 K in cells with raw condensed volume fraction between 0.45 and 0.55.
This demonstrates that the physical scalar mixture conductivity transfers
substantial heat without an artificial interface boost.

LowMach no longer parses or applies `thermal_bridge`, and the AP inputs no
longer configure it. The implicit temperature equation now always uses

```text
C_mix * (T_new-T_old)/dt = div(k_mix * grad(T_new)) + sources,
k_mix = alpha_g*k_g + sum_s(alpha_s*k_s).
```

Its harmonic face conductivity produces conservative two-point fluxes and
preserves temperature bounds for positive heat capacity and conductivity. The
generic tensor capability remains in `Operator::Diffusion`, but it is no longer
part of diffuse-interface thermal transport.

The next thermodynamic improvement should be conservative mixture enthalpy,
including sensible and latent enthalpy carried by phase change, followed by
interface-width convergence. Artificially increasing conductivity or directly
overwriting/extrapolating temperature should not be used to compensate for a
missing energy term.

Radiation cannot supply the missing pure-AP heat feedback. A deliberately
generous blackbody bound,
`sigma_SB*(T_flame^4-T_surface^4)`, is only 0.168 MW/m2 for the 1350 K AP flame
and 773 K surface reported by Gross and Beckstead at 20 atm. This is 1.35% of
their 12.45 MW/m2 calculated surface flux. Even a unit-emissivity, unit-view-
factor 3000 K source gives only 4.57 MW/m2; an actual micrometer-scale gas flame
has lower effective emissivity and view factor. Radiation is therefore not a
plausible replacement for conduction in the current pure-AP or nonmetalized
AP/HTPB calibration.

Aluminum changes this conclusion from negligible to potentially important.
Burning Al and Al2O3 particles provide hot particulate emitters. A 1991
experiment found radiation negligible in nonmetalized propellants but equal to
26% of total heat feedback for a 20%-aluminum formulation at 1 MPa. Later work
also reports downstream aluminum combustion increasing the surface burning rate
through radiative feedback, in competition with the inert aluminum heat sink.
The eventual AP/HTPB/Al model should therefore include radiation after aluminum
ignition and particle transport are represented, but the evidence supports a
significant correction rather than radiation as the general dominant feedback.

The thin `AP_gas` and `HTPB_gas` layers are not passive species that can be
removed by redirecting phase change to `Mono` and `Premixed`. The six-species
Rocfire mechanism represents four flames:

```text
AP_gas -> Mono
HTPB_gas -> Premixed
beta AP_gas + HTPB_gas -> (beta+1) Primary
gamma Mono + Premixed -> (gamma+1) Final.
```

Each pathway has independent Arrhenius and pressure-dependent kinetics. Direct
production of `Mono` and `Premixed` would disable the first three pathways. In
the pure-AP case, no reaction would remain because the final pathway requires
both `Mono` and `Premixed`. It would also remove the net 330 cal/g AP and 147
cal/g HTPB heat release currently attached to the first two reactions. The
generic phase-change mechanism supplies mass and low-Mach volume change but no
thermal source, so an input-only product substitution would silently discard
that energy.

An instantaneous source reconstruction from the AP/HTPB output at 8.0708 ms
gave area-averaged heat releases of 10.42, 0.89, 1.48, and 11.19 MW/m2 for the
AP-to-Mono, HTPB-to-Premixed, primary, and final pathways respectively. The
three precursor-dependent pathways therefore supplied about 53% of the
instantaneous reconstructed total. Their spatial support can be narrow while
their dynamical contribution remains large.

A reduced surface-flame closure is possible as a separate model. It must add
the eliminated heat at a physically chosen location, retain the effective
pressure/temperature response, and replace the primary flame. The last item is
not a simple local reduction: competition between self-reaction and lateral
AP/HTPB mixing is the mechanism that creates the primary diffusion flame and
particle-size effects. Keep the existing six-species Rocfire model as the
reference implementation; test any direct solid-to-product closure under a
distinct model name and compare heat flux, flame standoff, grooving, and
pressure response.

## Verified Results

The following results include the 2026-07-22 regressions and the scalar thermal
investigation on 2026-07-23 with the 2D clang build.

| Case | Configuration | Result |
| --- | --- | --- |
| `LMDrivenCavity` | Re=100, AMR levels 0-1, serial | Run and reference-profile check passed |
| `LMDrivenCavity` | Re=100, AMR levels 0-1, MPI 2 | Run and reference-profile check passed |
| `LowMachChemistry` | 29-reaction H2/O2, explicit, MPI 2, AMR | Product, positivity, pressure, and AMR checks passed |
| `LowMachChemistry` | Same case, implicit with 10x flow timestep | Product, positivity, pressure, and AMR checks passed |
| `input.lm.ap_htpb` | Physical scalar conduction, MPI 8, AMR levels 0-2 | Reached 6.0006 ms in 7154 steps |

Near-6-ms AP/HTPB metrics were:

```text
flow timestep            0.604 microseconds (advective limit)
temperature range        700.04--2941.12 K
maximum speed            4.66 m/s
pressure range           2.999--5.381 MPa
```

At 20.176 microseconds, serial and two-rank fields agreed to at most
approximately `2.3e-14` relative for velocity and `2.6e-15` for the checked
thermochemical fields. At 500 microseconds, temperature, pressure, velocity,
eta, and all partial densities were finite; covering-grid checkerboard metrics
were below `0.0027` after normalization by each field range.

At the AP centerline, the corrected `eta=0.5`, `eta=0.0033`, and 1000 K
locations moved from approximately `0`, `28.68`, and undefined micrometers at
initialization to `-1.87`, `27.24`, and `27.80` micrometers at 500 microseconds.
In the defective run the core moved to `-3.16` micrometers, but the dilute tail
remained at `38.45` micrometers and the 1000 K contour remained at `40.62`
micrometers. A longer corrected output remained finite through 3.09
milliseconds, with the three locations at `-11.59`, `17.56`, and `18.19`
micrometers respectively.

The checked-in AP/HTPB settings use `cfl=0.6` and
`dynamictimestep.max=2.5e-6`. A sweep through CFL 1.5 found that larger values
could remain finite but accumulated increasing operator-splitting error; CFL
1.5 also increased velocity and checkerboard measures. At 200 microseconds,
CFL 0.6 differed from a CFL 0.125 reference by 0.78% in temperature, 2.36% in
speed, and 1.06% in the aggregate gas-species field. A matched two-rank run to
500 microseconds took 16.37 seconds, compared with 37.53 seconds using the old
CFL 0.25 and 1-microsecond cap.

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

1. Evolve conservative mixture enthalpy rather than temperature so changing
   phase fractions carry sensible energy consistently.
2. Verify scalar thermal convergence as the phase-field width is reduced on
   planar and smooth curved interfaces.
3. Add phase-change latent enthalpy and interfacial heat feedback. Condensed heat
   capacities and conductivities are now species properties, but phase change
   still transfers mass without a corresponding thermal source.
4. Use one projected face-velocity field for both scalar fluxes and the
   projection. The current one-step volume discrepancy correction leaves about
   five percent pointwise volume error in the 100-microsecond AP/HTPB case.
5. Revalidate `input.lm.couette_solid` and `input.lm` for long times after the
   species, projection, and chemistry refactors. The most recent exhaustive
   validation focused on pure flow and rigid-phase combustion.
6. Quantify conservation across AMR regrids for each partial density, not only
   serial/MPI agreement at a fixed time.
7. Add a focused deformable-solid regression that checks reference-map
   deformation, stress, and mechanical response rather than relying on visual
   Couette output.
8. Add a rigid-effusion regression that checks condensed mass loss, gas product
   gain, and the prescribed rigid velocity.
9. Generalize beyond one deformable solid and one shared reference map if
   multiple mobile condensed bodies or phases are required.
10. Generalize phase change to permit positive eta rates for condensation or
   growth while preserving bounds and mass conservation.
11. Recalibrate the distinct AP and HTPB Arrhenius surface kinetics with scalar
   conduction after the conservative enthalpy equation is in place.
12. Revisit the explicit elastic timestep only after the spatial stress and
    reference-map discretizations have dedicated regressions. Damping cannot
    remove the elastic wave CFL in general.
13. Run and validate the complete method in 3D. Compilation alone is not a
    mechanics, chemistry, or AMR validation.
14. Keep the LowMach hot path compact. Diagnostics and temporary checks should
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
