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

At 3 milliseconds in the updated AP/HTPB output, temperature on the `eta=0.5`
contour is about 700.7 K over AP and 701.6 K over HTPB. Heating is visible
farther into the gas-side tail: on `eta=0.01`, the corresponding temperatures
are about 784.7 K and 789.9 K. An isolated conduction probe, with chemistry,
regression, and temperature advection disabled, heated the initially cold
`eta=0.9` solid to 713.9 K over AP and 716.0 K over HTPB in 3 milliseconds.
This verifies that the configured solid coefficients are active. The modest
surface heating in the coupled case is now attributable to missing interfacial
heat feedback rather than to gas properties or an eta cutoff.

The two configured phase-change mechanisms are also exactly equivalent:
mobility, surface energy, interface width, driving force, rate multiplier, and
temperature cutoff are identical, while activation temperature is zero. Their
volume regression rates must therefore be identical regardless of the
temperature field. Unlike the Flame `FullFeedback` model, LowMach currently
has no AP/HTPB-specific Arrhenius prefactors or activation temperatures and no
interfacial heat-feedback or phase-change enthalpy source. Grooving is not
expected from the checked-in input until these material couplings are added.

### AP monopropellant pressure calibration

Figure 5a of `main.pdf` gives approximate pure-AP experimental rates of zero
through 1.5 MPa and 2.65, 3.8, 5.5, and 7.9 mm/s at 2, 3, 4.5, and 6 MPa. The
paper uses a temperature-only Arrhenius mobility with a 11000 K activation
temperature and no direct pressure factor.

The effectively one-dimensional `input.lm.ap_monopropellant` case was run at
1.5, 2, 3, and 6 MPa. With the current temperature-independent kinetics, the
computed rates were 3.754, 3.776, 3.772, and 3.776 mm/s. The corresponding
late `eta=0.5` temperatures were only 700.10, 700.26, 700.70, and 702.72 K.
Normalizing the paper's 11000 K activation law to the 3 MPa rate changed the
2--6 MPa prediction only from 3.77 to 3.88 mm/s. A second 4050 K candidate
similarly gave 3.73 to 3.78 mm/s.

Formally fitting an Arrhenius law to the center temperatures and experimental
rates requires an activation temperature near 205000 K and a prefactor of
order `1e129`. That fit is neither physical nor numerically usable: temperature
in the gas-side diffuse tail reaches 725--901 K, so the fitted mobility would
become enormous there and destroy the interface profile. A localized ignition
pulse and temperature cutoffs of 701--702 K were also tested; both 2 and 6 MPa
cases settled back to almost the same regression rate.

The present model therefore cannot be calibrated to the experimental pressure
curve by changing phase-field parameters alone. The missing constraint is a
pressure-dependent thermal feedback that produces a meaningful surface-
temperature response. The paper supplies this through its mass-flux/pressure
surface-heat-flux model; a fully coupled calculation must instead recover the
same feedback from gas reaction, interface heat transfer, and phase-change
enthalpy. Do not introduce pressure directly into the phase-change mobility to
force this fit.

The apparently inverted gas-temperature trend has now been resolved. At 5 ms,
the 2 MPa case has a maximum temperature of 1337.9 K approximately 114
micrometers above its moving `eta=0.5` surface, while the 6 MPa case has a
maximum of 1056.6 K approximately 74 micrometers above the surface. The
pressure-dependent Rocfire rate moves the high-pressure reaction zone into the
diffuse solid/gas mixture: 46.2% of the effective 6 MPa heat release occurs at
`eta>0.01` and 7.6% at `eta>0.5`, compared with 7.3% and 1.5% at 2 MPa. The
heat-release-weighted eta changes from 0.0195 to 0.1018. Consequently, more of
the high-pressure reaction heats the mixed condensed/gas heat capacity: the
heat-release-weighted condensed share of local heat capacity is 53.7% at 6 MPa
and 22.3% at 2 MPa. The gas temperature peak is therefore lower and closer to
the surface.

The conductive heat flux into AP was reconstructed using the exact thermal
coefficient used by the implicit solve,
`k_mix = alpha_g*k_g + eta*k_AP`, followed by the same harmonic cell-to-face
average. At the face crossing `eta=0.5`, `k_mix*dT/dy` is 0.00417 MW/m2 at 2
MPa and 0.0641 MW/m2 at 6 MPa. Thus the high-pressure case has 15.4 times more
interface heat feedback despite its lower maximum temperature. Near
`eta=0.01`, the corresponding fluxes are 0.431 and 2.19 MW/m2. The effective
Rocfire heat release integrated through the domain is 13.33 and 11.71 MW/m2.

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

By contrast, the values above are mixed-coefficient fluxes at the center of a
roughly 20-micrometer diffuse interface on a 3.125-micrometer mesh. They are not
the same observable as the sharp gas-side Gross-Beckstead flux. Equation (15a)
and Table 2 of our paper fit the Gross-Beckstead calculation as
`q = 1e7*(0.46*p + 0.42) W/m2`, giving 13.4 MW/m2 at 2 MPa and 31.8 MW/m2 at 6
MPa, but the apparent factors of 3200 and 500 relative to the current
`eta=0.5` values cannot be interpreted as a validated heat-transfer deficit.
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

The final 3 MPa monopropellant output quantifies the issue in the present
diffuse model. At the gas-side `eta=0.01` contour, `T=772.4 K`, `u_y=0.596 m/s`,
`rho_g=12.03 kg/m3`, `k_g=0.0711 W/(m K)`, and
`alpha_g=4.71e-6 m2/s`. Thus `alpha/u=7.9 micrometers`. The distance from
`eta=0.5` to `eta=0.01` is 23.1 micrometers, giving `Pe=2.92`. The local mass
flux `rho_g*u_y=7.16 kg/m2/s` is consistent with the approximately
`rho_AP*r_b=7.35 kg/m2/s` supplied by regression, so the outward blowing is
physical. Conduction can oppose it over several micrometers, but thermal
feedback is strongly attenuated over the much wider diffuse tail. Gross and
Beckstead's sharp-interface AP result rises from 773 K to about 1350 K within 8
micrometers at 20 atm, a scale comparable to the calculated penetration length.

This points to an interface-resolution/coupling problem rather than the wrong
physical heat-transfer mechanism. A sharp or asymptotically consistent diffuse
formulation must place the representative surface within the gas thermal
preheat length and transmit the gas-side Fourier flux into the condensed energy
balance without forcing it through a many-penetration-length mixed layer.

The interface is not primarily a low-conductivity layer. At the sampled
`eta=0.572` cell in the same final profile, `rho_AP=1115.5 kg/m3` and
`rho_g=5.22 kg/m3`. The AP volumetric heat capacity is approximately
`1.45e6 J/(m3 K)`, while the gas contributes only about `0.007e6 J/(m3 K)`.
The mixed conductivity is approximately `0.24 W/(m K)`, which is greater than
the gas conductivity, but the resulting diffusivity is only about
`1.6e-7 m2/s`, roughly thirty times below the gas value. A localized
conductivity boost can therefore transmit heat through the finite-width
thermal mass, and it remains conservative when implemented inside
`div(k*grad(T))`.

Before using such a boost as the model, correct the temperature-advection
closure. LowMach currently applies `-u*grad(T)` to the common temperature at
every mixture point. For a stationary rigid solid and moving gas under local
thermal equilibrium, the corresponding energy equation is instead

```text
C_mix * dT/dt + rho_g*cp_g*u_g*grad(T) = div(k_mix*grad(T)) + Q.
```

At the sampled midpoint, `rho_g*cp_g/C_mix` is approximately 0.005. The current
equation therefore transports the condensed thermal inertia at the gas
velocity and can remove heat from the interface much too strongly. This is a
more fundamental inconsistency than the conductivity interpolation.

After that correction, a useful diagnostic is
`k_eff=k_mix+k_bridge*4*eta*(1-eta)`. The preferable production form is a
normal-only tensor bridge,
`K=k_mix*I+k_bridge*4*eta*(1-eta)*n*n`, so heat crosses the diffuse layer
without smoothing tangential AP/HTPB temperature structure. Any bridge must be
tested over multiple phase-field widths and calibrated to remove width
dependence. Directly overwriting or extrapolating the physical temperature is
not conservative and should not be used as the energy update.

The bridge amplitude should explicitly reference the configured phase-field
width `ell`; using `grad(eta)` directly for its magnitude is unnecessarily
noisy and becomes singular where the phase is uniform. The gradient is useful
for constructing the interface normal. The width scaling depends on the
constraint imposed on the artificial layer:

```text
fixed crossing time tau:       k_bridge ~ C_int * ell^2 / tau
fixed interface Peclet Pe:     k_bridge ~ C_int * U * ell / Pe
fixed contact resistance R:    k_bridge ~ ell / R.
```

Consequently, there is no unique gradient-based correction until the desired
sharp-interface limit is stated. For the present goal, a small prescribed
interface Peclet number is the closest expression of making the diffuse layer
transparent to heat while gas blows away from the surface. The enhancement
then decreases linearly with width as the sharp-interface limit is approached.

The positive coefficient is handled by the backward-Euler MLMG diffusion solve
and adds no explicit diffusion CFL restriction. The initial scalar bridge
confirmed that normal heat transfer could be restored, but it also held the
entire AP/HTPB interface at nearly one temperature. Supplying both
`thermal_bridge.width` and `thermal_bridge.peclet` now gives the full tensor

```text
K = k_mix*I
  + 4*eta*(1-eta)*C_mix*|u dot n|*width*(n tensor n)/peclet,
n = grad(eta)/|grad(eta)|.
```

The physical mixture conductivity remains isotropic. Only the artificial
bridge is projected onto the reconstructed condensed-interface normal. The
feature is inactive unless both inputs are present, and parsing rejects
incomplete, nonpositive, or non-condensed configurations.
`input.lm.ap_monopropellant` and `input.lm.ap_htpb` use `width=20 um`, matching
their Allen-Cahn epsilon, and `peclet=0.01`.

`Operator::Diffusion` now accepts either scalar or full tensor cell mobility.
The total tensor is transferred to each face using the matrix harmonic mean
`2*(K_lo^-1+K_hi^-1)^-1`. This preserves positive definiteness and reduces to
the former scalar harmonic average for grid-aligned isotropic coefficients.
The diagonal face block is supplied to `MLABecLaplacian`; the extended operator
adds every off-diagonal face flux to `Fapply`, uses the full residual in a
weighted-Jacobi smoother, and returns the same full tensor flux through
`FFlux`. Consequently, composite-AMR reflux and residual evaluation use the
same conservative discretization. Every face-tensor row is also coarsened
through the multigrid and AMR hierarchy.

Strong directional ratios use eight pre/post smoothing sweeps, and aligned
cases can use AMReX semicoarsening. The requested `1e-11` relative tolerance is
unchanged, and failed solves are not accepted. `Kxy`, `Kxz`, and `Kyz` are now
represented rather than projected out.

A planar AP sweep through bridge Peclet numbers off, 1, 0.1, and 0.01 produced
late `eta=0.5` temperatures of 700.70, 701.46, 709.58, and 768.10 K. At
`Pe=0.01`, the reconstructed center flux is 5.57 MW/m2, compared with 0.012
MW/m2 unbridged. Every case completed to 5 ms and retained the same 22.12 um
`eta=0.9--0.1` width. The temperature-independent phase-change settings also
left the fitted regression rate unchanged at 3.772 mm/s, as expected.

The final directional planar AP case completed to 5.0016 ms in 2349 steps with
a 2.197 us final timestep. Its final temperature differs from the former scalar
bridge by at most 0.009 K, which is the expected equivalence for a planar
grid-aligned normal.

The final matrix-harmonic tensor rerun gives `T(eta=0.5)=768.098 K` at 5 ms,
matching the established 768.10 K planar value. Its transverse spread is below
`1e-9 K`.

The final directional AP/HTPB case completed to 3.0007 ms on two MPI ranks and
AMR levels 0--2 in 1395 steps, with a 2.138 us final timestep. Every composite
solve converged through repeated regrids. On the interpolated `eta=0.5`
surface, the scalar bridge produced 763.67--765.74 K with a 0.64 K standard
deviation. The directional bridge produced 718.98--825.00 K with a 29.65 K
standard deviation while retaining a similar mean temperature. The AP/HTPB
mean surface-temperature separation increased from 1.37 K to 24.02 K. Thus the
normal transfer remains active without artificially equilibrating the complete
interface laterally.

`input.lm.ap_htpb_inclined` provides the off-diagonal regression geometry. Its
two long interfaces have slopes `+/-1`, so the former diagonal projection
reduced the artificial bridge to `(k_bridge/2)*I`. At 3 ms that approximation
gave an `eta=0.5` range of 735.19--740.45 K and only 0.07 K AP/HTPB mean
separation. The full tensor MPI-2, AMR 0--2 run completed to 3.0006 ms in 2762
steps with strict convergence. It gives 712.16--835.93 K and 29.05 K AP/HTPB
mean separation. Both versions move the interface by the same 14.8 micrometers,
isolating the change to tangential thermal transport.

LMDrivenCavity serial/MPI AMR and LowMachChemistry explicit/implicit
regressions all pass after this change. The unresolved shared-temperature
advection closure is unchanged by the bridge implementation.

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

The following results were checked on 2026-07-22 with the 2D clang build.

| Case | Configuration | Result |
| --- | --- | --- |
| `LMDrivenCavity` | Re=100, AMR levels 0-1, serial | Run and reference-profile check passed |
| `LMDrivenCavity` | Re=100, AMR levels 0-1, MPI 2 | Run and reference-profile check passed |
| `LowMachChemistry` | 29-reaction H2/O2, explicit, MPI 2, AMR | Product, positivity, pressure, and AMR checks passed |
| `LowMachChemistry` | Same case, implicit with 10x flow timestep | Product, positivity, pressure, and AMR checks passed |
| `input.lm.ap_htpb` | MPI 2, AMR levels 0-2 | Reached 501.653 microseconds in 223 steps |

Final AP/HTPB metrics were:

```text
flow timestep            1.99 microseconds (advective limit)
maximum temperature      2262.96 K
maximum speed            1.42 m/s
minimum partial density  -3.0e-16 (roundoff)
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

1. Correct common-temperature advection to transport gas enthalpy at the gas
   velocity rather than transporting the complete condensed thermal inertia at
   that velocity.
2. Verify full-tensor thermal convergence as the phase-field width is reduced
   on smooth curved interfaces.
3. Add phase-change enthalpy and interfacial heat feedback. Condensed heat
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
11. Add calibrated AP and HTPB Arrhenius surface kinetics. The current identical
   rate multipliers intentionally accelerate regression for transport testing
   and cannot produce grooving.
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
