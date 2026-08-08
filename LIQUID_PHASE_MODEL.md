# LowMach diffuse-interface and phase-change model

## Refactor log

### 2026-08-05: baseline and acceptance criteria

- Preserved the existing uncommitted equilibrium phase-change work and user
  input/test edits as the starting point for this refactor.
- Replayed the production aluminum checkpoint with individual mechanisms
  disabled.  The timestep collapse persisted when aluminum melting and
  vaporization were disabled, but disappeared when the AP/HTPB phase-change
  mechanisms were disabled.  The old dimensionless-rate Allen--Cahn
  decomposition model is therefore the identified source of the roughly
  `1e-9 s` timestep, rather than the equilibrium aluminum projection or the
  explicit capillary-wave limit.
- The refactor will keep partial densities as the authoritative conserved
  state, use no level-set or sharp-interface representation, and evaluate
  candidate diffuse-interface transport algorithms against mass, reduction,
  equilibrium, and runtime gates before selecting a default.
- Every numerical experiment, rejected model, calibrated dimensional input,
  and material change is recorded in this section as it is made.

### Starting implementation

The implementation described below predates the refactor.  Sections are
updated in place as each staged change passes its tests; obsolete behavior is
not retained as a compatibility layer.

### 2026-08-05: interfacial API and internal naming

- Renamed the model class to `MultiphaseInterface` and the implementation
  header to `MultiphaseInterface.H`.
- Replaced ambiguous internal names with
  `interfacial_volume_fraction`, `interfacial_chemical_potential`, and
  `interfacial_model`.  The public diagnostics remain
  `liquid_species_eta_<species>`, `rigid_species_eta_<species>`, and
  `gas_eta`, because those fields directly describe material occupancy.
- Replaced the old `liquid.*` input hierarchy with the forward API:
  `interface.enabled`, `interface.thickness`,
  `interface.surface_tension.*`, and `interface.liquid_solid.*`.
  `surface_tension` now applies only to liquid--liquid and liquid--gas pairs.
  Each enabled liquid--solid pair requires a positive regularization and
  exactly one of a contact angle or signed surface-energy difference.
- Built `lowmach-2d-clang++` and ran all 10 `LMLiquidCapillary` and
  `LMPhaseChange` cases.  All simulations and post-checks passed without
  permissive unused-input handling.

### 2026-08-05: interface-transport candidates

- Implemented and exercised three transport candidates on a 32-by-32
  periodic droplet translated for 500 steps: conserved partial-density
  advection, conservative Allen--Cahn profile restoration, and pairwise
  Cahn--Hilliard diffusion.
- Minmod advection conserved liquid mass but increased the integrated mixed
  volume by 23.1%.  The bounded Superbee scheme reduced that change to 5.0%
  without a separate phase equation.
- Conservative Allen--Cahn at `D=1e-3 m^2/s` generated a `1.23e-4` boundedness
  error, and at `D=1e-2 m^2/s` eventually caused the projection to fail.  At
  `D=1e-5 m^2/s` it passed mass, positivity, reduction, and long-transport
  checks, but changed the profile no more accurately than Superbee advection.
  It remains an explicitly selected option; `advected` remains the default
  because it is faster and equally accurate in the accepted regime.
- The explicit fourth-order Cahn--Hilliard candidate destabilized the pressure
  solve on the second global step at a mobility large enough to affect the
  interface.  A useful implementation requires a coupled implicit mixed
  chemical-potential solve that the current scalar diffusion operator cannot
  provide.  The failed candidate and its input mode were removed rather than
  retaining a nominal option that is either unstable or ineffective.

### 2026-08-05: dimensional phase-change closures

- Removed inferred Allen--Cahn/prescribed-speed phase-change behavior and the
  numerical maximum-volume-rate cap.  Every phase change now declares one of
  `equilibrium_enthalpy`, `arrhenius_surface_flux`, or `hertz_knudsen`.
- Equilibrium melting is an availability-constrained enthalpy projection at a
  specified transition temperature.  It consumes exactly the sensible
  superheat permitted by latent heat and can melt an entire eligible particle
  in one implicit thermal step without a phase-field timestep restriction.
- Surface kinetics use the diffuse pair measure
  `|eta_g grad(eta_0) - eta_0 grad(eta_g)|`, which integrates to unit area for
  a binary interface and vanishes at a buried condensed--condensed boundary.
- AP and HTPB now use the dimensional law
  `j=j_ref (p/p_ref)^n exp[-T_a(1/T-1/T_ref)]`, referenced at 3 MPa and 600 K
  with `n=0.4` and a common nominal regression speed of 6 mm/s.
- Aluminum evaporation now uses a locally coupled Hertz--Knudsen flux with a
  Clausius--Clapeyron saturation pressure anchored at 2743 K and 101325 Pa.
  The current aluminum accommodation coefficient is 0.04.  Latent cooling,
  vapor partial-pressure feedback, availability, and forward-only unresolved
  oxidation enthalpy are iterated locally without constraining the global
  timestep.
- Coarse/fine Hertz--Knudsen tests produced `1.64359e-2` and `1.65086e-2` mass
  units of vapor (0.44% difference), conserved total mass, cooled the liquid,
  and kept vapor production away from the buried liquid--solid boundary.

### 2026-08-05: split-volume consistency and end-to-end validation

- Corrected the pressure source so the mixture-volume defect removes volume
  changes already represented by split chemistry, diffusion, and phase-change
  operators before the exact split dilatations are added.  The previous form
  counted the same local phase-change expansion twice.
- Replaced the abrupt hot-gas aluminum validation with gradual surface
  heating.  The abrupt initial temperature jump had projected mixed
  solid--gas cells directly to the melting temperature, causing a nonphysical
  startup expansion and a timestep near `1e-9 s` before the actual heating
  problem began.
- Made the long aluminum validation open at the domain boundary.  A closed,
  constant-pressure periodic box cannot accommodate the large specific-volume
  increase from dense liquid aluminum to vapor without driving a domain-scale
  velocity and an ordinary advective CFL restriction.
- The 0.5 ms validation heats an initially solid aluminum particle from 900 K,
  melts more than 99% of it, retains a liquid core, forms a measurable exterior
  vapor layer, excludes vapor from the buried solid--liquid interface, and
  conserves aluminum after accounting for material leaving the open domain.

### 2026-08-05: diffuse-surface kinetics and final validation

- A 1.5 ms AP/HTPB surface regression initially narrowed from 22.1 micrometers
  to 6.3 micrometers.  The cause was not a timestep instability: evaluating
  the strongly temperature-dependent Arrhenius flux independently in every
  mixed cell made the hot edge of one diffuse surface regress faster than its
  cold edge.
- A first correction reconstructed temperature at an `eta_0=0.5` contour.
  Although local, that construction imposed an unnecessary sharp-interface
  interpretation on diffuse partial densities and was removed.
- The final 0.5 ms three-phase aluminum run completed in 501 steps without a
  phase-change timestep collapse.  Its timestep remained between `1.0e-7 s`
  and `1.86e-6 s`; all solid melted, 99.99% of the retained aluminum was
  liquid, and the maximum exterior vapor fraction was `1.26e-4`.
- All modified production inputs were parsed with every input consumed.  The
  final C++ unit suite and long `LMRFSandwich`, `LMThreePhaseAluminum`, and
  coarse/fine `LMHertzKnudsenPhaseChange` checks passed.
- Kinetic temperature is now the normalized diffuse moment
  `K_ell*(delta_0g T)/K_ell*delta_0g`.  `K_ell` integrates along the local
  diffuse normal over `[-atanh(0.8)*ell/2,atanh(0.8)*ell/2]`, the 0.1--0.9
  band for `eta=0.5[1-tanh(2x/ell)]`, with
  `ell=interface.thickness`.  Five-point Gauss--Legendre quadrature evaluates
  that physical interval using multilinear interpolation of the diffuse
  fields.  Its physical sample locations therefore scale with `ell` rather
  than a fixed number of cells.  LowMach allocates the required per-level
  ghost extent from `ell/dx` and fills the moment field across coarse/fine
  boundaries.  The implicit
  latent/coupled-enthalpy response uses the matching
  `K_ell*(delta_0g^2/C)/K_ell*delta_0g` moment.  Both are normalized sums of
  resolved diffuse cells; no fixed cell stencil, contour, signed distance,
  extrapolation, level set, or sharp surface is constructed.
- Hertz--Knudsen vapor partial pressure remains local to the gas
  thermodynamic state.  Averaging that feedback allowed exponentially small
  diffuse gas tails at a buried liquid--solid boundary to act as vapor
  nucleation sites; local partial pressure suppresses that nonphysical growth
  while the temperature driving the surface flux remains diffuse-averaged.
- With the final moment formulation, the compact three-phase regression melts
  all solid, develops a 2374.9 kg/m^3 liquid core, retains 99.99% liquid, and
  produces an appreciable `1.26e-4` exterior vapor fraction.  Coarse/fine
  Hertz--Knudsen vapor masses are `1.64359e-2` and `1.65086e-2`, with both
  buried-interface checks passing.
- A 100-step AP/HTPB timing is 2.23 s for the physical-normal formulation,
  compared with 2.00 s for the earlier local kinetic implementation and
  1.75 s for the branch baseline.  Skipping thermochemical reconstruction in
  zero-interface cells and bypassing the general AMR fill copy on a
  single-level hierarchy keep the added physical averaging bounded, although
  the diffuse kinetic model still has measurable overhead.
- The exact final configuration completed the full 1.5 ms sandwich in 3745
  steps and 79.27 s.  Its extreme-front and lateral-mean rates were 11.79 and
  8.20 mm/s, respectively, with a self-sustaining 2603 K GrossModel flame.
- `interface.thickness` is now parsed independently of capillary enablement
  and is required for every kinetic phase-change mechanism.  All supplied
  kinetic inputs state the physical thickness used by their initial diffuse
  profile; a missing value is a configuration error rather than an implicit
  grid-dependent default.
- An explicit contour/extrapolation audit found no remaining reconstruction in
  LowMach phase change or liquid interfacial thermodynamics.  The separate
  Eulerian-solid reference-map algorithm still extrapolates the deformation
  map from a configurable `eta_core` (default 0.5) into the low-solid-fraction
  band.  That operation conditions an otherwise undefined material map; it is
  not used to locate a liquid interface or evaluate phase-change kinetics.

### 2026-08-05: branch performance and consistency check

- A clean `origin/eulerian-solids` executable was built separately and run on
  the identical gas-only `LMDrivenCavity` problem for three 50-step samples.
  Current runtimes were 0.59, 0.65, and 0.62 s; baseline runtimes were 0.59,
  0.64, and 0.59 s.  The current mean was 2.2% higher, with fully overlapping
  run-to-run ranges, so the new inactive liquid/interfacial paths do not cause
  a significant gas-only wall-time regression.
- Gas-only results agree to roundoff: relative L2 differences were
  `1.31e-16` and `2.62e-16` for the two velocity components, zero for pressure
  and temperature, and `4.64e-17` for density.
- A branch-native 100-step `LMRFSandwich` sample took 2.00 s on the current
  branch and 1.75 s on the baseline.  Those cases do not have identical input
  physics: the current dimensional regression model advanced 31% farther in
  simulated time.  Its cost per simulated time was therefore about 13% lower;
  the raw 14% per-step difference is not evidence of added stiffness.

### 2026-08-07: flame-temperature transport and chemistry subcycling

- Retained the calibrated constant-heat-capacity GrossModel unchanged.  A
  proposed temperature-dependent heat capacity was removed rather than using
  thermodynamic retuning to suppress the aluminum-case ignition transient.
- Gas temperature transport now advects volumetric sensible enthalpy with the
  same conservative operator as the gas partial densities and converts that
  balance back to temperature.  The periodic thermal-contact regression
  conserves sensible enthalpy to its `2e-5` relative tolerance and avoids the
  heat created by independently mixing temperature and density.
- In the production aluminum case, the pre-ignition timestep of
  `1.683608338e-7 s` equals the explicit capillary estimate to the printed
  precision.  That limiter uses the configured maximum interfacial stiffness,
  minimum reference density, and finest-cell spacing; it does not represent a
  chemical timescale.
- Fresh, matched runs through 0.321 ms compared one and four implicit
  Backward-Euler chemistry substeps per Strang half-step.  Four substeps raised
  the 0.31 ms maximum temperature from 5578.8 K to 5752.5 K (`+3.11%`), raised
  the 99.9th percentile from 5449.1 K to 5575.4 K (`+2.32%`), and raised the
  maximum speed from 49.83 to 53.60 m/s (`+7.57%`).  The minimum global step
  fell from `2.86e-8 s` to `2.62e-8 s`.  Fixed local chemistry subcycling was
  therefore rejected as a flame-spike mitigation: it removes some
  backward-Euler damping and resolves a sharper reaction/expansion impulse,
  while leaving transport and pressure projection at the global-step cadence.
- Added an optional source-based chemistry timestep informer.  It reports the
  minimum timescale for a fractional temperature rise or major-reactant
  depletion, using `chemistry.timestep.reactant_mass_fraction_floor` to keep
  trace reactants from setting the global step.  `chemistry.timestep.mode` is
  `off`, `report`, or `limit`; only `limit` participates in dynamic timestep
  selection.  The corresponding allowed change is configured with
  `chemistry.timestep.max_fractional_change`.  This is explicitly an accuracy
  indicator for split implicit chemistry, not an explicit stability CFL.
- Aluminum-flame restart samples show that this is primarily a flame-growth
  restriction rather than a permanent chemical CFL.  With a 0.1 fractional
  change, the candidate was `8.59e-11 s` at 0.310 ms and `2.53e-10 s` at
  0.350 ms, versus accepted hydrodynamic steps of `7.36e-8 s` and `6.28e-8 s`.
  At 6.300 ms, the candidate had relaxed to `2.98e-7 s`, above the accepted
  `1.0e-7 s` step.  A hard chemical limit would therefore dominate initial
  flame development, but it need not dominate a mature flame indefinitely.
- LowMach dimensional timestep bounds, refinement criteria, gravity, and
  initial/boundary field values now use the unit-aware input path.  Bare SI
  values remain valid.  CFL values, phase fractions, and the two normalized
  chemistry timestep controls remain intentionally nondimensional.

# Current model

## Conserved state and phase reconstruction

Every gas, liquid, and solid material is stored as a partial density
\(\widetilde\rho_k\).  Partial densities are the only authoritative material
state.  The mechanical volume fractions used for plots and interfacial forces
are reconstructed from them:

$$
\eta_k=\operatorname{clamp}(\widetilde\rho_k/\rho_k^\ast,0,1)
\quad\hbox{for condensed species},
$$

$$
\eta_g=1-\sum_{k\in\mathrm{condensed}}\eta_k.
$$

Overfilled condensed fractions are normalized only in the derived mechanical
state; the conserved masses are not silently modified.  The output names
`liquid_species_eta_<species>`, `rigid_species_eta_<species>`, and `gas_eta`
all follow this same reconstruction.  Gas species share one aggregate
mechanical gas phase, while retaining their individual thermochemical partial
densities.

All phases use the common LowMach velocity, pressure, and temperature.  The
model is therefore a one-fluid diffuse-interface model, not a particulate
slip or unresolved-film model.

## Interfacial free energy

For a fluid pair \(i,j\), where each pair is liquid--liquid or liquid--gas,

$$
F_{ij}=\sigma_{ij}\int\left[
\frac{3\ell}{4}|\eta_i\nabla\eta_j-\eta_j\nabla\eta_i|^2
+\frac{12}{\ell}\eta_i^2\eta_j^2\right]dV.
$$

The coefficients are calibrated so an isolated binary interface has energy
\(\sigma_{ij}\) per area.  These are the only coefficients named
`surface_tension`.

A liquid--solid boundary needs both a positive profile regularization and a
signed physical preference relative to exposed solid.  Its energy is

$$
F_{LS}=\kappa_{LS}F_0
+2(\Delta\gamma_{LS}-\kappa_{LS})
\int h(\eta_L)\delta_\ell(\eta_S)dV,
$$

where \(F_0\) is the normalized pair functional, \(\kappa_{LS}>0\),
\(h(q)=q^2(3-2q)\), and

$$
\delta_\ell(\eta_S)=
\sqrt{|\nabla\eta_S|^2+r^2}-r,
\qquad r=10^{-12}/\ell.
$$

Because \(2\int_0^1h(q)dq=1\), the energy of a flat binary liquid--solid
interface is exactly

$$
\Delta\gamma_{LS}=\gamma_{SL}-\gamma_{SG}.
$$

The positive \(\kappa_{LS}\) maintains a bounded diffuse interface even when
the referenced physical surface-energy difference is negative.  It is a
numerical regularization with physical units, not a second independently
measured interfacial energy.  The normalization removes the former
factor-of-two ambiguity.

The user supplies either `surface_energy_difference` or a contact angle
measured through the liquid, never both.  Young's equation gives

$$
\Delta\gamma_{LS}=-\sigma_{LG}\cos\theta.
$$

Solid--gas energy is the reference zero, and no solid--solid or solid--gas
pair force is added.  All pair and surface terms are reduction-consistent:
terms involving an absent phase vanish.

## Chemical potential and capillary projection

The code evaluates the variational derivatives of the full free energy and
constructs

$$
\mathbf f_{\mathrm{cap}}=\sum_i\mu_i\nabla\eta_i.
$$

`capillary_acceleration` is evaluated on faces with the same inverse density
and gradient stencil used by the variable-coefficient pressure projection.
This balanced-force placement lets pressure cancel a static capillary load and
avoids the ambiguous former name `acceleration`.

Capillarity remains explicit.  The capillary-wave estimate is

$$
\Delta t_\sigma=C_\sigma
\sqrt{\rho_{\min}\Delta x^3/\Gamma_{\max}},
$$

where \(\Gamma_{\max}\) conservatively includes fluid surface tensions and
liquid--solid surface-correction stiffness.  An implicit capillary method
would require a coupled nonlinear phase/velocity/pressure solve; the existing
scalar implicit diffusion operators cannot provide that coupling.  No nominal
or partially implicit option is retained without such a solver.

## Interface transport

The default `interface.transport.type=advected` transports every partial
density conservatively with the normal LowMach advection operator.  The
bounded Superbee limiter is the accepted long-transport configuration.

`conservative_allen_cahn` is an optional conservative profile-restoration
operator.  It acts on partial densities after transport and preserves material
mass, but it is explicit and useful only at a sufficiently small configured
diffusivity.  It is not used by default because the accepted stable setting
did not improve the 500-step profile error over Superbee advection.

An explicit Cahn--Hilliard candidate was implemented and rejected.  Its
fourth-order timestep restriction destabilized the pressure solve before a
mobility large enough to improve the interface became useful.  A future
Cahn--Hilliard option should use a genuinely coupled implicit mixed
phase/chemical-potential solve.

No part of the current model uses a level set, signed-distance field,
reinitialization, geometric VOF reconstruction, or sharp-interface jump.
Topology changes occur naturally through the diffuse partial-density fields.

## Phase change

`PhaseChange` transfers conserved mass from one condensed `phase0` material to
one condensed material or a fixed-composition gas bundle.  It does not own the
interface mechanics.  Each mechanism explicitly selects one closure:

- `equilibrium_enthalpy`: availability-constrained local enthalpy projection
  at `transition_temperature`; used for rapid solid--liquid equilibrium.
- `arrhenius_surface_flux`: dimensional surface regression
  \(j=j_\mathrm{ref}(p/p_\mathrm{ref})^n
  \exp[-T_a(1/T-1/T_\mathrm{ref})]\); used for AP/HTPB surfaces.
- `hertz_knudsen`: liquid evaporation/condensation driven by the difference
  between Clausius--Clapeyron saturation pressure and local vapor partial
  pressure.

Kinetic fluxes are converted to volumetric sources using

$$
\delta_{0g}=|\eta_g\nabla\eta_0-\eta_0\nabla\eta_g|.
$$

This integrates to unit area across a binary exposed surface and vanishes at
a buried condensed--condensed interface.  It therefore evaporates the outside
of a molten shell without creating vapor along its solid core.

Arrhenius and Hertz--Knudsen rates use a diffuse interfacial temperature

$$
T_\Gamma=\frac{K_\ell*(\delta_{0g}T)}
                 {K_\ell*\delta_{0g}},
$$

Here, \(K_\ell\) is five-point Gauss--Legendre integration along the local
diffuse normal over the physical 0.1--0.9 interface band.  The quadrature
locations scale with `interface.thickness`; multilinear interpolation makes
them independent of a particular cell stencil.  Averaging temperature before
evaluating the exponential kinetic law prevents a hot diffuse tail from
dominating the rate.  Because the moment is normalized by the same positive
surface measure, it remains bounded by the sampled diffuse temperatures.

The locally implicit thermal response uses

$$
R_T=\frac{K_\ell*(\delta_{0g}^2/C_V)}
           {K_\ell*\delta_{0g}},
\qquad
T_\Gamma(j)=T_\Gamma^n-j\,\Delta t\,\Delta h\,R_T.
$$

The same diffuse measure therefore defines localization, temperature, and
latent/coupled-enthalpy feedback.  No `eta` contour, signed distance, or
surface position enters the governing calculation.  Vapor partial pressure
is evaluated from the local gas partial density and gas volume fraction; it
is a bulk gas thermodynamic state rather than a reconstructed surface value.

The Hertz--Knudsen closure locally iterates flux, vapor pressure, latent
cooling, coupled heat, and temperature.  Mass transfer is bounded by available
material.  This local implicit coupling removes a numerical phase-field rate
restriction, although the physical gas expansion can still produce velocity
and reduce the advective timestep.

`latent_heat` is positive for endothermic forward conversion.  The optional
signed `coupled_enthalpy_change` represents an unresolved process tied only to
forward conversion: negative is exothermic and positive is endothermic.  In
the current aluminum demonstration, inert aluminum vapor is retained while a
negative coupled enthalpy approximates a selected fraction of unresolved
oxidation heat.  This shortcut is explicitly not aluminum combustion
chemistry.

## Time advancement and volume constraint

A global step performs advection and explicit sources, implicit conduction
and viscosity when enabled, optional conservative interface restoration,
equilibrium phase projections, and then kinetic surface phase changes.  The
updated partial densities are projected back onto the equation of state and
mixture-volume constraint before the next step.

Split chemistry, diffusion, and phase change report their exact integrated
volume changes.  Since the state already contains those changes, the residual
mixture-volume defect first subtracts their contributions; the projection then
adds each exact split dilatation once.  This avoids treating physical
solid--liquid or liquid--gas expansion twice.

For reactive cases, `chemistry.timestep.mode=report` evaluates a local
source-based chemical timescale without changing the step.  `limit` multiplies
that timescale by `chemistry.timestep.max_fractional_change` and includes the
result in dynamic timestep selection.  The normalized rate is the larger of
the fractional temperature-rise rate and the destruction rate of each
reactant above `chemistry.timestep.reactant_mass_fraction_floor`.  Since the
chemistry solve is implicit, this is an optional splitting-accuracy control;
it is not required for nonlinear stability.

Both `chemistry.timestep.max_fractional_change` and
`chemistry.timestep.reactant_mass_fraction_floor` are nondimensional fractions.
All associated times reported in diagnostics, as well as
`dynamictimestep.min` and `dynamictimestep.max`, are in seconds internally and
accept unit-bearing time inputs.

## Input contract

A minimal interfacial configuration is:

```text
interface.enabled = 1
interface.thickness = 8.0_um
interface.surface_tension.AlLiquid__gas = 0.85_N/m
interface.liquid_solid.enabled = 1
interface.liquid_solid.regularization.AlLiquid__AlSolid = 0.85_N/m
interface.liquid_solid.contact_angle.AlLiquid__AlSolid = 30.0_deg
```

The contact-angle line may instead be
`interface.liquid_solid.surface_energy_difference.AlLiquid__AlSolid` with
units of N/m.  Supplying both, omitting both for an active pair, using a
nonpositive regularization, or defining `surface_tension` for a solid pair is
an input error.  Legacy `liquid.*` and dimensionless phase-change rate inputs
are intentionally unsupported.  `interface.thickness` is also required when
kinetic phase change is active with capillary forces disabled, because it
defines the physical normal-integration interval rather than only a capillary
coefficient.

## Validation and acceptance criteria

The model is covered at three levels:

- Unit tests check pair-energy normalization, the signed liquid--solid surface
  correction, Young's-equation conversion, zero bulk force, flat-interface
  force balance, and reduction consistency.
- `LMLiquidCapillary` checks static and advected droplets, liquid--liquid and
  liquid--solid interfaces, contact angles, mass conservation, positivity,
  and 500-step transport.
- `LMEquilibriumPhaseChange`, `LMPhaseChange`, and
  `LMHertzKnudsenPhaseChange` check enthalpy conservation, coupled-heat signs,
  coarse/fine vapor production, and exterior-only evaporation.
- `LMThreePhaseAluminum` heats a cold solid particle long enough to complete
  melting and produce exterior vapor.  Its accelerated laser flux is a
  validation device, not a production motor calibration.  The production
  AP/HTPB/aluminum input uses calibrated dimensional parameters, but its full
  multi-millisecond motor run is not substituted for this compact regression
  test.

Tests must consume every input; none may use `allow_unused`.

## Scope and known limitations

- All phases share one resolved velocity.  Unresolved phase slip, oxide-shell
  fracture, and subgrid lubrication are not modeled.
- Thin liquid or oxide layers must span enough cells to resolve the diffuse
  profile; three to five cells is a practical minimum.
- Capillarity is explicit and scales approximately as \(\Delta x^{3/2}\).
- Pairwise interface work scales as \(O(N^2)\) in the number of mechanical
  phases, while phase reconstruction and transport scale as \(O(N)\).
- The saturation-pressure closure is pure-component Clausius--Clapeyron.  It
  does not yet include aluminum chemistry, dissolved oxide, curvature
  correction, or a chamber-specific nonideal vapor-pressure model.
- Strong evaporation in a low-Mach constant-pressure domain needs an open
  outflow or sufficient expansion volume.  A closed periodic box is not a
  physically meaningful long-time evaporation experiment.

Primary implementation locations are
`src/Model/PhaseField/MultiphaseInterface.H`,
`src/Model/Mechanism/PhaseChange.H`, `src/Integrator/LowMach.cpp`, and
`src/Operator/PressurePoisson.cpp`.
