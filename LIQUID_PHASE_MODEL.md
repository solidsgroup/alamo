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

- Renamed the constitutive class to `MultiphaseFreeEnergy`; its current
  implementation is `src/Model/Capillarity/MultiphaseFreeEnergy.H`.
- Replaced ambiguous internal names with
  `interfacial_volume_fraction`, `interfacial_chemical_potential`, and
  `capillary_free_energy`.  The public diagnostics remain
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

### 2026-08-10: unresolved liquid-surface heating

- The branched aluminum history showed that increasing the Hertz--Knudsen
  accommodation coefficient could not initiate appreciable evaporation while
  the highly conductive melt remained near the 933.45 K melting plateau.
- Added `interfacial_heat_source`, a prescribed unresolved heat flux localized
  by the diffuse liquid--gas surface measure.  It follows a moving liquid,
  selects no sharp contour, changes no species mass, and is distinct from a
  conservative `interphase_reaction`.
- The branched demonstration now deposits 100 MW/m2 on exposed molten aluminum
  and no longer attaches oxidation heat to the amount already evaporated.
  This removes the former vapor--heat bootstrap and avoids counting the same
  unresolved oxidation energy in two closures.
- The focused regression heats only an exposed liquid--gas interface, leaves
  a buried liquid--solid boundary unheated, and conserves component mass to
  machine precision.  The existing interphase-reaction, Hertz--Knudsen, and
  complete three-phase aluminum regressions remain unchanged and pass.

### 2026-08-10: pressure-outlet consistency and backflow stability

- Traced a spurious outlet impulse to the algebraic mixture-volume correction.
  A small gas equation-of-state defect was divided by the global timestep and
  retained at a prescribed-pressure face, converting density drift into a
  strong inward pressure correction.  Fresh states now reconcile gas partial
  density once with the low-Mach pressure, temperature, composition, and
  available condensed volume.  Runtime valid cells remain conservative.
- Pressure-outlet ghost cells now contain an equation-of-state-consistent gas
  reservoir and do not extrapolate condensed material back into the domain.
  Explicitly prescribed component Dirichlet data still take precedence.
  The algebraic mixture-drift correction fades over the configured physical
  interface thickness at an open boundary; periodic and closed faces are
  unchanged.  Capillarity instead uses a zero unresolved exterior traction.
- Added the energy-stable open-boundary traction
  `0.5 rho min(u.n,0) u` wherever a pressure outlet with extrapolated velocity
  has local backflow.  It cancels the incoming kinetic-energy flux without a
  velocity threshold, a magnitude cap, or a prohibition on physical
  backflow.  Explicit velocity boundary data bypass this treatment.
- `LMPressureOutletVolume` reduces the former 443.55 m/s startup suction to
  less than `1e-12 m/s`.  `LMPressureOutletBackflow` removes about 20% of the
  kinetic energy carried by an imposed incoming oblique stream instead of
  recycling its tangential momentum.  `LMOpenBoundaryPhase` still conserves
  liquid mass while a diffuse phase leaves an open boundary, and the complete
  `LMThreePhaseAluminum` transition still passes.

### 2026-08-12: conservative capillary momentum coupling

- A matched branched-aluminum restart isolated the persistent loss of upward
  translation and increasing liquid rotation to capillarity, rather than the
  pressure outlet or a downward bulk-gas recirculation.  Over 20 microseconds,
  the former facewise `sum(mu_i grad(eta_i))` discretization reduced liquid
  vertical speed by about 0.27 m/s and increased its angular rate, while the
  identical capillary-disabled restart retained its upward translation.
- Replaced that nonconservative face body force with the finite-volume
  divergence of the symmetric Korteweg stress derived from the same pair and
  signed liquid--solid free energies.  Every interior traction is now shared
  by its neighboring cells, face tractions are averaged down before AMR
  divergence, and stress symmetry prevents capillarity from supplying an
  internal torque.
- A physical domain boundary now uses the natural zero unresolved exterior
  capillary traction.  This works for open, closed, and periodic geometries
  without the former distance-dependent capillary fade.  The energy-stable
  pressure-outlet backflow treatment is unchanged.
- Added checks for stress symmetry, discrete net force, discrete net torque,
  static-droplet translation and rotation, translating-droplet velocity, and
  capillary behavior across an AMR coarse/fine interface.
- Repeating the matched 250--270 microsecond production restart retained an
  upward liquid speed of 1.075 m/s.  The capillary-disabled reference was
  0.997 m/s, while the former force fell to 0.731 m/s.  Liquid angular rate was
  5.20e3 1/s versus 4.06e3 1/s disabled and 5.48e3 1/s formerly.  The remaining
  liquid rotation can change as surface tension reshapes the nonspherical
  melt, but the conservative symmetric stress cannot inject net mixture
  angular momentum.  Accepted steps remained approximately
  0.09--0.17 microseconds instead of collapsing from force-generated speed.
- The tensor is evaluated once per cell and then averaged to faces.  This
  avoids recomputing every phase-pair gradient on both sides of every face,
  while producing the same synchronized conservative tractions.

### 2026-08-13: resolved branch initialization and capillary audit

- Aluminum-scale static and branched-particle runs exposed an odd/even mode
  that the low-density static regression did not.  They motivated a complete
  audit of the capillary stress and pressure-projection stencils.  The
  pressure-face force and zero-mode-removal experiments performed during that
  audit were subsequently removed; the accepted conservative formulation is
  documented in the 2026-08-14 entry below.
- The branched melt layer formerly used a hard minimum of three signed circle
  distances and subtracted three core indicators independently.  Their
  derivative kinks and overlap holes supplied unresolved curvature that
  surface tension correctly but violently tried to remove.  Both branch
  inputs now form the outer and inset-core envelopes with the same smooth
  union and subtract those two envelopes.  Its 20 micrometer smoothing length
  equals `interface.thickness`, so the smallest initialized curvature is
  actually resolved by the diffuse model.
- In the quasi-steady branch case, peak particle-region vorticity at 0.94
  microseconds fell from `7.74e5` to `4.55e4 1/s`; peak local speed fell from
  5.92 to 2.23 m/s.  Through 5 microseconds the capillary case retained an
  upward liquid speed within 0.5% of the capillary-disabled reference.
- Conservative Allen--Cahn restoration now parameterizes the physical
  relaxation time and derives \(D=\ell^2/\tau\).  Its diffusion is one
  composite backward-Euler solve; only the conservative compression flux is
  lagged.  This removes the former explicit diffusion subcycling and its
  grid-dependent stability limit.  The branched production inputs retain
  `advected`, since enabling profile relaxation did not reduce their startup
  vorticity.
- Added an aluminum-scale static regression with 2375/4 kg/m3 liquid/gas
  densities, 0.85 N/m surface tension, an 80 micrometer radius, and a 20
  micrometer diffuse thickness.  It checks the Laplace jump, liquid mass,
  checkerboard amplitude, bulk translation, rotation, momentum, and
  interfacial vorticity.  Static, AMR, no-capillary, 500-step advected, and
  500-step implicit-restoration variants all pass.

### 2026-08-14: conservative momentum and compatible projection

- Capillary momentum coupling is the finite-volume divergence of one symmetric
  Korteweg stress.  Each interior or periodic face has one synchronized
  traction, so capillary impulses telescope and the symmetric constitutive
  stress supplies no internal torque.  There is no force-mean subtraction,
  velocity-mean correction, affine pressure gradient, or zero-mode filter.
- The former nodal projection constrained a divergence stencil different from
  the arithmetic face velocity used by finite-volume transport.  Conservative
  momentum advection then retained a spurious
  \(\mathbf u\nabla\cdot\mathbf u\) contribution.  The projection now
  constrains that same arithmetic face velocity with a face-centered
  variable-coefficient operator and applies its face correction back to the
  cell velocity.  Periodic boundaries remain purely periodic: synchronization
  only makes both copies of the periodic seam represent the same geometric
  face.
- The algebraic mixture-volume defect divided by the timestep was removed
  from the projection source.  Only volume changes reported by physical
  chemistry, phase-change, and diffusion operators enter the dilatation.
  A periodic domain cannot support a nonzero integrated volume source; such a
  problem needs a physical expansion boundary rather than an artificial mean
  correction.
- Phase transfer changes source and destination partial densities by equal and
  opposite amounts at each cell and uses their common one-fluid velocity, so
  it produces no separate momentum source.  Latent and coupled enthalpy terms
  are applied locally with the transferred mass.
- `include_viscosity=0` now disables viscosity in both explicit and implicit
  paths.  The explicit flux is the symmetric Newtonian stress rather than a
  component-wise scalar Laplacian.
- The driven-cavity regression completes with its established reference
  velocity and vorticity profiles.  Static aluminum capillary tests at a
  594:1 density ratio retain negligible bulk translation and angular rate,
  and phase-change/reaction tests check mass, linear momentum, angular
  momentum, and resolved thermal-energy balances.
- A fresh 100-step developed-flame branched-particle run remained vertically
  dominated without a mean-flow correction: gas volume-mean velocity changed
  from `(-0.0174, 1.7023) m/s` to `(-0.0293, 1.8831) m/s`, and final gas RMS
  velocity was `(0.2488, 1.9139) m/s`.
- The later branch-particle jets were traced to a discrete coupling regression,
  not to aluminum mass creation or the signed liquid--solid wetting term.  The
  refactor applied capillary acceleration at cell centers before reconstructing
  the projection faces; that left a density-weighted odd/even component outside
  the pressure constraint.  The same refactor also replaced the established
  centered-cell stress average with a separate exact face-normal gradient,
  changing the shortest-wave discrete operator for thin films.  The stress is
  again averaged consistently to faces, then converted with the projection's
  face mobility and added to that exact face predictor.  With the original
  178.85 ns capillary timestep restored, the branched case retained a maximum
  liquid fraction of 0.856 and 8.22 m/s through 7.25 microseconds instead of
  reaching liquid fraction 1.94 and 51.9 m/s by 6.98 microseconds.
- A subsequent long branch run completed step 647 but aborted while forming
  the chemical-potential plot field.  An exponentially small solid-gradient
  tail was nonzero while its cubed norm underflowed in the curvature
  denominator.  The liquid--solid surface delta is now the differentiable
  thickness-scaled functional
  `sqrt(|grad(eta_s)|^2+(alpha/ell)^2)-alpha/ell`, with the nondimensional
  `alpha` exposed as `interface.liquid_solid.surface_delta_regularization`.
  Its free energy, stress, and both chemical potentials use the same
  regularization.  A unit regression reproduces the former underflow-scale
  tail and requires finite stress and chemical potential.
- Matched regularized and unregularized outputs were identical through 60
  microseconds, ruling out `surface_delta_regularization` as the source of the
  renewed particle vorticity.  A capillary-disabled run remained at
  `2.66e4 1/s` peak particle-region vorticity and 1.98 m/s through 5.16
  microseconds, whereas liquid--gas tension produced `4.40e5 1/s` and 4.10
  m/s by 4.93 microseconds.  Over the same interval, the approximately
  `3.4e-6 J/m` reduction relative to the no-capillary liquid--gas free energy
  appeared as `3.2e-6 J/m` additional kinetic energy.  This identifies the
  remaining localized vorticity as resolved relaxation of the deliberately
  non-equilibrium three-lobed free surface, rather than a force source or the
  liquid--solid tail regularization.
- The unusually quiet historical branch run was not a valid reference: its
  nonlinear rigid solve reapplied the implicit Brinkman penalty on every
  iteration and therefore added iteration-count-dependent damping.  Restoring
  that behavior would suppress the physical rounding motion.  Conversely,
  adding a pressure-balancing isotropic gauge reduced local parasitic motion,
  but it changed only the pressure representative and did not reduce the
  physical shape-relaxation force.  That experimental gauge was removed from
  the accepted formulation so capillarity remains one pressure-reduced stress
  without an auxiliary scalar field.
- A phase-weighted-potential experiment using
  `-sum_i eta_i grad(mu_i)` reduced the 4.93
  microsecond dilute-liquid peak vorticity from `4.40e5` to `1.89e5 1/s`, but
  it also generated `1.74e-3 m/s` of bulk translation in the asymmetric
  periodic wetting test, so that nonconservative discretization was rejected.
  A pressure gauge must not be used as numerical damping.

### 2026-08-14: selectable capillarity and interface kinetics

- Moved the common multiphase free energy and the input-selected interface
  evolution models into `Model::Capillarity`.  Compile-time tuple dispatch
  provides `direct_surface_tension`, `conservative_allen_cahn`, and
  `singly_degenerate_cahn_hilliard` without virtual calls in cell kernels.
- Kept the established normalized pair free energy for all three options.
  An experimental multiplicative tail-degeneracy narrowed the resolved stress
  layer and raised the aluminum-density static-droplet parasitic speed from
  8.05% to 17.8% of the capillary velocity, so it was rejected.  Tail/profile
  control is instead performed by a conservative phase-field evolution law.
- Conservative Allen--Cahn uses one composite backward-Euler diffusion solve
  with a lagged conservative compression flux.  Singly-degenerate
  Cahn--Hilliard uses pair mobility `4 eta_i eta_j`; its stiff fourth-order
  linear response is factored into two positive composite Helmholtz solves,
  while mobility and nonlinear chemical potential remain lagged.
- Both profile models conserve every liquid integral.  Their constitutive
  density remap retains cell momentum before the projection, so it cannot add
  net linear or angular momentum.  The periodic 500-step Allen--Cahn test has
  `1.91e-16 m/s` normalized momentum drift, and the Cahn--Hilliard static test
  has `1.90e-20 m/s` drift while decreasing surface-plus-kinetic energy.
- Uniform-grid and AMR static liquid--gas and liquid--solid tests retain their
  Laplace response and conservative force/torque balance.  No phase-fraction
  force threshold, level set, contour projection, or mean-flow correction was
  introduced.  An open-boundary Allen--Cahn test also retains liquid mass to
  the reported eight decimal places while advection carries the interface
  away from the pressure outlet.
- A materially active Cahn--Hilliard test (`tau=0.01 s` over `8e-4 s`) exposed
  the expected small overshoot of a polynomial free energy: the liquid
  fraction reached `1.000296` even though its integral was conserved.  Both
  profile models now finish with a conservative Gibbs-simplex projection.
  It projects the local liquid vector into the solid-limited available volume,
  then restores each pre-relaxation liquid integral by a multiplicative
  liquid--gas redistribution within that liquid's existing diffuse support.
  The active uniform and AMR tests now retain liquid volume to roundoff, have
  phase sum error at `2.22e-16`, and remain in `[0,1]` without a phase cutoff,
  contour, or sharp-interface reconstruction.
- Production testing also found why the earlier unprojected Cahn--Hilliard
  state could abort the next thermal solve.  A liquid undershoot of only
  `-5.63e-4`, multiplied by aluminum's `237 W/m/K` conductivity, made the
  assembled thermal conductivity `-0.0364 W/m/K`.  Condensed heat capacity,
  conductivity, and viscosity are now evaluated from nonnegative admissible
  material content while the conserved partial density remains visible.
  The corrected branched aluminum case advances beyond the former failure at
  its unchanged `1.79e-7 s` global step; with simplex restoration its liquid,
  gas, and five solid fractions sum to one within `2.22e-16`.

# Current model

## Conserved state and phase reconstruction

Every gas, liquid, and solid material is stored as a partial density
\(\widetilde\rho_k\).  Partial densities are the only authoritative material
state.  The mechanical volume fractions used for plots and interfacial forces
are reconstructed from them:

$$
\eta_k=\widetilde\rho_k/\rho_k^\ast
\quad\hbox{for condensed species},
$$

$$
\eta_g=1-\sum_{k\in\mathrm{condensed}}\eta_k.
$$

With direct surface tension, the diagnostic liquid and solid fractions are
not independently clamped or normalized:
`liquid_species_eta_<species>` and
`rigid_species_eta_<species>` remain direct reconstructions of the conserved
partial densities, so transport errors stay visible.  The private
constitutive phase vector used by the multiphase free energy is projected onto
the Gibbs simplex by clipping negative material fractions and, only when their
sum exceeds one, scaling all condensed fractions by the same factor.  Its gas
entry is the nonnegative complement and is reported as `gas_eta`.  This keeps
the free-energy model inside its admissible phase space without altering any
conserved partial density or hiding the raw liquid/solid diagnostics.  Gas
species share one aggregate mechanical gas phase while retaining their
individual thermochemical partial densities.

When either optional profile-relaxation model is selected, its update is
followed by a conservative Gibbs-simplex projection.  This changes the local
partial-density distribution, as any phase-field relaxation must, but exactly
restores every liquid integral recorded before relaxation.  Its redistribution
weight is `eta_liquid eta_gas`, so a correction remains on the existing
diffuse liquid--gas support and cannot seed liquid into a pure-gas region.
This admissibility operation is unrelated to the constitutive-only projection
above: the former makes the evolved phase state physical and conservative;
the latter only supplies a safe phase vector to local free-energy kernels.

Fluid, liquid, and deformable phases use the common LowMach velocity, pressure,
and temperature.  A freely moving rigid species is transported by the affine
translation and rotation obtained from its own conserved mass and momentum,
and exchanges momentum with the one-fluid field through the implicit Brinkman
constraint.  No level set or independent sharp-interface velocity is used.

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
\delta_{\ell,\alpha}(\eta_S)=
\sqrt{|\nabla\eta_S|^2+(\alpha/\ell)^2}-\alpha/\ell.
$$

Because \(2\int_0^1h(q)dq=1\), the energy of a flat binary liquid--solid
interface approaches

$$
\Delta\gamma_{LS}=\gamma_{SL}-\gamma_{SG}.
$$

Here `interface.liquid_solid.surface_delta_regularization` supplies the
positive nondimensional \(\alpha\) (default `1.0e-12`).  It makes the diffuse
surface-area functional differentiable when \(\nabla\eta_S=0\), while the
subtraction leaves its energy exactly zero in uniform bulk phases.  The
resulting \(O(\alpha)\) change to the integrated interfacial energy is
negligible at the default value.  The same regularized functional is varied
for both chemical potentials and used in the capillary stress; it is not a
force cutoff or a phase-fraction threshold.

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

## Chemical potential and capillary momentum coupling

The chemical potentials are the variational derivatives of the full free
energy.  Their fluid-interface terms provide capillary momentum coupling and
the full values are available as diagnostics:

$$
\mu_i=\frac{\delta F}{\delta\eta_i}.
$$

For the pair free-energy density \(f_{ij}\), define
\(\mathbf A_{ij}=\eta_i\nabla\eta_j-\eta_j\nabla\eta_i\).  The momentum
equation uses the pressure-reduced Korteweg stress

$$
\mathbf T_{ij}=-\frac{3\ell c_{ij}}{2}
\mathbf A_{ij}\otimes\mathbf A_{ij}.
$$

The omitted isotropic part of the canonical stress is a pressure gauge.  It is
not evaluated as a second capillary force or added to the reported mechanical
pressure.  The phase dependence is already contained in \(\mathbf A_{ij}\),
which vanishes when either member of a pair is absent and decays smoothly in a
diffuse tail without a volume-fraction threshold.

For the signed liquid-covered-solid correction
\(C_{LS}=\Delta\gamma_{LS}-\kappa_{LS}\),

$$
\mathbf T^{\mathrm{surface}}_{LS}=2C_{LS}h(\eta_L)
\left[\left(s_\alpha-\frac{\alpha}{\ell}\right)\mathbf I-
\frac{\nabla\eta_S\otimes\nabla\eta_S}{s_\alpha}\right],
\qquad
s_\alpha=\sqrt{|\nabla\eta_S|^2+(\alpha/\ell)^2}.
$$

Phase gradients and stresses are evaluated with the same
centered cell stencil used by the established diffuse-interface operator,
then neighboring stresses are averaged to their shared face.

The shared traction is synchronized across periodic/patch seams and averaged
down at AMR interfaces.  Its conservative divergence is converted with the
projection face mobility and added to the exact face predictor that pressure
constrains.  Interior stress impulses telescope exactly, and the constitutive
stress is symmetric.  No phase receives a separate capillary body force.

The pressure projection constrains the same arithmetic finite-volume face
velocity used by transport.  This avoids the conservative-form error produced
when projection and advection use different divergence stencils.  It does not
alter the mean velocity or impose an affine pressure component.  Rigid
fixed-point iterations reuse the same predictor and do not accumulate the
capillary traction more than once per iterate.

At a physical domain boundary the unresolved exterior capillary traction is
zero.  Periodic boundary faces share exactly the same synchronized traction.
No level-set surface, outlet-distance fade, velocity cap, or hard-coded force
cutoff is used.

Capillarity remains explicit, so the dynamic timestep retains the established
resolved-wave estimate

$$
\Delta t_\sigma \le C_{mathrm{CFL}}
\sqrt{\frac{\rho_{\mathrm{ref,min}}\,\Delta x^3}
{c_{\max}}}.
$$

Here, \(c_{\max}\) conservatively includes fluid surface tensions and the
magnitude of liquid--solid surface-correction stiffness, and
\(C_{\mathrm{CFL}}\) is the normal `cfl` input.  The projection removes the
pressure-like part of the capillary impulse on those same faces; it does not
make physical, shape-changing capillary motion implicit.

An implicit momentum-capillary method would require a coupled nonlinear
phase/velocity/pressure solve; the existing scalar implicit diffusion
operators cannot provide that coupling.  The optional phase-field kinetics
below treat their stiff interface-profile response implicitly, but do not make
the physical capillary-wave force implicit.

## Interface transport

`interface.model.type=direct_surface_tension` is the default and least
expensive option.  It transports partial densities with the normal conservative
LowMach advection operator and applies the common Korteweg stress, with no
additional profile kinetics.

`conservative_allen_cahn` adds the locally conservative profile flux

$$
\partial_t\eta=\nabla\cdot\left[D\left(\nabla\eta-
\frac{4\eta(1-\eta)}{\ell}\mathbf n\right)\right],
\qquad D=\ell^2/\tau.
$$

Its diffusion is one composite backward-Euler solve and its conservative
face-compression flux is lagged.  It restores the declared tanh thickness
without a contour or signed-distance reconstruction.

`singly_degenerate_cahn_hilliard` instead advances each liquid with

$$
\partial_t\eta_i=\nabla\cdot\sum_{j\ne i}
M_0\,4\eta_i\eta_j\nabla(\mu_i-\mu_j),
\qquad M_0=\frac{\ell^3}{\sigma_{\rm ref}\tau}.
$$

The single mobility zeros suppress diffusion in either pure phase while
retaining finite interfacial mobility.  Liquid--liquid pair fluxes are equal
and opposite, physical boundaries use zero chemical flux, and AMR faces share
one averaged flux.  A linear energy stabilization factors the stiff
fourth-order gradient response into two positive composite Helmholtz solves.
The nonlinear chemical potential and degenerate mobility are lagged, so this
is a linearly implicit energy-stabilized update rather than a fully coupled
nonlinear Cahn--Hilliard solve.

Both phase-field options update the conserved liquid partial densities,
preserve each liquid integral, and retain local mixture momentum through the
constitutive remap before projection.  They introduce no explicit
fourth-order stability limit, although `relaxation_time` still controls the
physical/numerical rate and temporal accuracy.

Because the common polynomial pair energy has finite derivatives at the
simplex boundary, neither linearly implicit profile update is mathematically
bound-preserving by itself.  Their shared conservative Gibbs projection
enforces nonnegative liquid and gas volume after relaxation.  It first uses
the volume not occupied by solids as the local simplex capacity, then restores
each liquid's pre-relaxation integral through the diffuse weight
`eta_liquid eta_gas`.  No user tolerance, phase-presence threshold, contour,
or global mean-flow correction enters this operation.

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
cases that use it, inert aluminum vapor can be retained while a negative
coupled enthalpy approximates a selected fraction of unresolved oxidation
heat.  This shortcut is explicitly not aluminum combustion chemistry.

## Prescribed unresolved interfacial heating

`interfacial_heat_source` supplies a positive prescribed `heat_flux` on the
exposed diffuse boundary of one named liquid:

$$
\dot q'''=\dot q'' A(T)
|\eta_g\nabla\eta_L-\eta_L\nabla\eta_g|.
$$

The bounded activation

$$
A(T)=\frac{1}{2}\left[1+\tanh\left(
\frac{T-T_a}{\Delta T_a}\right)\right]
$$

uses dimensional `activation_temperature` and `activation_width`.  It lets the
hot gas-side portion of the diffuse boundary initiate unresolved heating while
making the source negligible on a cold melt, without imposing a discontinuous
temperature cutoff.

Its volume integral converges to heat flux times liquid--gas area.  It vanishes
in bulk phases and at a buried liquid--solid boundary, follows the advected
partial-density field, and requires neither a level set nor an `eta` contour.
The mechanism changes no species.  It is therefore appropriate only for an
energy-producing process deliberately omitted from the resolved chemistry;
when reactants and products are represented, `interphase_reaction` provides
the conservative model instead.  Positive `heat_flux` adds thermal energy,
following the usual source-flux convention rather than the signed enthalpy-
change convention used by `coupled_enthalpy_change`.

The activation is smooth, bounded, and saturates at `heat_flux`, so it creates
no unbounded Arrhenius feedback.  It is integrated by the Runge--Kutta source
update, while thermal diffusion remains implicit when configured.  Its gas
thermal expansion is included in the low-Mach projection source.

## Time advancement and volume constraint

A global step performs advection and explicit sources, implicit conduction
and viscosity when enabled, optional conservative interface restoration,
equilibrium phase projections, and then kinetic surface phase changes.  The
partial densities and mixture momentum remain conservative Runge--Kutta state;
the pressure solve constrains the transport face velocity rather than
algebraically changing material mass.

Chemistry, diffusion, and phase change report their physical integrated volume
changes.  The projection includes each corresponding dilatation once.  It does
not divide a residual equation-of-state or mixture-volume defect by the
timestep.  In a fully periodic domain the volume integral of the prescribed
dilatation must be zero.  Net thermal or phase expansion therefore requires an
open physical boundary; the solver does not alter the periodic mean velocity
or add a compensating pressure gradient.

Every local phase transfer removes mass from its source and adds exactly the
same mass to its destination.  Because phases share the one-fluid velocity,
there is no separate interphase momentum impulse.  Latent heat and optional
coupled enthalpy are applied locally with the transfer, while conservative
sensible enthalpy is reconstructed using the updated phase masses and heat
capacities.

The explicit viscous flux is the symmetric Newtonian stress.  The optional
implicit momentum-diffusion solve currently treats the component-wise
Laplacian; for constant viscosity its omitted cross term is a pressure
gradient.  A fully coupled variable-viscosity vector solve remains a future
extension.

At a prescribed-pressure boundary, outflow transports the conservative
interior state.  Backflow receives an equation-of-state-consistent gas ghost
state, excludes unprescribed condensed inflow, and uses an energy-stable
momentum traction so extrapolated velocity cannot inject kinetic energy.

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
interface.model.type = direct_surface_tension
interface.surface_tension.Al_liquid_gas = 0.85_N/m
interface.liquid_solid.enabled = 1
interface.liquid_solid.surface_delta_regularization = 1.0e-12
interface.liquid_solid.regularization.Al_liquid_Al_solid = 0.85_N/m
interface.liquid_solid.contact_angle.Al_liquid_Al_solid = 30.0_deg
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

`interface.liquid_solid.surface_delta_regularization` is nondimensional and
must be positive when liquid--solid capillarity is enabled.  It is divided by
`interface.thickness` internally, so the differentiable surface-delta model
scales with the declared diffuse thickness rather than the mesh spacing.

The direct model needs no additional input.  Conservative Allen--Cahn is
selected with

```text
interface.model.type = conservative_allen_cahn
interface.model.conservative_allen_cahn.relaxation_time = 20.0_us
```

The variational alternative is

```text
interface.model.type = singly_degenerate_cahn_hilliard
interface.model.singly_degenerate_cahn_hilliard.relaxation_time = 20.0_us
```

Both relaxation times are dimensional.  Diffusivity or mobility is derived
from the declared interface thickness rather than tuned against the grid.

## Validation and acceptance criteria

The model is covered at three levels:

- Unit tests check pair-energy normalization, the signed liquid--solid surface
  correction, Young's-equation conversion, zero bulk force, flat-interface
  force balance, reduction consistency, Korteweg-stress symmetry, and
  discrete force/torque conservation.
- `LMLiquidCapillary` checks static and advected droplets, liquid--liquid
  interfaces, mass conservation, positivity, zero capillary
  translation/rotation, surface-plus-kinetic-energy behavior, AMR face
  consistency, and 500-step transport without loss of bulk velocity.  Its
  aluminum-scale variant also
  checks a 594:1 density ratio, physical aluminum surface tension, the Laplace
  jump, checkerboard amplitude, and parasitic interfacial vorticity.
- `LMLiquidSolidCapillary` applies strong wetting to an asymmetric liquid-coated
  freely moving solid and checks that internal capillarity excites resolved
  shape relaxation without generating bulk translation or excessive net
  rotation on both uniform and adaptive meshes.
- `LMEquilibriumPhaseChange`, `LMPhaseChange`, and
  `LMHertzKnudsenPhaseChange` check component and total mass, linear and angular
  momentum, enthalpy conservation, coupled-heat signs, coarse/fine vapor
  production, and exterior-only evaporation.  The closed periodic kinetic
  test is volume-neutral; open-boundary tests cover net phase expansion.
- `LMInterphaseReaction` checks component transfer, total mass, linear and
  angular momentum, and sensible energy for a zero-heat reaction.
- `LMInterfacialHeatSource` checks exposed-interface localization, zero heating
  at a buried liquid--solid boundary, and exact component-mass invariance.
- `LMPressureOutletVolume`, `LMPressureOutletBackflow`, and
  `LMOpenBoundaryPhase` check equation-of-state reconciliation, kinetic-energy
  stability during local backflow, and diffuse-phase transport across an open
  boundary.
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
- LowMach evolves conservative sensible enthalpy rather than a conservative
  total-energy variable.  Phase-change energy and short-time
  surface-plus-kinetic-energy balances are regression tested, but viscous
  dissipation is not returned to the thermal equation.
- Implicit viscosity is currently component-wise.  Strongly
  variable-viscosity cases that require the full cross-component Newtonian
  stress need a coupled vector diffusion operator.
- Pairwise interface work scales as \(O(N^2)\) in the number of mechanical
  phases, while phase reconstruction and transport scale as \(O(N)\).
- The saturation-pressure closure is pure-component Clausius--Clapeyron.  It
  does not yet include aluminum chemistry, dissolved oxide, curvature
  correction, or a chamber-specific nonideal vapor-pressure model.
- Strong evaporation in a low-Mach constant-pressure domain needs an open
  outflow or sufficient expansion volume.  A closed periodic box is not a
  physically meaningful long-time evaporation experiment.

Primary implementation locations are
`src/Model/Capillarity/MultiphaseFreeEnergy.H`,
`src/Model/Capillarity/Capillarity.H`,
`src/Model/Mechanism/PhaseChange.H`,
`src/Model/Mechanism/InterfacialHeatSource.H`,
`src/Integrator/LowMach.cpp`, and `src/Operator/PressurePoisson.cpp`.
