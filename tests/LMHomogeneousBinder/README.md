# Homogeneous AP/binder decomposition

Set `<mechanism>.phase_change.homogeneous = true` to represent subgrid AP
inside a rigid binder species. The default is false. Both constituent species
must supply their **pure** density, specific heat, and conductivity. The AP
species may also describe resolved particles; its properties are unchanged.

```text
mechanisms.names = binder_regression
binder_regression.type = phase_change
binder_regression.phase_change.phase0 = binder_solid
binder_regression.phase_change.phase1 = binder_gas
binder_regression.phase_change.kinetics = arrhenius_surface_flux
binder_regression.phase_change.reference_pressure = 1_bar
binder_regression.phase_change.pressure_exponent = 0
binder_regression.phase_change.reference_temperature = 1000_K
binder_regression.phase_change.homogeneous = true
binder_regression.phase_change.homogeneous.ap_solid = AP_solid
binder_regression.phase_change.homogeneous.total_mass_fraction = 0.85
binder_regression.phase_change.homogeneous.resolved_mass_fraction = 0.65
```

Supply the six constituent kinetic/heat parameters under `homogeneous`:

| Parameter | Meaning | Units |
| --- | --- | --- |
| `binder_pre_exponential_speed` | Binder Arrhenius prefactor A | m/s |
| `ap_pre_exponential_speed` | AP Arrhenius prefactor A | m/s |
| `binder_activation_temperature` | Binder E/R | K |
| `ap_activation_temperature` | AP E/R | K |
| `binder_heat_release` | Binder condensed decomposition Q | J/kg |
| `ap_heat_release` | AP condensed decomposition Q | J/kg |

These require a calibration appropriate to the surface-flux discretization.
The numerical values in `input` are artificial regression-test parameters.
Allen–Cahn rate multipliers from other branches are not speeds and must not
be copied into these inputs. Calibration and sweep results are kept outside
the committed tests.

Both total and resolved fractions refer to the whole propellant's mass.
The AP fraction within the blend is `w = (total - resolved)/(1 - resolved)`.
Alternatively specify `homogeneous.mass_fraction = w` and omit both other
fraction keys. For total 0.85 and resolved 0.65, `w = 4/7`.

The AP volume fraction is `v = w*rho_b / (w*rho_b + (1-w)*rho_AP)`.
The effective properties follow:

- `rho = 1 / ((1-w)/rho_b + w/rho_AP)` (additive constituent volumes).
- cp and Q are mass-weighted, since both are specified per unit mass.
- `ln(A) = (1-v)*ln(A_b) + v*ln(A_AP)` and
  `E/R = (1-v)*(E/R)_b + v*(E/R)_AP`.
- Conductivity solves Chen's unsquared physical relation
  `k-k_AP = (1-v)*(k_b-k_AP)*(k/k_b)^(1/d)`, with the build dimension d.
  The root lies between the two constituent conductivities.

All blend properties are resolved before any mechanism caches densities.
The effective density replaces the binder reference density throughout the
solver. Density initial and boundary conditions must consequently use the
**blend density times its volume-fraction profile**. Subgrid AP must not also
appear in the resolved AP density field. Restart densities must represent
the same composition.

The blend emits its entire mass into one gas species, `binder_gas` above.
This species already represents binder plus homogenized fine AP in Gross's
four-flame model. No fraction is rerouted to `AP_gas`. A separate resolved-AP
mechanism uses `phase0 = AP_solid` and `phase1 = AP_gas`. At w=1 the
homogeneous material has AP solid properties but retains the configured
lumped gas product; a reactive pure-AP validation therefore uses the separate
AP mechanism. A frozen-chemistry Chen endpoint can use either gas label.

The ordinary Arrhenius surface-flux law also accepts `pre_exponential_speed`
instead of `reference_mass_flux`, including for resolved AP:

```text
j(T,P) = rho*A*exp(-(E/R)/T)*(P/P_ref)^n
j_ref = rho*A*exp(-(E/R)/T_ref)
```

The existing diffuse surface measure supplies the volumetric source. No
factor involving mesh spacing or interface thickness is folded into A.
Temperature remains coupled to the transient energy equation.
Irreversible solid-to-gas consumption uses an upwind front gradient and a
frozen density stencil, so it can advance into initially pure solid without
leaving oscillatory remnants or racing neighboring device threads. Stefan
and recoil fluxes retain the centered phase-pair gradient. Interface width
and measured rates still require grid-convergence checks.

Negative Q is endothermic. Homogeneous setup converts it to the branch's
`coupled_enthalpy_change = -Q`, applied once during implicit mass transfer.
Do not also supply `latent_heat`, `coupled_enthalpy_change`,
`pre_exponential_speed`, `reference_mass_flux`, or `activation_temperature`
for the homogeneous mechanism: those values are derived. Gas-phase
`gross_model` reaction heats remain gas-only. Use frozen chemistry for the
prescribed-flux Chen calibration.

The checks cover both pure endpoints, a mixture, mechanism reordering,
conductivity, dimensional recession speed, mass conservation, gas-product
routing, and the sign and magnitude of decomposition heat. The unit tests
also check the implicit mass/heat update at multiple temperatures. Mixing
occurs only on the host during input setup; the device mechanism remains
the existing fixed-size, statically dispatched surface-flux model.
