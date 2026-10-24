# Homogeneous binder phase change

Set `<mechanism>.phase_change.homogeneous = true` to represent unresolved
particles inside a rigid binder species. The default is `false`.
The usual `in`, `out`, `rate_multiplier`, and `activation_temperature`
describe the pure binder. Its species density, specific heat, and conductivity
also remain the **pure constituent** input values. Add:

```text
HTPB_pyrolysis.phase_change.homogeneous = true
HTPB_pyrolysis.phase_change.homogeneous.ap_solid = AP_solid
HTPB_pyrolysis.phase_change.homogeneous.total_mass_fraction = 0.85
HTPB_pyrolysis.phase_change.homogeneous.resolved_mass_fraction = 0.65
HTPB_pyrolysis.phase_change.homogeneous.ap_gas = AP_gas
HTPB_pyrolysis.phase_change.homogeneous.ap_rate_multiplier = 2.75e6
HTPB_pyrolysis.phase_change.homogeneous.ap_activation_temperature = 3145.0_K
HTPB_pyrolysis.phase_change.homogeneous.binder_heat_release = -300.0_cal/g
HTPB_pyrolysis.phase_change.homogeneous.ap_heat_release = -100.0_cal/g
```

The AP species supplies the pure AP density, specific heat, and thermal
conductivity. It may also occur as resolved particles with its own ordinary
phase-change mechanism; their densities and kinetics are unchanged. An AP
species must be present in `species.names` even for a fully homogenized pack,
but it need not have a density IC or a separate mechanism.

Both total and resolved fractions are fractions of the **whole propellant's
mass**. The fraction within the binder blend is
`w = (total - resolved)/(1 - resolved)`: the 85%/65% example therefore gives
`w = 4/7`, not `0.20`. Alternatively, specify
`homogeneous.mass_fraction = 0.5714285714285714` and omit both total/resolved
keys. These inputs describe composition; they do not generate or measure the
resolved particle geometry.

The implementation follows `homogeneous_binder_model.pdf`, with the heat
capacity rule from Chen, Buckmaster, Jackson, and Massa, *Homogenization issues
and the combustion of heterogeneous solid propellants*, Proceedings of the
Combustion Institute 29 (2002), 2923–2929:

- Density adds constituent volumes:
  `rho_blend = 1 / ((1-w)/rho_B + w/rho_AP)`.
  Equivalently, density is volume-weighted. Specific heat and signed heat
  release, both expressed per unit mass, use mass-weighted arithmetic averages.
- `t = w*rho_B/(w*rho_B + (1-w)*rho_AP)` is the AP volume fraction in the blend.
- The conductivity solves equations (9)/(11) on the physical interval between
  the constituent conductivities, using the build's dimension (2D or 3D).
- The rate multiplier uses a geometric average with volume fraction `t`, and
  activation temperature `E/Ru` uses an arithmetic average with `t`. Both act
  on the existing Allen–Cahn operator, retaining its mobility and interface
  parameters. Constituent rate multipliers must use the same phase-field
  calibration; these dimensionless multipliers are not speeds in m/s.
- Transfer removes blend mass and releases binder/AP gas products in
  fractions `1-w` and `w`. `homogeneous.ap_gas` identifies the AP gas product;
  the ordinary `out` products describe the binder. Gas expansion uses each
  product's molecular weight.

LowMach installs the effective properties throughout the solid-volume and
thermal calculations. **Density initial and boundary conditions must describe
the blended material**, not pure HTPB: for the example and constituent
densities 920 and 1950 kg/m³, use `rho_blend = 1317.7334732423922 kg/m³`
times the desired blend volume-fraction profile. Restart densities must also
use this composition. AP-containing blend is represented by the one binder
species, so do not also initialize its subgrid AP as resolved `AP_solid`.

Temperature remains coupled to LowMach's transient energy equation. The
mass-weighted surface heat is applied with the actual phase transfer, and its
thermal expansion is included in the projection. Positive heat release heats
the mixture; negative heat release cools it. The notes' relation
`Ts = T0 + Q/c + q/(rho*c*r)` describes the steady planar limit under a
prescribed heat flux; it is not imposed as a second temperature equation in
this transient solver. No prescribed-flux steady closure is added.

Rocfire's gas reactions already include `chemistry.model.rocfire.qsolid`.
If those values are retained, set the corresponding homogeneous heat-release
inputs to zero. To place the blend's solid heat at phase transfer, account for
it only once: for example, separate its gas products/chemistry from resolved
AP if resolved AP still relies on Rocfire's lumped solid heat. This feature
does not change the Rocfire reaction heats automatically.

The accompanying input uses artificial kinetics and heats for a one-step
regression test. It checks effective conductivity, conserved mass, AP gas
production, and the temperature change due to phase-transfer heat.
