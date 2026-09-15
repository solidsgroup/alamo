# A/E-only calibration: algebraic feasibility and physical regime

This is an analysis-only assessment of the corrected user request. No
Arrhenius parameter file or simulation is created here. Binder thermal
properties are fixed at rho=920 kg/m³, cp=2418.29 J/(kg K),
Q=−1255200 J/kg (−300 cal/g), k=0.13 W/(m K), T0=300 K.
All AP properties remain fixed. Targets are the previously digitized Chen
Figure 4 pure-binder solid-curve endpoints, not experimental measurements.

The energy balance determines a required surface temperature from each
target before Arrhenius parameters enter:

    Ts = T0 + [q/(rho*r) + Q]/cp.

| q [cal/(cm² s)] | Target r [cm/s] | Required Ts [K] |
|---:|---:|---:|
| 200 | 0.801937820 | 249.968430 |
| 500 | 1.755380052 | 316.621838 |
| 1000 | 3.148772971 | 378.202748 |

Holding q=500 out, the q=200 and 1000 targets give the unique mathematical
two-point Arrhenius fit r=A exp[−(E/R)/Ts]:

    E/R = ln(r1000/r200)/(1/T200 − 1/T1000) = 1008.342770 K,
    A = r200 exp[(E/R)/T200] = 0.452931711 m/s.

The resulting full energy/kinetic closure predicts q=500 at
Ts=311.190109 K and r=1.773362187 cm/s, about +1.0244% from the held-out
Chen target. This is an analytic prediction, not a measured solver result.
The three inferred temperatures/rates are not exactly collinear in
log(r) versus 1/T: the 200/500 pair implies E/R=930.234738 K and the
500/1000 pair implies 1136.254826 K. Therefore one constant A/E pair
cannot exactly match all three with these fixed thermal properties.

The low-flux endpoint requires Ts approximately 50 K below the initial/deep
solid temperature. If the intended hot-surface pyrolysis regime requires
Ts≥300 K, its maximum energy-limited rate at q=200 is
q/[rho*(−Q)]=0.724637681 cm/s. Reaching the Chen target at Ts=300 K
would require q≥221.334838 cal/(cm² s). Changing A/E cannot remove this
energy-bound conflict while cp/Q and the prescribed net flux remain fixed.

However, a below-T0 surface is not algebraically or thermodynamically
forbidden by this energy balance: endothermic regression can draw sensible
heat from cooling incoming solid, giving a reversed solid-side temperature
gradient. That is a different regime from the original hot pyrolysis
interpretation, and the extremely low inferred activation temperature is
an effective fit rather than established physical kinetics. Whether to
accept that regime is a modeling decision; this assessment does not silently
change cp/Q, the target, or the meaning of the incident heat flux.

Two distinct positive inferred temperatures identify A and E/R only
conditional on the fixed thermal properties and surface-flux definition.
Neither the article curve fit nor its small residual establishes physical
material properties, experimental validation, or agreement of mixed
predictions. Chen graphical resolution is approximately 0.005 cm/s.

## Transfer to the original input's temperature cutoff

The standalone study omits `temperature_cutoff`, retaining the 0 K default
in `src/Model/Mechanism/PhaseChange.H`. That implementation returns a zero
kinetic factor unless temperature is strictly above the cutoff.
`input.lm.ap_htpb` instead sets
`HTPB_pyrolysis.phase_change.temperature_cutoff = 360.0_K`.
Both q=200 and q=500 targets require temperatures below that value when
the thermal data are fixed. Consequently, inserting fitted A/E into that
original input while keeping its cutoff cannot reproduce the same targets
in this steady heat-flux closure. No cutoff or production input is changed;
the fit is specific to the executed standalone configuration.
