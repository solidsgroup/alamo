# Mixing and runtime diagnostics

## Main recommendation: use a volume-additive mixture density

With AP volume fraction t, the volume-additive density is `rho=(1-t)*rho_binder+t*rho_AP`. Equivalently, for AP mass fraction w, `1/rho=(1-w)/rho_binder+w/rho_AP`. Keep cp and Q mass-weighted. This makes `rho*cp=(1-t)*rho_binder*cp_binder+t*rho_AP*cp_AP`, with the corresponding volume-weighted energy density for Q.

The current implementation in `src/Model/Mechanism/PhaseChange.H` instead uses an arithmetic mass-weighted density. At t=0.4 it gives 1523.153 kg/m³ instead of 1332.000 kg/m³. This raises both rho*cp and rho*Q by 14.351%, altering the energy available per unit regressed volume and lowering the coupled surface temperature.

The calculation below changes **only density in the analytic closure**, retaining the frozen A/E, cp, Q, T0 and AP parameters. At q=500, t=0.4, the predicted surface temperature increases from 843.600 to 854.740 K. The analytic rate error against Chen falls from -12.163% to -1.004%.

| AP volume fraction | Current analytic error [%] | Volume-additive density analytic error [%] |
|---|---|---|
| 0.0 | +1.409 | +1.409 |
| 0.2 | -10.486 | -0.090 |
| 0.4 | -12.163 | -1.004 |
| 0.6 | -10.243 | -1.546 |
| 0.8 | -6.684 | -1.933 |
| 1.0 | -2.134 | -2.134 |

These are diagnostic analytic predictions, not new LowMach runs with a changed density rule. No production code, mixing rule or constituent parameter was changed. The correction would need a consistent implementation in phase change, thermal capacity, initial density and reference/deck generation, followed by a new sweep. The earlier article audit records a density-weighting ambiguity in Chen's text versus its plotted blend results; this diagnostic does not establish which exact properties the authors used. See `../../analysis/REFERENCE_NOTES.md`.

Pure AP is also an unfitted endpoint and is underpredicted by about 3.7–3.9% in the present solver sweep. After resolving density mixing, I recommend calibrating AP A/E separately against the pure-AP endpoints using fixed, supported AP thermal properties and the same held-out strategy. Density mixing alone cannot change a pure endpoint. Only after those checks should residual blend errors motivate a different effective kinetic interpolation; Chen's Eq. (15) supplies a target from the two pure rates, with Eq. (14) supplying its coupled temperature.

## Surface temperature is an output of the prescribed-flux problem

Chen equations (12)–(14) determine r and Ts together:

`q = rho*r*[cp*(Ts-T0)-Q]`, `r = A*exp[-(E/R)/Ts]`.

Prescribing q while solving these equations does not require prescribing Ts independently. Changing density changes that energy balance, then changes Ts and the Arrhenius rate. Thus the observed mixing error is mediated by temperature, but is not evidence that the problem needs an independently prescribed Ts. The measured LowMach/analytic discrepancy in the completed sweep is only -1.760% to -1.482%, much smaller than the largest error against Chen. Incident-flux and interface-width effects remain a separate numerical check. Conductivity controls the thermal profile and relaxation time, but it cancels from this ideal steady surface balance, so changing k to fit steady rates would obscure the issue.

Source: [Chen et al. (2002)](https://doi.org/10.1016/S1540-7489(02)80357-1), user-supplied PDF, printed pp. 2926–2927. The density alternative is our conservation-based diagnostic, not an asserted quotation of their implemented rule.

## Q units, averaging and gas chemistry

The supplied binder Q is −66 cal/g = −276144 J/kg; the retained AP value is −100 cal/g = −418400 J/kg. `PhaseChange.H` parses both with `Unit::Energy()/Unit::Mass()`, averages them as `Q=(1-w)*Q_binder+w*Q_AP`, then applies heat as `-mass_change*Q`. During regression mass_change is negative, so these negative Q values remove heat. The reference generator uses the same mass-specific values and averaging. The gas molecular weight and `system.amount=kmol` do not turn this Q into a molar quantity.

At 40% AP by volume, w=0.585585586, giving Q=-359447.063063 J/kg = -85.909910 cal/g. This Q averaging is appropriate for additive mass-specific constituent heats and should be retained. With volume-additive density, `rho*Q=(1-t)*rho_binder*Q_binder+t*rho_AP*Q_AP`; the diagnostic verifies this equality for all 18 cases. The present inconsistency is the density used to convert Q into heat per regressed volume, rather than the mass weighting of Q itself.

For a constituent Q given in J/mol, use mole/amount fractions: `Q_molar=sum(x_i*Q_molar_i)`, where `x_i=(t_i*rho_i/M_i)/sum(t_j*rho_j/M_j)`. Volume weighting is equivalent only when the constituents have equal molar density rho_i/M_i. An alternative is to convert each molar Q to J/kg using its constituent molar mass, then apply mass weighting. [IUPAC amount-fraction definition](https://goldbook.iupac.org/terms/view/A00296) and [volume-fraction definition](https://goldbook.iupac.org/terms/view/V06643) distinguish these weights. These formulas assume additive constituent heats without an extra heat of mixing.

Every executed sweep and duration-check deck uses `chemistry.model.type = frozen` and only the binder-regression phase-change mechanism. Gas-phase reactions therefore add no heat; the prescribed gas-slab heat source and the endothermic solid-to-gas Q are the configured thermal sources. The execution audits verify these chemistry settings.

## Rate fitting and runtime sufficiency

The 18 sweep cases use six thermal relaxation times, `tau=delta/r=k/(rho*cp*r²)`. Rates are deterministic planar interface-position slopes, fit over the final 40% (last 2.4 tau), and checked against independent integrated solid-volume loss. The diagnostic is convergence to a steady rate, rather than statistical sampling of randomly packed particles or experimental noise.

Across the completed sweep, maximum absolute first-half/second-half late-window speed drift is 0.02695%. Recomputing every case with the last 20%, 30%, 40% and 50% of its data changes the measured rate by at most **0.01655%**, at `q200_t1.000_e0.25_n8_p1e+08_binder_q66_sweep` with a 20% fitting window. The windows use 8–21 saved interface positions. `fit_window_sensitivity.csv` records all 72 fits. Their slope standard errors measure numerical fit scatter; they are not confidence intervals for the physical model or Chen's data.

The representative worst-error composition (t=0.4, q=500) was rerun for twelve relaxation times, with identical constituent properties, interface width, spatial resolution, timestep ceiling, heating prescription and other solver options. Duration increased from 0.0055380031 to 0.011075943 s. The domain extends six thermal lengths deeper so the final surface retains the same eight-thermal-length separation from the cold boundary; thus this is a combined duration/deep-boundary check, not a same-domain restart.

The measured rate changes from **1.07967646 to 1.07842996 cm/s**, a **-0.11545%** change. Late-window drift is -0.00889% in the six-time run and -0.07646% in the twelve-time run. Incident-flux error also changes from -0.72595% to -0.80981%. Both rate measurements use the final 40% of their respective durations. This supports duration sufficiency for diagnosing a roughly 13.5% blend discrepancy, but does not establish 0.01% accuracy in the asymptotic rate; only this one case has been tested with doubled duration. `duration_comparison.json` and `duration_execution_audit.json` retain the comparison and execution checks.

Small window/runtime sensitivity supports using the reported steady rates. It does not establish spatial convergence, remove the approximately 1.5–2% numerical bias, or imply comparable precision in a figure with about 0.005 cm/s stroke resolution.
