# Binder cp restoration, Arrhenius recalibration and Figure 4 sweep

Binder cp is restored from **2130 to 2418.29 J/(kg K)**. Binder Q remains **−66 cal/g = −276144 J/kg**. The new solver-calibrated binder pair is **A={binder_A:.10g} m/s**, **E/R={binder_activation:.10g} K**, or **E={binder_E_kJ:.8f} kJ/mol = {binder_E_kcal:.8f} kcal/mol**. Previously A={previous_A:.10g} m/s and E/R={previous_activation:.10g} K. Only A and E/R are fitted; cp is a prescribed value.

Across 18 Figure 4 points, mean absolute error versus Chen changes from **{old_mean:.4f}% to {new_mean:.4f}%**, and maximum absolute error from **{old_max:.4f}% to {new_max:.4f}%**. The 12 unfitted interior mixtures change from **{old_mixture_mean:.4f}% to {new_mixture_mean:.4f}%** mean absolute error. The held-out pure-binder q=500 point has error **{heldout_error:+.4f}%**.

![Figure 4 comparison](figure4_comparison.png)

Each LowMach square is labeled with signed percent error `100*(r_LowMach/r_Chen−1)`. Solid curves are vector-extracted Chen Figure 4 Eq. 15 curves; dashed curves are the new analytic closure with mixed properties and a common coupled surface temperature. Faint crosses show the preceding accepted AP-calibrated sweep at binder cp=2130. Circles and asterisks retain Chen's 2D and 3D DNS points. Errors above refer to the solid curves; they are not statistical error bars.

{stats_table}

{point_table}

## Parameters and calibration

{unit_table}

The physical rate law is `r=A exp[−(E/R)/Ts]`. A is a speed prefactor; the input generator converts it to the phase-field multiplier appropriate to each interface width. `activation_temperature_K` stores E/R, not E. Conversions use R=8.31446261815324 J/(mol K) and 1 kcal=4.184 kJ. `parameter_comparison.csv` includes all constituent properties, including unchanged values and the original engineering baseline.

Binder density remains 920 kg/m³ and conductivity remains 0.213 W/(m K). All AP properties and its prior calibration history are unchanged: density=1950 kg/m³, cp=1297.9 J/(kg K), conductivity=0.4186 W/(m K), Q=−100 cal/g, A={AP_A:.10g} m/s, E/R={AP_activation:.10g} K, E={AP_E_kJ:.8f} kJ/mol={AP_E_kcal:.8f} kcal/mol. T0 remains 300 K. Three fresh pure-AP controls reproduce the previous rates within relative tolerance 1e−9.

The established two-endpoint calibration uses pure binder at q=200 and 1000 cal/(cm² s). For each target rate, the fixed balance gives `Ts=T0+[q/(rho*r)+Q]/cp`; then a two-point solve gives ln(A) and E/R. The initial analytic pair is A={direct_A:.10g} m/s and E/R={direct_activation:.10g} K. After running those two cases, each target is divided by its measured solver/analytic rate ratio and the Arrhenius inversion is repeated. Fresh endpoint runs verify that corrected pair. `freeze_record.json` records acceptance before the q=500 holdout and mixture runs were generated. The two accepted calibration cases supply two of the 18 sweep points; the other 16 are new frozen-parameter runs.

The q=500 binder point and all mixtures are excluded from fitting. The prior AP calibration is retained. These effective kinetics are conditional on the thermal properties, discretization and heating surrogate; they are not independently measured material constants. The abandoned cp-and-Q restoration is excluded from this study. This cp-only study has surface temperatures ranging from {min_Ts:.3f} to {max_Ts:.3f} K across the 18 simulations, all above T0.

Matching the pure endpoints does not fix the interior mixture response. The cp change alters the relation between heat input, surface temperature and rate. After refitting A and E, volume mixing of ln(A) and E/R with a single coupled surface temperature need not reproduce Chen's Eq. 15 harmonic interpolation of pure-component rates at their respective surface temperatures. The mixture errors therefore assess this combined mixing closure as well as the remaining numerical bias.

Density is volume-additive. Specific heat and Q, both specified per unit mass, are mass-weighted using mass fractions derived from the pure constituent densities. This conserves `rho*cp` and `rho*Q` under volume mixing. ln(A) and E/R retain volume weighting; conductivity retains the existing Chen two-dimensional rule. Gas chemistry is frozen, so there is no gas-reaction heat. The prescribed gas-slab heat input and endothermic phase-change Q remain the thermal sources. Surface temperature is solved through thermal/kinetic coupling.

## Runtime and numerical checks

Every standard sweep case reaches six thermal relaxation times, with ell/delta=0.25, eight cells per ell, 32 transverse cells and dt scale=0.2. Pressure=100 MPa, gas conductivity=100 W/(m K), gas cp=1000 J/(kg K), and temperature advection is disabled, as in the preceding accepted sweep. Rates are slopes fitted to raw eta=0.5 front positions over the final 40% of each run. Maximum late speed drift is {max_drift:.5f}%. Using the final 20%, 30% or 50% instead changes rates by at most {max_window:.5f}%. Independent volume-loss rates agree within {max_volume:.5f}%.

At q=500 and 40% AP, a twelve-relaxation-time repeat changes the measured speed from {short_rate:.9f} to {long_rate:.9f} cm/s, or **{duration_change:+.5f}%**, with late drift {long_drift:+.5f}%. The simulated durations are {short_duration:.9g} and {long_duration:.9g} s. The deeper lower boundary preserves final cold-boundary separation, so this checks duration and domain depth together at one representative point. It does not replace a standard sweep point.

Solver/analytic rate differences span {min_solver:+.4f}% to {max_solver:+.4f}%; incident-flux errors span {min_flux_error:+.4f}% to {max_flux_error:+.4f}%. Rate-fit scatter is not a statistical uncertainty estimate. The duration and window checks assess steadiness, not mesh convergence or physical-model accuracy. Chen's curve stroke resolution is approximately 0.005 cm/s; extra digits preserve reproducibility. The source consists of model curves and DNS, not experimental measurements: [Chen et al. (2002)](https://doi.org/10.1016/S1540-7489(02)80357-1).

## Reproduction and provenance

Cheaper `gpt-5.6-luna` agents launched the simulations and recorded actual terminal exits. The root/current model performed calibration, raw-output extraction, audits, statistics and plotting. All accepted cases have matching input/executable hashes, successful receipts, finalized solver logs, full requested durations, matching frozen parameters and independently recomputed closures. No production solver code or executable changed for this follow-up.

`parameters_baseline.json` preserves the prior accepted parameter file byte for byte. `fixed_properties.json` fixes the new binder cp and all retained properties. `parameters_frozen.json`, `freeze_record.json`, `parameter_comparison.csv`, `comparison.csv`, `error_summary.csv`, execution audits, `fit_window_sensitivity.csv`, `duration_comparison.json` and `provenance.json` retain the results and checks. Earlier reports and figures remain intact. See [REPRODUCE.md](REPRODUCE.md) for commands.
