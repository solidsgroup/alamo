# Pure-binder Arrhenius calibration with fixed user properties

For the existing LowMach calibration setup, the final fitted values are **A_binder = 24.608333 m/s** and **E_binder/R = 5568.8506 K**, equivalent to **E_binder = 46.302 kJ/mol** (11.06644 kcal/mol). The physical speed law is `r = A_binder exp[-(E_binder/R)/Ts]`.

The direct analytic fit, before correcting for the finite interface and heating surrogate, gives **A_binder = 23.2122709 m/s**, **E_binder/R = 5602.09895 K**, and **E_binder = 46.578442 kJ/mol**. Use this pair for Chen's ideal steady surface balance; use the final pair to reproduce the executed LowMach setup. Numerical correction is conditional on that setup, and is not a separate material measurement.

## Fixed properties and source

Density is 920 kg/m³ (0.92 g/cm³), conductivity 0.213 W/(m K), specific heat 2130 J/(kg K), and Q = −276144 J/kg (−66 cal/g). These are exactly the user's supplied values; only A and E/R were fitted. Initial/deep-solid temperature remains 300 K as in the existing study. All AP values remain those in `reference/pure_htpb.json`. Pure-binder calibration alone does not validate predictions for mixtures.

Targets are the t=0 endpoints of Figure 4 in Chen et al. (2002), printed p. 2927 / PDF page 5, extracted from the original vector paths in `/home/esandall/Downloads/homogeneous_model_paper.pdf`. Its SHA-256 matches the existing digitization provenance. [Publisher record](https://doi.org/10.1016/S1540-7489(02)80357-1). These are model/DNS curves, not experimental pure-binder measurements. Fit fluxes are 200 and 1000 cal/(cm² s); 500 is held out throughout fitting.

## Method and results

Chen equations (12)–(14) give `q = rho*r*[cp*(Ts-T0) - Q]`. With q converted using 1 cal/(cm² s) = 41840 W/m², each fitting rate fixes `Ts = T0 + [q/(rho*r) + Q]/cp`. Solve the two equations `ln(r) = ln(A) - (E/R)/Ts`. Initial fitting temperatures are 702.8475 and 848.4380 K. Conductivity fixes the thermal profile length `delta = k/(rho*cp*r)`; it does not enter this steady endpoint inversion.

The existing `calibrate_pure_htpb.py` now accepts `--fixed-parameters` so the immutable baseline can be supplied explicitly. Subsequent stages divide each target rate by its measured LowMach/analytic rate ratio and repeat the same two-equation inversion. Both endpoint runs must have late speed drift below 0.2%; final endpoint errors must be below 1%. The held-out case is generated only after the parameters are frozen.

| q [cal/(cm² s)] | Use | Chen [cm/s] | Initial analytic [cm/s] | Final LowMach [cm/s] | LowMach error [%] | LowMach Ts [K] |
|---|---|---|---|---|---|---|
| 200 | fit | 0.801938 | 0.801938 | 0.801944 | +0.0008 | 698.15 |
| 500 | held out | 1.755380 | 1.753089 | 1.753105 | -0.1296 | 774.15 |
| 1000 | fit | 3.148773 | 3.148773 | 3.148832 | +0.0019 | 842.81 |

The initial analytic held-out error is -0.1305%; final LowMach held-out error is -0.1296%. All rates measured from LowMach are obtained from raw eta=0.5 interface-position slopes over the final 40% of each run, with a separate volume-loss check.

| Stage | A [m/s] | E/R [K] | q | LowMach [cm/s] | Target error [%] |
|---|---|---|---|---|---|
| stage00 | 23.2122709 | 5602.099 | 1000 | 3.100144 | -1.5444 |
| stage00 | 23.2122709 | 5602.099 | 200 | 0.789963 | -1.4933 |
| stage01 | 24.6083330 | 5568.851 | 1000 | 3.148832 | +0.0019 |
| stage01 | 24.6083330 | 5568.851 | 200 | 0.801944 | +0.0008 |

## Numerical interpretation and uncertainty

The inherited setup uses ell/delta=0.25, eight cells per ell, 32 transverse cells, dt scale=0.2, and six thermal relaxation times. Heating is supplied by a fixed gas slab, gas conductivity 100 W/(m K), gas cp=1000 J/(kg K), pressure 100 MPa, frozen chemistry and temperature advection disabled. It is a prescribed-flux condensed-phase surrogate. Final incident-flux error ranges from -0.638% to -0.526%; maximum late speed drift is 0.0183%. Gas thermal storage ranges from 0.477% to 0.573% of input. Volume-loss and interface-position rates differ by at most 0.0005%. These numerical effects are partly absorbed by the final fitted pair; this is not a grid-convergence study.

One printed stroke corresponds to about 0.0050 cm/s. Perturbing both fitting endpoints independently by ±one stroke gives an analytic-fit envelope A=17.917–30.637 m/s and E/R=5390.5–5828.9 K. This is graphical sensitivity, not a statistical confidence interval; extra stored decimals support reproducibility only.

## Using the parameters

`parameters_frozen.json` contains the final constituent parameters for `run_study.py`; `parameters_00.json` contains the direct analytic fit. `activation_temperature` in LowMach is E/R in kelvin, not E in J/mol. Physical A in m/s must be normalized for the chosen phase field: for the study's lambda=mobility=w1=1, w12=2, kappa=3*ell² configuration, `rate_multiplier = A/(1.5*ell)`. The generator performs that conversion for each input. A is not directly interchangeable with `rate_multiplier`, and changes to the phase-field parameters require their corresponding conversion.

For Chen's existing volume-fraction mixing rule, `ln(A_blend)=(1-t)*ln(A_binder)+t*ln(A_AP)` and `(E/R)_blend=(1-t)*(E/R)_binder+t*(E/R)_AP`. Supply pure constituent values once. No AP/blend data enter this fit.

`comparison.csv`, `comparison.pdf`, `simulation_summary.json`, `stage_history.csv`, `execution_audit.json` and `provenance.json` record the results. Previous rejected studies and production input files are preserved. See `REPRODUCE.md` for commands.
