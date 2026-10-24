# STOPPED: pure-HTPB constituent rerun of the Chen heat-flux sweep

The user stopped this parameter study because its agreement with Chen Figure 4 became worse. Nine baseline flux/composition combinations completed; the two q=500, t=0.5 refinements were stopped and their partial outputs are excluded. No new simulations are requested by this report. The completed baselines differ from their own analytic Eq. 14+19 closure by -3.90% to -2.25%. This documents numerical consistency for the selected parameters, without a completed refinement assessment.

The user's new requested direction is calibration to experimental pure-HTPB data. Chen Figure 4 supplies model/DNS endpoints, while the supplied Gross paper's propellant experimental comparisons concern AP/HTPB mixtures. Neither is a directly interchangeable experimental pure-HTPB calibration target. The specific pure-HTPB dataset and its conditions must be identified before fitting; no experimental calibration is performed here.

The new material set is an AP-free HTPB constituent engineering baseline, with AP mixed in once through the homogeneous model. It is not fitted to Chen Figure 4 or independently validated for a particular cured HTPB formulation. Consequently, changed agreement with that figure measures changed material assumptions as well as numerical error. The prior fitted study remains intact in the parent analysis directory.

## Parameters and provenance

| Constituent | ρ [kg/m³] | cp [J/(kg K)] | k [W/(m K)] | Q [J/kg] | A [m/s] | E/R [K] |
|---|---|---|---|---|---|---|
| HTPB | 920 | 2418.29 | 0.13 | -1.2552e+06 | 10.36 | 7500 |
| AP | 1950 | 1297.9 | 0.4186 | -418400 | 948 | 11000 |

T0=300 K. `../../reference/pure_htpb.json` is the exact parameter record. Densities, heat capacities, and conductivities retain the original input's provisional pure-constituent values. Chen (2002), Figs. 5–6 supplies the Arrhenius endpoints. Gross et al. (2013), section 4.2 supplies the selected endothermic heats, −300 cal/g for HTPB and −100 cal/g for AP (1 cal/g=4184 J/kg). The sources are the user-supplied PDFs; no additional scientific references were used.

In the prior Figure-4-fitted realization, HTPB cp=1255.24 J/(kg K), Q=−196590 J/kg and k=0.271382 W/(m K). The new HTPB heat capacity is larger and its decomposition heat sink much stronger. At the same incident heat flux, those changes lower the rate; conductivity also changes the thermal length and numerical setup. AP parameters change much less. Neither dataset identifies uniquely the physical properties used by Chen.

At q=500 for pure HTPB, an analytic-only change of cp alone lowers the prior 1.75525 cm/s rate to 1.07690; changing Q alone gives 1.01619; changing both gives 0.74864 cm/s. Both thermal changes matter. Conductivity does not enter this steady surface energy/kinetic balance directly, although it changes the resolved thermal profile and its numerical errors.

The independent analytic calculation uses q=ρr[cp(Ts−T0)−Q], r=A exp[−(E/R)/Ts]. For a blend, AP mass fraction is w=ρ_AP t/[ρ_AP t+ρ_B(1−t)]; ρ, cp and Q use the existing arithmetic mass-weighted rule, while log A and E/R use volume fraction t. Eq. 15 instead harmonically averages the two pure-component rates at equal q. The density-convention ambiguity documented in `../REFERENCE_NOTES.md` remains; it is not resolved by changing the binder properties.

The change also reverses the pure-endpoint ordering. At q=500, the new analytic HTPB rate is 0.74864 cm/s and AP is 0.85292 cm/s, while the prior HTPB rate was 1.75525 cm/s. New HTPB therefore regresses more slowly than AP under the same prescribed heat flux. The new t=0.5 Eq. 14+19 rate, 0.71847 cm/s, is below both pure endpoints. This interior depression comes from the chosen combined thermal/kinetic mixing closure, conditional on its arithmetic mass-weighted density rule; it is not a feature forced by Eq. 15, whose harmonic mean remains between the pure rates. With unequal constituent densities, arithmetic mass weighting of density raises blend volumetric heat capacity relative to the simple volume-weighted constituent capacity, contributing to slower mixed regression. These conditional model predictions should not be interpreted as independently validated physical trends or a correction to Chen's article.

An **analytic-only density sensitivity** holds cp, Q, A and E/R fixed at the same blend values and substitutes volume-weighted density. At q=500, t=0.5 this gives 0.80257 cm/s, between the pure endpoints, compared with 0.71847 cm/s for the executed rule. The arithmetic density is 12.88% larger here. This isolates a substantial density contribution to the interior depression without attributing it solely to material kinetics. `density_sensitivity.csv` retains all three fluxes. No simulations use this alternative; the source and all inputs retain the user-specified arithmetic mass-weighted density. This comparison does not identify which convention Chen used.

## Measured comparison

All rates are cm/s. New closure values are recomputed independently in `../compare_pure_htpb.py` and checked against every case's recorded targets. New LowMach rates come from raw AMReX eta=0.5 surface-position fits, not initialized rates. Prior rates are the retained raw-output-derived original summary. Chen values come from the previously extracted PDF vector curves; graphical resolution is approximately 0.005 cm/s, not a confidence interval.

| q [cal/(cm² s)] | t | New LowMach | New Eq. 14+19 | Error [%] | Prior LowMach | Chen Fig. 4 | New vs Chen [%] |
|---|---|---|---|---|---|---|---|
| 200 | 0.0 | 0.31444 | 0.32719 | -3.90 | 0.78089 | 0.80194 | -60.79 |
| 200 | 0.5 | 0.30185 | 0.30986 | -2.59 | 0.44513 | 0.50951 | -40.76 |
| 200 | 1.0 | 0.35702 | 0.36544 | -2.30 | 0.36440 | 0.37280 | -4.23 |
| 500 | 0.0 | 0.72144 | 0.74864 | -3.63 | 1.70534 | 1.75538 | -58.90 |
| 500 | 0.5 | 0.70008 | 0.71847 | -2.56 | 1.01884 | 1.16483 | -39.90 |
| 500 | 1.0 | 0.83331 | 0.85292 | -2.30 | 0.85133 | 0.87151 | -4.38 |
| 1000 | 0.0 | 1.34279 | 1.38992 | -3.39 | 3.04723 | 3.14877 | -57.36 |
| 1000 | 0.5 | 1.31753 | 1.35117 | -2.49 | 1.89790 | 2.16474 | -39.14 |
| 1000 | 1.0 | 1.57783 | 1.61422 | -2.25 | 1.61308 | 1.64937 | -4.34 |

`parameter_comparison.pdf` and `.png` overlay the two material sets, analytic closures, simulations, and Chen's solid curves. `parameter_comparison.csv` includes both simulation-versus-paper and analytic-versus-paper errors, so numerical discrepancy is separated from the material/closure discrepancy. `analytic_curves.csv` retains both reconstructed closures and new Eq. 15 curves.

## Baseline numerical and heating checks

| Setup | ell/delta | dy [µm] | Rate [cm/s] | Closure error [%] | Late drift [%] |
|---|---|---|---|---|---|
| Baseline | 0.250 | 0.3689 | 0.70008 | -2.560 | +0.018 |

Baseline ell=delta/4, dy=ell/8 and delta=k/(ρ cp r_ref); the periodic transverse strip has 32 cells. Runs last six delta/r_ref with approximately forty output intervals. The fit uses the final 40% of physical time. Maximum absolute baseline late-window drift is 0.053%. Surface-position and integrated-solid-volume speeds differ by at most 0.001% across the nine completed cases. The two refinements are incomplete; no convergence or refined-rate claim is made for this parameter set.

Baseline eta=0.5 surface-temperature offsets from the sharp-interface analytic targets span -6.50 to +4.04 K. Distributed kinetics across a finite thermal/interface transition are not identical to an Arrhenius law evaluated at a single interpolated temperature. Fit standard errors describe positional scatter and do not represent total numerical or material uncertainty.

| Setup | Incident flux error [% q] | Gas storage [% q] | Snapshot residual [% q] | Steady-solid residual [% q] |
|---|---|---|---|---|
| Baseline | -1.319 | +0.382 | +2.586 | +2.197 |

The input-only fixed gas slab supplies integrated heat; gas-side conductive flux is measured at eta=10⁻⁶ using harmonic face conductivity. The artificial conducting layer retains k=100 W/(m K), cp=1000 J/(kg K), R=319.787 J/(kg K), P=100 MPa and no temperature advection. Gas storage and incident-flux errors quantify the actual delivered heating. The thermal residual uses source, latent heat from measured volume loss, lower-boundary loss and reconstructed thermal storage. It is a snapshot quadrature diagnostic; the separate steady-solid residual replaces solid storage by measured speed times the thermal-profile gradient. Sparse output can bias moving-interface storage even when solver steps are small. These diagnostics do not establish independence from surrogate conductivity or pressure.

Across the nine baselines, incident-flux errors span -2.19% to -1.18% of prescribed heat input, and gas storage spans 0.35% to 0.65%. These measured surrogate effects are part of the remaining speed bias; they should not be treated as material disagreement with the paper.

## Execution evidence and scope

Every selected baseline run has returncode=0, an input SHA-256 matching both execution receipt and case metadata, and the exact new material dictionary. Each solver log records AMReX finalization after reaching the requested duration. The last raw plot is within one final solver step of that time and at least 99.9% of the requested duration; some output writers retain the preceding step rather than the final step. `execution_audit.json` records commands, hashes, durations and independent-closure checks; `provenance.json` hashes the two source PDFs, parameter sets, original simulation summary and digitized curves. All nine completed receipts match the current executable hash: `a01a7fffe7f5b778537a877430d5e0f0e3d0907659de415ca46ab5d7f3fd475b`. `stopped_cases.json` identifies the two excluded partial refinements. New extraction caches, summaries and figures are confined to this directory. The old run outputs, summaries and report are retained.

This is a planar condensed-phase heat-flux study using the existing executable, frozen gas chemistry, and one product gas. It does not reproduce heterogeneous DNS or validate reconciliation with Gross's premixed gas chemistry, composition-dependent gas reactions, flame heat feedback, or experimental pressure-dependent burning rates. Provisional cp/k and the changed constituent basis require separate physical validation before using this as a predictive combustion model.

## Reproduce postprocessing

From the repository root, with numpy/scipy/matplotlib/yt available:

```bash
python validation_homogeneous_chen2002/analysis/compare_study.py --tag _pure_htpb --output validation_homogeneous_chen2002/analysis/pure_htpb --refresh
python validation_homogeneous_chen2002/analysis/compare_pure_htpb.py
```

`simulation_summary.csv` / `.json` contain the nine completed baselines; `timeseries_*.csv` retain raw-field measurements. `numerical_diagnostics.pdf` / `.png` show the baseline key-case thermal checks. The original `write_report.py` must not be run on mixed material datasets. These commands only regenerate analysis of retained completed outputs; they do not resume the stopped simulations.
