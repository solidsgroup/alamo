# Standalone Chen et al. (2002) validation study

The selected configuration uses **binder cp=2130 J/(kg K)**, **Q=−66 cal/g**,
and the matching binder/AP Arrhenius calibration from the
[pure-AP recalibration study](ap_endpoint_calibration/REPORT.md).
Its 18-point mean absolute error versus Chen is **0.1314%**, with maximum **0.2790%**.
The user selected this physically plausible constant-cp model after comparing
the cp=2418.29 alternative. See the [selection rationale](CURRENT_PARAMETERS.md),
[current parameters](parameters_current.json), [selection record](current_selection.json),
and [selected Figure 4 plot](ap_endpoint_calibration/figure4_comparison.png).
The [regression GIF](animations/propellant_regression_q500_ap40.gif) shows the
selected 40% AP, q=500 simulation; see its [rendering notes](animations/README.md).
For future study runs, explicitly pass `--parameters validation_homogeneous_chen2002/parameters_current.json`
and `--dt-scale .2` to `run_study.py`; its historical default parameter file is unchanged.

The archived [binder cp restoration and Figure 4 sweep](binder_cp2418_calibration/REPORT.md)
restores binder cp to 2418.29 J/(kg K), retains Q=−66 cal/g and all AP parameters,
and refits binder A and E/R using q=200 and 1000, with q=500 held out.
The 18-point mean absolute error versus Chen increases from 0.1314% to 0.6164%;
the maximum is 1.2143%. The doubled-duration check changes its rate by −0.1343%.
See the [annotated single-panel figure](binder_cp2418_calibration/figure4_comparison.png),
[frozen parameters](binder_cp2418_calibration/parameters_frozen.json),
[complete parameter comparison](binder_cp2418_calibration/parameter_comparison.csv),
and [reproduction commands](binder_cp2418_calibration/REPRODUCE.md).

The preceding [pure-AP recalibration and Figure 4 sweep](ap_endpoint_calibration/REPORT.md)
fits AP A and E/R at q=200 and 1000, with q=500 held out. It retains the binder
calibration, all thermal properties, volume-additive density and frozen gas
chemistry. See the [frozen parameters](ap_endpoint_calibration/parameters_frozen.json)
and [reproduction commands](ap_endpoint_calibration/REPRODUCE.md).

The preceding [density-corrected Figure 4 sweep](binder_q66_calibration/density_correction/REPORT.md)
uses volume-additive density with the same frozen pure-constituent calibration.
Specific heat and Q remain mass-weighted. See its
[reproduction commands](binder_q66_calibration/density_correction/REPRODUCE.md).
The earlier binder-calibration reports below describe the previous arithmetic
mass-weighted density implementation; their outputs are preserved for comparison.

The preceding pure-binder calibration used the user's specified rho=920 kg/m³,
k=0.213 W/(m K), cp=2130 J/(kg K), and Q=−66 cal/g, fitting only A and E/R.
See [the calibration report](binder_q66_calibration/REPORT.md) and
[reproduction commands](binder_q66_calibration/REPRODUCE.md). It uses the
same two-endpoint fitting and solver-correction method, with q=500 held out.
The follow-up [frozen-constituent Figure 4 sweep](binder_q66_calibration/sweep/SWEEP_REPORT.md)
compares six AP volume fractions at each of the three heat fluxes and reports
pointwise errors against Chen's solid curves.
The [mixing and runtime diagnostics](binder_q66_calibration/sweep/DIAGNOSTICS.md)
evaluate volume-additive density without refitting kinetics and check rate
sensitivity to the fitting window and a longer representative run.

The earlier 2026-09-09 [AP-free HTPB literature audit](analysis/HTPB_LITERATURE_REVIEW.md)
documents physical-property sources and the conflict under the older −300 cal/g
baseline. The stopped and rejected stages below remain historical records;
the new user-specified thermal baseline is isolated in `binder_q66_calibration/`.

This directory is intentionally untracked and is not part of the regression
test suite. It holds source provenance, physical reference calculations,
LowMach input decks, simulation outputs,
and the final comparison report.

Target: Figure 4 of Chen, Buckmaster, Jackson, and Massa, *Homogenization
issues and the combustion of heterogeneous solid propellants*, Proceedings
of the Combustion Institute 29 (2002), 2923–2929. The prescribed surface heat
fluxes are 200, 500, and 1000 cal/(cm² s). Reference equations are (12)–(19).

The study uses the existing `bin/lowmach-2d-clang++` executable. Heating is
configured only through `heat_source.ic.expression`: a fixed slab in the gas
supplies a prescribed integrated power, with an insulated upper boundary.
The gas conducts that power to the moving surface; postprocessing must check
the incident flux and sensitivity to gas heat storage and interface width.
This is a condensed-phase thermal surrogate, not a gas-flame calculation.
The original study used the existing solver. The density-correction follow-up
updates the homogeneous density rule in production source and its conservation
tests, as documented in the density-correction report linked above.

The original fitted study uses only the supplied Chen article as a scientific reference. It does not
tabulate all constituent properties. The article-analysis files distinguish
explicit parameters from combinations inferred from Figure 4. Runs using
those fits are article-constrained reconstructions, not independent validation
of the pure-material endpoints. Assumed density ratios must be stated.

The initially created solver copy, custom executable and build script were
withdrawn at the user's request; recoverable files are in
`/tmp/alamo-rejected-flux-adapter-pjwIjk/`. They are not used by this study.

Work allocation: cheaper models ran the isothermal calibration and pressure-
independent heat-flux/composition sweep; a more capable model independently
reads the raw outputs and assesses agreement with the source paper.

## Reproduce

Use the existing solver containing the homogeneous phase-change feature;
there is no study-specific compilation step. From the repository root,
with a Python environment containing numpy, scipy, matplotlib and yt:

```bash
python validation_homogeneous_chen2002/run_study.py --fractions 0 .5 1 --fluxes 200 500 1000 --tag _nx32 --run
python validation_homogeneous_chen2002/run_study.py --fractions .5 --fluxes 500 --dt-scale .1 --tag _nx32_small_dt --run
python validation_homogeneous_chen2002/run_study.py --fractions .5 --fluxes 500 --width-ratio .125 --dt-scale .1 --tag _thin_small_dt --run
python validation_homogeneous_chen2002/analysis/compare_study.py
```

The generator preserves successful cases and refuses to overwrite a deck
with different settings; use a new `--tag` for another configuration.
Alternatively run any generated deck directly with
`bin/lowmach-2d-clang++ path/to/input`. The recorded runs used this direct
invocation, because the coding environment's MPI permission applies to the
solver executable rather than a Python subprocess wrapper.

The 32-cell periodic transverse strip permits sufficient multigrid coarsening
for this planar problem. Its physical width does not prescribe the surface
speed. The gas-layer heat capacity is valid for the ideal-gas EOS; high
conductivity reduces its storage effect. Temperature advection is disabled
for the condensed-phase thermal surrogate. Failed narrow-strip and invalid-
heat-capacity probes are retained for provenance, not used as validation.

See `analysis/STUDY_REPORT.md` for executed cases, numerical sensitivity,
model limitations, and the comparison with the article. Case-level results
are in `analysis/simulation_summary.csv` and `analysis/comparison.pdf`.

## Pure-HTPB constituent rerun

**Stopped at the user's request.** Nine baselines completed; the two
refinements were stopped and their partial outputs are excluded. The
commands below document the setup and are not an instruction to resume it.

`reference/pure_htpb.json` defines a separate engineering baseline with an
AP-free HTPB constituent. The Arrhenius endpoints remain from Chen; the
selected phase-change heats are −300 cal/g for HTPB and −100 cal/g for AP
from section 4.2 of the supplied Gross et al. (2013) paper. Density, heat
capacity and conductivity retain provisional original-input constituent
values. The homogeneous model mixes these pure constituents once. This
parameter set is not fitted to Chen Figure 4 or validated for a specific
cured binder formulation.

Generate the nine baseline cases and two key-case refinements separately:

```bash
python validation_homogeneous_chen2002/run_study.py --parameters validation_homogeneous_chen2002/reference/pure_htpb.json --fractions 0 .5 1 --fluxes 200 500 1000 --tag _pure_htpb_nx32
python validation_homogeneous_chen2002/run_study.py --parameters validation_homogeneous_chen2002/reference/pure_htpb.json --fractions .5 --fluxes 500 --dt-scale .1 --tag _pure_htpb_nx32_small_dt
python validation_homogeneous_chen2002/run_study.py --parameters validation_homogeneous_chen2002/reference/pure_htpb.json --fractions .5 --fluxes 500 --width-ratio .125 --dt-scale .1 --tag _pure_htpb_thin_small_dt
```

Execute the generated inputs with the existing executable as above (or add
`--run` if the environment allows solver subprocesses). Postprocess only
this tagged set into its own directory:

```bash
python validation_homogeneous_chen2002/analysis/compare_study.py --tag _pure_htpb --output validation_homogeneous_chen2002/analysis/pure_htpb --refresh
python validation_homogeneous_chen2002/analysis/compare_pure_htpb.py
```

The final `analysis/pure_htpb/STUDY_REPORT.md` and
`analysis/pure_htpb/parameter_comparison.pdf` compare new raw-output-derived
rates, an independently recomputed new closure, prior simulations and the
digitized Chen curves. `execution_audit.json` records completion and matching
input hashes for the nine completed baselines; `stopped_cases.json` records
the excluded refinements. Original fitted runs and analysis
artifacts are retained; avoid running the unfiltered original report writer
after both material sets exist. Numerical consistency with the new closure
is separate from agreement with Chen. Frozen chemistry and imposed heating
do not validate reconciliation with Gross premixed gas chemistry.

The user's subsequent requested direction is fitting experimental pure-HTPB
data. Chen model/DNS endpoints and Gross's AP/HTPB mixture experimental
comparisons do not by themselves identify that calibration target; the
specific pure-HTPB experimental dataset and conditions remain to be supplied
or selected. No experimental calibration is included in this stopped study.

## Pure-binder calibration to the confirmed Chen target

**Rejected by the user: this attempt adjusted cp/Q, but the requested fit
is to Arrhenius A/E.** The two cp/Q calibration stages are retained as
provenance; their fitted parameters are not the accepted result. The frozen
cp/Q prediction sweep was stopped. Commands below document that rejected
workflow and are not instructions to resume it.

The user subsequently selected Chen Figure 4's pure-binder model/DNS
endpoint as the calibration target. This resolves the target choice for
this separate follow-up; it does not turn those curves into experimental
pure-HTPB measurements.

The archived `calibration_pure/rejected_cp_Q_fit.py` adjusted only binder heat capacity and decomposition
heat using q=200 and 1000 cal/(cm² s), with binder density, conductivity,
Arrhenius parameters and cold temperature fixed. All AP parameters remain
exactly those in `reference/pure_htpb.json`. The q=500 pure-binder point is
held out; mixed-composition data are excluded from fitting. Solver-informed
corrections are conditional on this numerical setup and heating surrogate.

The retained rejected-fit file is `calibration_pure/rejected_cp_Q_parameters.json`;
the original stage dictionary remains in `calibration_pure/parameters_01.json`. After one
numerical correction, cp=1227.562955 J/(kg K) and Q=−199369.444793 J/kg
(−47.65044 cal/g); both fitting endpoint errors are below 0.002%. These tiny
fit residuals do not imply comparable precision in Chen's graphical data.
No held-out or mixture result is used to change the frozen parameters.

Parameter stages are preserved in `calibration_pure/parameters_*.json`;
independent raw-output measurements and execution audits are isolated under
`analysis/calibrated_htpb/stage*/`. After parameters are frozen, the two
completed fitting endpoints are reused alongside seven new runs for the
held-out pure-binder point and six AP/blend predictions. The final report
and overlay are `analysis/calibrated_htpb/STUDY_REPORT.md` and
`analysis/calibrated_htpb/comparison.pdf`. They preserve the distinction
between fitting, held-out assessment and unfitted predictions.

For any completed fitting stage, independently extract its tag and audit
the selected parameter file with:

```bash
python validation_homogeneous_chen2002/analysis/compare_study.py --tag _htpb_cal00 --output validation_homogeneous_chen2002/analysis/calibrated_htpb/stage00
python validation_homogeneous_chen2002/analysis/report_calibration.py --stage-only --summary validation_homogeneous_chen2002/analysis/calibrated_htpb/stage00/simulation_summary.json --parameters validation_homogeneous_chen2002/calibration_pure/parameters_00.json
```

Final reporting for the rejected cp/Q workflow is disabled. Its final
prediction sweep was interrupted and must not be presented as complete.

No analysis command starts simulations. The original reconstruction and
stopped engineering-baseline artifacts remain intact.

The corrected requested direction holds all thermal properties at
`reference/pure_htpb.json` and fits only binder Arrhenius A/E to the Chen
pure-binder target. The q=200 target implies a surface temperature below
300 K with those fixed heats/capacities, so the physical regime requires
explicit attention before treating an algebraic fit as a predictive model.
See `analysis/arrhenius_htpb/FEASIBILITY.md`; no Arrhenius simulation or
parameter file is created by that feasibility calculation.

## Arrhenius-only fitting workflow

**Rejected and stopped by the user.** The fitted low-temperature regime
is not acceptable for the intended physical pyrolysis model. Existing
stages and outputs are preserved, but no frozen Arrhenius parameter set
is published, no AP/composition sweep proceeds, and no further extraction
or calibration advancement is requested. See
`analysis/arrhenius_htpb/STATUS.md`. Commands below are historical workflow
documentation, not instructions to resume it.

The current `calibrate_pure_htpb.py` adjusts only binder A and E/R, asserting
that every thermal property and every AP property remains fixed. Parameter
stages are retained in `arrhenius_calibration/`. The q=200 and 1000
pure-binder rates define the fit; q=500 remains held out, and AP/blend data
are not fitted. The below-T0 surface regime is reported explicitly.
The standalone cutoff remains its default 0 K. The original
`input.lm.ap_htpb` uses a 360 K cutoff, which is incompatible with the low-
and middle-flux target temperatures under these fixed thermal properties.
The fitted coefficients therefore do not directly transfer to that input
with its cutoff intact; no cutoff is changed by the study.

Independently extract and audit a completed fitting stage with:

```bash
python validation_homogeneous_chen2002/analysis/compare_study.py --tag _arrhenius_cal00 --output validation_homogeneous_chen2002/analysis/arrhenius_htpb/stage00
python validation_homogeneous_chen2002/analysis/report_arrhenius.py --stage-only --summary validation_homogeneous_chen2002/analysis/arrhenius_htpb/stage00/simulation_summary.json --parameters validation_homogeneous_chen2002/arrhenius_calibration/parameters_00.json
```

The separate `report_arrhenius.py` combines final-stage endpoints with
seven frozen-parameter runs using `--summary`, `--frozen-summary` and
`--parameters`. It writes `analysis/arrhenius_htpb/STUDY_REPORT.md`, the
comparison plot, stage history and execution audit only after all nine
final measurements are available. It never changes parameters or launches
simulations, and it does not re-enable the rejected cp/Q report.
