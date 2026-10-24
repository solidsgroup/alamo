# Reproduce pure-AP recalibration and the Figure 4 sweep

Run from the repository root using Python with numpy, scipy, matplotlib and
yt. This execution uses `/home/esandall/Software/anaconda3/bin/python` and the
unchanged `bin/lowmach-2d-clang++` containing volume-additive mixture density.
The executable requires MPI initialization permission. No solver rebuild is
needed for this parameter-only follow-up.

`parameters_baseline.json` is an exact copy of the previously frozen binder
calibration. All binder properties and AP thermal properties remain fixed.
Only AP A and E/R change. Calibration scripts refuse to overwrite parameter
files; choose a new filename/tag when repeating a fit.

First compute the direct analytic pair for reference, then fit using the
already measured pure-AP solver/analytic ratios at q=200 and 1000:

```bash
python validation_homogeneous_chen2002/ap_endpoint_calibration/calibrate.py \
  --parameters validation_homogeneous_chen2002/ap_endpoint_calibration/parameters_baseline.json \
  --output validation_homogeneous_chen2002/ap_endpoint_calibration/parameters_analytic.json
python validation_homogeneous_chen2002/ap_endpoint_calibration/calibrate.py \
  --parameters validation_homogeneous_chen2002/ap_endpoint_calibration/parameters_baseline.json \
  --measurements validation_homogeneous_chen2002/binder_q66_calibration/density_correction/simulation_summary.json \
  --output validation_homogeneous_chen2002/ap_endpoint_calibration/parameters_00.json
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/ap_endpoint_calibration/parameters_00.json \
  --fractions 1 --fluxes 200 1000 --dt-scale .2 --tag _ap_cal00 --run --timeout 3600
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _ap_cal00 --output validation_homogeneous_chen2002/ap_endpoint_calibration/stage00
python validation_homogeneous_chen2002/ap_endpoint_calibration/freeze.py \
  --parameters validation_homogeneous_chen2002/ap_endpoint_calibration/parameters_00.json \
  --summary validation_homogeneous_chen2002/ap_endpoint_calibration/stage00/simulation_summary.json
```

Freezing requires two successful, fully completed pure-AP runs, endpoint
errors below 1%, and late-window drift below 0.2%. If another solver correction
is needed, repeat `calibrate.py` with the latest parameter file and measured
summary, using a fresh stage filename and run tag. Freeze that accepted stage.
The q=500 point and all interior mixtures remain excluded from fitting.

After freezing, run the 16 remaining sweep cases (the two fresh accepted
calibration runs provide the other two endpoints) and one duration check:

```bash
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/ap_endpoint_calibration/parameters_frozen.json \
  --fractions 0 .2 .4 .6 .8 --fluxes 200 1000 --dt-scale .2 \
  --tag _ap_frozen --run --timeout 3600
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/ap_endpoint_calibration/parameters_frozen.json \
  --fractions 0 .2 .4 .6 .8 1 --fluxes 500 --dt-scale .2 \
  --tag _ap_frozen --run --timeout 3600
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/ap_endpoint_calibration/parameters_frozen.json \
  --fractions .4 --fluxes 500 --dt-scale .2 --relaxations 12 \
  --tag _ap_frozen_long12 --run --timeout 3600
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _ap_frozen --output validation_homogeneous_chen2002/ap_endpoint_calibration/sweep
python validation_homogeneous_chen2002/ap_endpoint_calibration/analyze.py \
  --summaries validation_homogeneous_chen2002/ap_endpoint_calibration/stage00/simulation_summary.json \
              validation_homogeneous_chen2002/ap_endpoint_calibration/sweep/simulation_summary.json
```

The recorded execution uses cheaper `gpt-5.6-luna` agents to invoke the binary
directly on generated case inputs, capture stdout and record actual successful
terminal exits in `run.json` with input/executable hashes. The current/root
model performs all calibration and numerical postprocessing. The commands
above execute simulations sequentially. No agent changes parameters or fits
regression rates during launching.

`analyze.py` audits all 19 runs, retains the duration check separately from the
18-point sweep, and writes the single-panel PNG/PDF with percentage errors
beside each new LowMach point. It also writes pointwise errors, all-point and
mixture-only error summaries, fitting-window checks and provenance. The
previous density-corrected sweep remains in `../binder_q66_calibration/density_correction/`.
