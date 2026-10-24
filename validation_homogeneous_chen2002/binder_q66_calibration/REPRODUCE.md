# Reproduce the fixed-property pure-binder calibration

Run from the repository root in an environment with numpy, scipy, matplotlib
and yt. The recorded runs use `/home/esandall/Software/anaconda3/bin/python`
and the existing `bin/lowmach-2d-clang++`; there is no custom solver build.
Use `python` below for that interpreter. The standalone solver needs MPI
initialization permissions, even for one process.

Existing parameter files are deliberately never overwritten. To rerun an
inversion, write to a fresh filename and compare it with the recorded stage.
Existing successful solver runs are retained by `run_study.py --run`.
Changed inputs require a new case tag. Analysis commands do not run simulations.

1. The immutable baseline is `fixed_properties.json`. Fit q=200 and 1000:

```bash
python validation_homogeneous_chen2002/calibrate_pure_htpb.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/fixed_properties.json \
  --fixed-parameters validation_homogeneous_chen2002/binder_q66_calibration/fixed_properties.json \
  --output validation_homogeneous_chen2002/binder_q66_calibration/parameters_00.json
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_00.json \
  --fractions 0 --fluxes 200 1000 --dt-scale .2 --tag _binder_q66_cal00 --run
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _binder_q66_cal00 --output validation_homogeneous_chen2002/analysis/binder_q66/stage00
python validation_homogeneous_chen2002/binder_q66_calibration/analyze.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_00.json \
  --summary validation_homogeneous_chen2002/analysis/binder_q66/stage00/simulation_summary.json
```

2. Correct only A/E using the two measured solver/analytic ratios, then rerun:

```bash
python validation_homogeneous_chen2002/calibrate_pure_htpb.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_00.json \
  --fixed-parameters validation_homogeneous_chen2002/binder_q66_calibration/fixed_properties.json \
  --measurements validation_homogeneous_chen2002/analysis/binder_q66/stage00/simulation_summary.json \
  --output validation_homogeneous_chen2002/binder_q66_calibration/parameters_01.json
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_01.json \
  --fractions 0 --fluxes 200 1000 --dt-scale .2 --tag _binder_q66_cal01 --run
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _binder_q66_cal01 --output validation_homogeneous_chen2002/analysis/binder_q66/stage01
python validation_homogeneous_chen2002/binder_q66_calibration/analyze.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_01.json \
  --summary validation_homogeneous_chen2002/analysis/binder_q66/stage01/simulation_summary.json
```

3. After fitting endpoint errors are below 1% and late speed drift is below
0.2%, freeze the parameter file unchanged as `parameters_frozen.json`.
Run the q=500 held-out case with that file:

```bash
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_frozen.json \
  --fractions 0 --fluxes 500 --dt-scale .2 --tag _binder_q66_frozen --run
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _binder_q66_frozen --output validation_homogeneous_chen2002/analysis/binder_q66/heldout
python validation_homogeneous_chen2002/binder_q66_calibration/analyze.py --final \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_frozen.json \
  --summary validation_homogeneous_chen2002/analysis/binder_q66/stage01/simulation_summary.json \
            validation_homogeneous_chen2002/analysis/binder_q66/heldout/simulation_summary.json
```

During this execution the solver was invoked directly on each generated
`runs/CASE/input`, redirecting output to `runs/CASE/stdout.log`. Return codes
were read from completed process sessions, then recorded in `run.json` with
input and binary SHA-256 digests. Those receipts state that elapsed time is
an orchestration upper bound, not a solver performance measurement. The
standard `--run` route above writes equivalent execution receipts itself.

The final audit verifies unchanged thermal/AP properties, the two fitting
fluxes, matching parameter dictionaries, input/binary hashes, successful
solver finalization, full simulated duration, and an independently recomputed
analytic closure. It combines the corrected fitting endpoints with the one
frozen held-out case. `comparison.csv` records both analytic and LowMach rates.

## Frozen volume-fraction sweep

The requested Figure 4 sweep uses AP volume fractions 0, 0.2, 0.4, 0.6, 0.8
and 1 at all three fluxes. Reuse the three completed pure-binder cases and
generate 15 new cases from the frozen file:

```bash
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_frozen.json \
  --fractions .2 .4 .6 .8 1 --fluxes 200 500 1000 \
  --dt-scale .2 --tag _binder_q66_sweep --run
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _binder_q66_sweep --output validation_homogeneous_chen2002/analysis/binder_q66/sweep
python validation_homogeneous_chen2002/binder_q66_calibration/analyze_sweep.py \
  --summary validation_homogeneous_chen2002/binder_q66_calibration/simulation_summary.json \
            validation_homogeneous_chen2002/analysis/binder_q66/sweep/simulation_summary.json
```

The recorded sweep runs at most five independent single-process cases at
once; the `--run` command above executes sequentially. Both routes use the
same generated inputs. No parameter is adjusted after examining the sweep.
The report, CSV errors and overlay are saved in `binder_q66_calibration/sweep/`.

## Runtime and density diagnostics

`diagnostics.py` recomputes rates over the last 20%, 30%, 40% and 50% of all
18 cases, and evaluates the analytic closure with volume-additive density.
That alternative-density calculation does not change solver inputs or refit
kinetics. The separate twelve-relaxation-time run checks duration sensitivity
at the worst-error blend composition:

```bash
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_frozen.json \
  --fractions .4 --fluxes 500 --dt-scale .2 --relaxations 12 --tag _binder_q66_long12 --run
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _binder_q66_long12 --output validation_homogeneous_chen2002/analysis/binder_q66/long12
python validation_homogeneous_chen2002/binder_q66_calibration/diagnostics.py \
  --long-summary validation_homogeneous_chen2002/analysis/binder_q66/long12/simulation_summary.json
```

The longer run uses a deeper domain to preserve the cold-boundary separation
after the additional regression. It does not replace a point in the original
18-case comparison. `sweep/DIAGNOSTICS.md` reports the measured rate change,
fit-window sensitivity and recommendations for the mixture rules.

Every executed input sets `chemistry.model.type = frozen` and lists only the
binder-regression phase-change mechanism. The audit verifies both settings:
no gas-phase reaction heat is added to the prescribed heat flux. Endothermic
phase-change Q remains part of the surface thermal balance.
