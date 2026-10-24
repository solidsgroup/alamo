# Reproduce binder cp-only recalibration and the Figure 4 sweep

Run from the repository root with Python containing numpy, scipy, matplotlib
and yt. The recorded environment is `/home/esandall/Software/anaconda3/bin/python`.
Use the unchanged `bin/lowmach-2d-clang++` with volume-additive density.
Its MPI initialization requires the existing solver execution permission.

The baseline is an exact copy of `../ap_endpoint_calibration/parameters_frozen.json`.
`fixed_properties.json` restores binder cp to 2418.29, keeps Q=-276144 J/kg,
retains all other thermal properties and AP kinetics, and archives the earlier
binder fit history. Only binder A and E/R are refitted. Existing files and runs
are preserved: use fresh output filenames and tags for a new replication.

For stages 00 (direct fit) and 01 (solver correction):

```bash
python validation_homogeneous_chen2002/calibrate_pure_htpb.py \
  --parameters validation_homogeneous_chen2002/binder_cp2418_calibration/fixed_properties.json \
  --fixed-parameters validation_homogeneous_chen2002/binder_cp2418_calibration/fixed_properties.json \
  --output validation_homogeneous_chen2002/binder_cp2418_calibration/parameters_00.json
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_cp2418_calibration/parameters_00.json \
  --fractions 0 --fluxes 200 1000 --dt-scale .2 --tag _cp2418_cal00 --run --timeout 3600
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _cp2418_cal00 --output validation_homogeneous_chen2002/binder_cp2418_calibration/stage00
python validation_homogeneous_chen2002/binder_cp2418_calibration/stage.py \
  --parameters validation_homogeneous_chen2002/binder_cp2418_calibration/parameters_00.json \
  --summary validation_homogeneous_chen2002/binder_cp2418_calibration/stage00/simulation_summary.json
python validation_homogeneous_chen2002/calibrate_pure_htpb.py \
  --parameters validation_homogeneous_chen2002/binder_cp2418_calibration/parameters_00.json \
  --fixed-parameters validation_homogeneous_chen2002/binder_cp2418_calibration/fixed_properties.json \
  --measurements validation_homogeneous_chen2002/binder_cp2418_calibration/stage00/simulation_summary.json \
  --output validation_homogeneous_chen2002/binder_cp2418_calibration/parameters_01.json
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_cp2418_calibration/parameters_01.json \
  --fractions 0 --fluxes 200 1000 --dt-scale .2 --tag _cp2418_cal01 --run --timeout 3600
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _cp2418_cal01 --output validation_homogeneous_chen2002/binder_cp2418_calibration/stage01
python validation_homogeneous_chen2002/binder_cp2418_calibration/stage.py \
  --parameters validation_homogeneous_chen2002/binder_cp2418_calibration/parameters_01.json \
  --summary validation_homogeneous_chen2002/binder_cp2418_calibration/stage01/simulation_summary.json --freeze
```

Freeze only after the two completed endpoints pass the execution audit,
errors are below 1%, drift is below 0.2%, and temperatures are above T0.
If another correction is needed, use the latest measurements with a fresh
stage filename and tag. Neither q=500 nor any mixture enters the fit.

Generate/run the 16 other sweep points plus a duration check after freezing:

```bash
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_cp2418_calibration/parameters_frozen.json \
  --fractions .2 .4 .6 .8 1 --fluxes 200 1000 --dt-scale .2 --tag _cp2418_frozen --run --timeout 3600
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_cp2418_calibration/parameters_frozen.json \
  --fractions 0 .2 .4 .6 .8 1 --fluxes 500 --dt-scale .2 --tag _cp2418_frozen --run --timeout 3600
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_cp2418_calibration/parameters_frozen.json \
  --fractions .4 --fluxes 500 --dt-scale .2 --relaxations 12 --tag _cp2418_frozen_long12 --run --timeout 3600
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _cp2418_frozen --output validation_homogeneous_chen2002/binder_cp2418_calibration/sweep
python validation_homogeneous_chen2002/binder_cp2418_calibration/analyze.py \
  --summaries validation_homogeneous_chen2002/binder_cp2418_calibration/stage01/simulation_summary.json \
              validation_homogeneous_chen2002/binder_cp2418_calibration/sweep/simulation_summary.json
```

For the recorded direct launches, first generate the three sweep groups above
without `--run`, then run `python validation_homogeneous_chen2002/binder_cp2418_calibration/manifest.py`
before launching the binaries. This records input hashes, executable hash,
parameter hash and case settings in `launch_manifest.json`; the final audit
verifies that manifest.
The recorded execution uses cheaper `gpt-5.6-luna` agents to invoke the binary
directly on generated inputs, capture stdout, and write `run.json` only after
observing an actual terminal exit. The commands above run sequentially;
the root/current model performs all calibration and postprocessing.

The final analysis requires exactly 18 standard points and one separate
twelve-time repeat. It writes a single-panel PNG/PDF with signed errors beside
LowMach points, pointwise and aggregate errors, parameter units, fitting-window
sensitivity, duration sensitivity, and full execution/provenance audits.
