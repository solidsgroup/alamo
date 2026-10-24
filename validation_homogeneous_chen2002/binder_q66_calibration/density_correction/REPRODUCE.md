# Reproduce the sweep with corrected density averaging

From the repository root, rebuild the existing solver and unit executable:

```bash
make -j4 bin/lowmach-2d-clang++ bin/test-2d-clang++
bin/test-2d-clang++
bin/lowmach-2d-clang++ tests/LMHomogeneousBinder/input plot_file=/tmp/chen-volume-density-homogeneous-test
python tests/LMHomogeneousBinder/test /tmp/chen-volume-density-homogeneous-test
```

Use a Python environment with numpy, scipy, matplotlib and yt. The recorded
analysis uses `/home/esandall/Software/anaconda3/bin/python`. MPI must be able
to initialize for the executables. The unit suite and uniform one-step
regression verify volume, heat capacity, phase-change energy, product mass
conservation and thermal coupling.

The calibrated constituents are unchanged in
`binder_q66_calibration/parameters_frozen.json`. Generate/run all 18 points:

```bash
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_frozen.json \
  --fractions 0 .2 .4 .6 .8 1 --fluxes 200 500 1000 \
  --dt-scale .2 --tag _binder_q66_volume --run
python validation_homogeneous_chen2002/run_study.py \
  --parameters validation_homogeneous_chen2002/binder_q66_calibration/parameters_frozen.json \
  --fractions .4 --fluxes 500 --dt-scale .2 --relaxations 12 \
  --tag _binder_q66_volume_long12 --run
```

Existing successful matching runs are preserved. Changed inputs require a
new tag. In this execution, three `gpt-5.6-luna` agents each owned one heat-flux
group. They invoked the
binary directly on each case's input, captured stdout, observed exit codes,
and wrote `run.json` receipts with input/executable digests. The q=500 agent
also ran the longer-duration case. No launch agent extracted rates or
postprocessed numerical data. The commands above run sequentially.

The root/current model performs the raw-output extraction and comparison:

```bash
python validation_homogeneous_chen2002/analysis/compare_study.py \
  --tag _binder_q66_volume --output validation_homogeneous_chen2002/analysis/binder_q66/volume
python validation_homogeneous_chen2002/binder_q66_calibration/density_correction/analyze.py \
  --summary validation_homogeneous_chen2002/analysis/binder_q66/volume/simulation_summary.json
```

The tag includes the duration check. `analyze.py` separates the 18 six-time
sweep cases from the one twelve-time check, audits both, compares the new
results with the previous sweep, and writes the report, plot and CSV files
in this directory. All rates come from raw eta=0.5 surface-position fits;
the analytic reference does not determine the measured surface motion.
The plot includes Chen's vector-extracted solid curves and DNS markers.

The former executable and edited files are retained in `baseline_snapshot/`,
with their hashes. Earlier reports describe the previous arithmetic
mass-weighted density implementation; their original outputs remain intact.
The new generator's default is volume-additive density and records that rule
in every new case. The historical `compare_pure_htpb.closure` helper retains
its legacy default so earlier analyses remain interpretable; this report
explicitly selects `volume_density=True` for corrected cases.
