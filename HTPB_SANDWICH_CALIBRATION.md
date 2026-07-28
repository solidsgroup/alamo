# Calibrating the `fullfeedback` HTPB regression rate against sandwich data

Fits **one** `(pre_exponential, activation_temperature)` pair for
`HTPB_pyrolysis`'s `Model::PhaseField::FullFeedback` phase field, using the
AP/HTPB sandwich geometry in `input.lm.ap_htpb_fullfeedback`, against the
experimental sandwich burning-rate-vs-pressure curve in
`HTPB_sandwich_reg_rate_chorpening2000_summary.csv` (0.2-3.2 MPa). This
mirrors the AP monopropellant calibration pipeline documented in
`FULLFEEDBACK_CALIBRATION.md` -- read that first if you haven't; this reuses
its exact pattern (template input -> render -> pressure sweep -> regression
rate -> `least_squares`).

**`AP_decomposition`'s fullfeedback parameters are held fixed** at the
values already calibrated against pure AP self-deflagration data
(`pre_exponential = 0.0006072874693037711 1/Pa/s`,
`activation_temperature = 1905.3816062268543 K`, see
`input.lm.ap_monopropellant_fullfeedback`). Only `HTPB_pyrolysis`'s own
fullfeedback parameters are fit here.

## Why HTPB's phase field changed from `allencahn` to `fullfeedback`

The pre-existing `input.lm.ap_htpb_fullfeedback` modeled AP with
`fullfeedback` (Arrhenius, pressure-scaled mobility) but HTPB with
`allencahn` (constant mobility, no pressure dependence) -- an ad hoc,
uncalibrated placeholder (`mobility = 0.01 1/Pa/s`). To let HTPB's
regression rate be calibrated the same way AP's was, the template
(`input.lm.ap_htpb_fullfeedback.template`) switches `HTPB_pyrolysis` to
`phase_field.type = fullfeedback` too, with its own
`pre_exponential`/`activation_temperature` as the two free parameters. The
mechanism-level `HTPB_pyrolysis.phase_change.activation_temperature` is
zeroed (mirroring `AP_decomposition`'s own `= 0.0_K`) to avoid double-
counting the Arrhenius factor that now lives inside `fullfeedback` itself.

`sigma`/`epsilon`/`driving_force` are **not** fit parameters -- they are
interface-regularization terms tied to this case's mesh resolution
(`amr.max_level = 2` on a 32x32 base grid), kept identical between AP and
HTPB and identical to the pre-existing `input.lm.ap_htpb_fullfeedback`
values (`epsilon = 80 um`, coarser than the 1D monopropellant calibration's
`20 um` because this is a full 2D case and 80 um is adequately resolved by
this grid's finest AMR level).

## Target data: pure-binder sandwich, not the fine-AP/binder matrix

Earlier in this session, AP/PBAN "matrix" sandwich data (fine AP particles
mixed into the binder, from Chakravarthy et al. 2003/2004) was digitized
into `sandwich_reg_rate_chakravarthy2003.csv`. That is **not** used here: the
existing `input.lm.ap_htpb_fullfeedback` geometry is a *pure* AP layer next
to a *pure* HTPB layer (no AP particles mixed into the binder), so the
physically correct comparison is pure-binder sandwich data:

- **Chorpening, Knott & Brewster, "Flame Structure and Burning Rate of
  AP/HTPB Propellant Sandwiches," Proc. Combust. Inst. 28:847-853 (2000)** --
  Fig. 3, pure HTPB binder lamina (50-450 um wide), 0.2-3.2 MPa, reports
  `r_b ~ P^0.4`, roughly independent of binder width above 100 um.

Fig. 3 was digitized the same pixel-calibrated way as the earlier sandwich
CSV (axis ticks located precisely, marker blobs centroided, overlapping
markers split by row-width profiling where the plot allowed it) into two
files:

- `HTPB_sandwich_reg_rate_chorpening2000.csv` -- all 13 digitized points,
  with the binder-width bin(s) each point/merged-blob belongs to and notes
  on which points were averaged from overlapping markers.
- `HTPB_sandwich_reg_rate_chorpening2000_summary.csv` -- 4 representative
  `(pressure, rate)` points (one per pressure cluster in Fig. 3), each
  averaged across the >=70 um binder-width bins present at that pressure
  (the paper's own conclusion is that rate is width-independent above
  ~100 um, so this collapses the width dimension the model doesn't resolve).
  The isolated 50 um point at 1.45 MPa was excluded from its cluster's
  average -- the paper explicitly calls it anomalous, "thin binder"
  behavior trending toward the pure-AP rate, not representative of this
  case's ~200 um HTPB lamina. **This summary CSV is what
  `optimize_htpb_fullfeedback.py` fits against by default.**

**Caveat:** this pressure range (0.2-3.2 MPa) sits mostly *below* the
2.76-13.79 MPa range `AP_decomposition`'s own parameters were calibrated
against (`AP_reg_rate.csv`). Evaluating AP's fixed, already-calibrated law
down at 0.2-1.5 MPa is an extrapolation below its own calibration floor, on
top of the usual sandwich-vs-monopropellant modeling gap. Treat the
resulting HTPB fit as approximate for this reason.

## Second dataset: internal-group 100 um HTPB laminate

`HTPB_sandwich_reg_rate_group_100um.csv` adds a second, independent
AP/HTPB sandwich dataset from prior work in this group: an HTPB laminate of
thickness `1.0e-4 m` (100 um) embedded in an AP matrix, 0.8-4.0 MPa. Source:
a "Sandwich - Experiment" vs. "Sandwich - Model" plot from that prior work
(the "Model" curve is that earlier study's own fit, not used here -- only
the "Experiment" curve was digitized). The 7 points were read directly off
the plot's gridlines (no markers on the experiment curve to pixel-lock onto,
so this is an approximate reading, not a pixel-precise digitization like the
Chorpening CSV above):

```
Pressure (MPa), Reg Rate (mm/s)
0.8, 2.00
1.0, 2.20
1.5, 2.68
2.0, 3.00
2.5, 3.50
3.0, 4.00
4.0, 5.00
```

This is a **different physical geometry** from the Chorpening-fit case
above: 100 um total HTPB lamina vs. that case's ~200 um. Rather than
stretching one template across two lamina widths, there is a second input
deck, `input.lm.ap_htpb_fullfeedback_100um.template`, resized for this:

- `geometry.prob_lo/hi.x` narrowed from +-0.3 mm to +-0.25 mm, which (with
  the AP/HTPB partition unchanged at `|x|<0.2 mm`) halves each periodic
  HTPB strip from 0.1 mm to 0.05 mm -- joined into one 0.1 mm (100 um)
  lamina across the periodic boundary, same convention as the base
  template. AP's central width (0.4 mm) was left unchanged -- there is no
  independently reported AP matrix width for this dataset, so only HTPB's
  width was resized to match what *is* specified. Revisit if the real AP
  matrix width becomes available.
- `fullfeedback.epsilon` (both mechanisms) reduced from 80 um to 20 um: 80
  um is already close to the *base* template's 100 um HTPB strip width, and
  would exceed half the strip width here (50 um), risking the AP-HTPB and
  HTPB-gas diffuse interfaces overlapping. 20 um matches the AP
  monopropellant calibration's interface width and is well within this
  case's existing mesh resolution (`max_level=2` already resolves ~4.7 um
  cells here, the same interface-to-cell ratio validated in that case).

This deck was smoke-tested this session (see "Verification already done"
below) -- it runs stably and reaches a steady combined-front regression
rate of the right order of magnitude at 1.5 MPa with an arbitrary initial
guess. It has **not** been run through the real optimizer yet.

To fit against this dataset instead of Chorpening's, point
`optimize_htpb_fullfeedback.py` at both the new data and template:

```bash
python scripts/optimize_htpb_fullfeedback.py \
    --workdir /path/to/htpb_100um_calib_run \
    --data HTPB_sandwich_reg_rate_group_100um.csv \
    --template input.lm.ap_htpb_fullfeedback_100um.template \
    --fit-pressures 0.8 1.0 1.5 2.0 2.5 3.0 4.0
```

This produces an independent `(htpb_pre_exponential,
htpb_activation_temperature)` fit from the Chorpening one -- the two
datasets are **not** combined into a single fit, since they correspond to
different physical geometries (lamina width) that the model does not
otherwise parameterize. Compare the two resulting fits once both have been
run; a large discrepancy would suggest the model's width-independence
assumption (inherited from Chorpening's own finding that rate is
binder-width-independent above ~100 um) doesn't hold as well down at 100 um.

## Pieces

- `input.lm.ap_htpb_fullfeedback.template` -- the AP/HTPB sandwich input
  parameterized by operating pressure and HTPB's two fit parameters. AP's
  fullfeedback parameters are hardcoded (not templated). Placeholders are
  filled in by `scripts/render_htpb_sandwich_input.py`.
- `scripts/render_htpb_sandwich_input.py` -- fills in the template; derives
  the pressure-dependent product density the same way the original
  `input.lm.ap_htpb_fullfeedback` does (`rho = P / (319.787 * 700)`).
- `scripts/run_htpb_sandwich_pressure_sweep.sh` -- for one candidate HTPB
  parameter pair, renders and runs one sim per requested pressure **in
  parallel** (each via `mpirun -np 2`, matching the pre-existing sandwich
  input's run instructions), measures each with
  `scripts/regression_rate.py --rate-only --unit mm/s` (tracking the
  combined `rigid_eta` field -- the sum of both AP_solid and HTPB_solid
  volume fractions -- so the reported rate is the overall sandwich surface
  regression rate, not a per-species one), and writes a
  `pressure_mpa,reg_rate_mm_s` CSV.
- `scripts/optimize_htpb_fullfeedback.py` -- `scipy.optimize.least_squares`
  driver over `log10(htpb_pre_exponential)` and
  `htpb_activation_temperature`. By default runs in two stages: first a
  cheap bracket fit against only the lowest and highest `--fit-pressures`
  (2 sims/iteration), then a full fit against all 4 summary-CSV points
  starting from that result. Logs every iteration, then runs a full
  validation sweep and writes the calibrated input plus a sim-vs-experiment
  plot.

## Running the real calibration

```bash
python scripts/optimize_htpb_fullfeedback.py --workdir /path/to/htpb_calib_run
```

Useful flags:
- `--fit-pressures` -- pressures (MPa) used during optimization (default:
  all 4 points in `HTPB_sandwich_reg_rate_chorpening2000_summary.csv`, i.e.
  `0.213 0.467 1.498 3.098`; there is no larger held-out set for this
  dataset, so validation reruns the same 4 points).
- `--pre-exponential0` / `--activation-temperature0` -- initial guess
  (defaults `0.001` 1/Pa/s, `4000` K -- rough starting points, not derived
  from any prior fit; adjust if the optimizer struggles to converge).
- `--max-nfev` -- cap on objective evaluations for the full-pressure-set
  stage (default 30).
- `--bracket-max-nfev` -- cap on objective evaluations for the initial
  lowest+highest-pressure bracket stage (default 15); `--no-bracket-stage`
  skips straight to the full-pressure fit.
- `--min-time` -- seconds of simulated time to exclude from the start of
  each run before measuring the regression rate (default 0.0, no
  exclusion). At 3 MPa the front position is still transient for roughly
  the first half of the (now 8.0e-2 s) run, so pass e.g. `--min-time 0.04`
  to keep that startup transient out of the steady-state window
  `regression_rate.py` fits over.
- `--xtol` -- `least_squares` convergence tolerance (default 1e-3).
- `--lowmach-bin` / `--template` -- override paths if they differ on the
  machine running this (defaults assume
  `/home/mungerct/research/alamo/bin/lowmach-2d-hdf5-clang++` and the
  template next to this file).

### Before running for real: build and smoke-test

This machine's `bin/lowmach-2d-hdf5-clang++` may or may not still reflect
the current source tree (AMReX headers/libs live outside the repo and
aren't tracked in this worktree). Before the real optimization:

1. Build the `lowmach` target (see the notes in this repo's commit history
   / `FULLFEEDBACK_CALIBRATION.md` for the include/lib path fixes this
   build has needed before -- MPI, HDF5, a matching `CC`/`COMP` for the
   prebuilt AMReX library).
2. Run **one** rendered input directly (not through the sweep script) with
   a short `stop_time` override, e.g.:
   ```bash
   python scripts/render_htpb_sandwich_input.py \
       --template input.lm.ap_htpb_fullfeedback.template \
       --pressure-mpa 1.5 --htpb-pre-exponential 0.001 \
       --htpb-activation-temperature 4000 --out /tmp/smoke_input
   mpirun -np 2 ./bin/lowmach-2d-hdf5-clang++ /tmp/smoke_input \
       stop_time=2.0e-4_s plot_file=/tmp/smoke_output
   python scripts/regression_rate.py /tmp/smoke_output
   ```
   Confirm it runs without crashing, both AP_solid and HTPB_solid regress
   (not just AP -- check the plotfile has moving fronts on both sides of
   the sandwich), and `rigid_eta` produces a sensible front position. This
   case is 2D with a corrugated regressing surface (unlike the 1D-like AP
   monopropellant case), so also sanity-check the `--field rigid_eta` front
   position is well-defined (a single crossing per x-column) at a few
   pressures before trusting the automated sweep.
3. Only then run the real multi-hour optimization.

### Cost

Each sim is a full 2D run (`amr.max_level = 2`, `stop_time = 8.0e-2_s`,
`plot_dt = 1.0e-3_s` in the template -- raised from the original 3.0e-2 s/
5.0e-4 s after the 3 MPa case was found to still be in its startup
transient for roughly the first half of a 3.0e-2 s run; use `--min-time` to
exclude that transient from the regression-rate fit) and will be substantially more
expensive per-run than the 1D-like AP monopropellant sweep. With only 4 fit
pressures (vs. AP's 6), each `least_squares` iteration runs 4 sims in
parallel, but expect this to still be a multi-hour-to-multi-day background
job depending on the machine -- run it under `nohup`/`tmux`/similar.
Progress streams to stdout and to `<workdir>/iterations.jsonl` as it goes.

### Outputs (in `--workdir`)

- `iterations.jsonl` -- one JSON record per objective evaluation: iteration
  number, both HTPB parameters, per-pressure sim rates, residual norm.
- `best_fit.json` -- final `htpb_pre_exponential` / `htpb_activation_temperature`.
- `input.lm.ap_htpb_fullfeedback` -- calibrated input rendered at 3 MPa,
  ready to run directly.
- `validation/` -- the full 4-pressure sweep at the best-fit parameters.
- `htpb_fullfeedback_fit.png` -- simulated vs. experimental r(P) overlay
  (log-log, matching the style of the original Fig. 3).

## Verification already done (this session)

- The digitization pipeline (pixel-precise axis calibration, marker
  centroiding, overlap splitting) was applied to Fig. 3 the same way it was
  verified against a cross-check in the earlier Chakravarthy CSV (duplicate
  series appearing in two different figures agreeing to <2.5%); here, the
  13 raw digitized points were cross-checked against Table 1's independently
  reported numeric values (e.g. binder-width/pressure/rate triples like
  "160 um, 0.2 MPa, 1.11 mm/s") and agreed well.
- The template renders cleanly: `render_htpb_sandwich_input.py` was run
  directly (pure text substitution, no simulation) and confirmed to leave no
  unreplaced `@...@` placeholders, and to correctly hardcode AP's fixed
  fullfeedback parameters while substituting HTPB's.
- `render_htpb_sandwich_input.py`, `optimize_htpb_fullfeedback.py` pass
  `python3 -m py_compile`; `run_htpb_sandwich_pressure_sweep.sh` passes
  `bash -n`.
- **Not done** (Chorpening-geometry template, `input.lm.ap_htpb_fullfeedback.template`):
  the smoke test in the "Before running for real" section above has not
  been executed against this exact template, so there is no confirmation
  yet that `stop_time`/mesh choices are adequate for its ~200 um lamina. Do
  that first, following the same pattern used below for the 100 um deck.
- **Done** (100 um-lamina template, `input.lm.ap_htpb_fullfeedback_100um.template`):
  built `lowmach-2d-hdf5-clang++` with real AMReX HDF5 support (see
  `FULLFEEDBACK_CALIBRATION.md`'s build notes) and ran two short smoke
  tests at 1.5 MPa with the untuned initial guess
  (`htpb_pre_exponential=0.001`, `htpb_activation_temperature=4000`): a
  truncated `stop_time=2.0e-4_s` run to confirm no crash, then a longer
  `stop_time=2.0e-3_s`/`amr.plot_dt=2.5e-4_s` run to check front motion.
  Both completed with a stable timestep (no dynamictimestep collapse). The
  combined `rigid_eta` front reached a "steady" classification from
  `regression_rate.py` at 2.39 mm/s -- the right order of magnitude versus
  both HTPB datasets' ~2.1-2.7 mm/s near 1.5 MPa, for parameters that were
  never fit. Not yet run through the real optimizer, and per-species
  (AP-only vs. HTPB-only) front motion was not independently verified
  beyond the combined rate -- `regression_rate.py` only tracks the summed
  `rigid_eta` field, not per-species fronts.

## If the linear-pressure model can't match the curve shape

Same contingency as the AP calibration: `fullfeedback.pressure_dependence`
only supports a linear `P/reference_pressure` mobility term. The paper
reports `r_b ~ P^0.4` for the *sandwich as a whole* (not HTPB alone), so a
poor fit here doesn't necessarily mean HTPB's own pressure dependence is
wrong -- the sandwich rate is a nonlinear combination of both species'
local mobilities and the 2D flame-coupling effects the paper's own
discussion attributes to (not modeled by this phase-field approach at all).
If the fit is qualitatively poor, consider first whether `stop_time`/mesh
resolution issues (never checked this session) are the cause before
concluding the pressure-dependence form itself needs revisiting.
