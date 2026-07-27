# Calibrating the `fullfeedback` AP regression rate against experimental r(P)

Fits **one** `(pre_exponential, activation_temperature)` pair for
`Model::PhaseField::FullFeedback` to the whole experimental AP
regression-rate-vs-pressure curve in `AP_reg_rate.csv` (2.76-13.79 MPa), not a
single pressure point. `pressure_dependence` stays linear
(`P / reference_pressure`); no C++ code change is involved.

## Pieces

- `input.lm.ap_monopropellant_fullfeedback.template` — the monopropellant
  input parameterized by operating pressure and the two fit parameters.
  Placeholders are filled in by `scripts/render_input.py`, which also derives
  the pressure-dependent product density the same way the original
  `input.lm.ap_monopropellant_fullfeedback` does
  (`rho = P / (319.787 * 700)`).
- `scripts/run_pressure_sweep.sh` — for one candidate parameter pair, renders
  and runs one sim per requested pressure **in parallel**, measures each with
  `scripts/regression_rate.py --rate-only --unit cm/s`, and writes a
  `pressure_mpa,reg_rate_cm_s` CSV.
- `scripts/regression_rate.py --rate-only --unit cm/s <plotfile_dir>` — prints
  a single steady-state regression-rate number (or `nan` + nonzero exit if the
  run never reached a steady-burning window).
- `scripts/optimize_fullfeedback.py` — `scipy.optimize.least_squares` driver
  over `log10(pre_exponential)` and `activation_temperature`, with relative
  residuals against a representative pressure subset. Logs every iteration,
  then runs a full 16-pressure validation sweep and writes the calibrated
  input plus a sim-vs-experiment plot.

## Running the real calibration

```bash
python scripts/optimize_fullfeedback.py --workdir /path/to/calib_run
```

Useful flags:
- `--fit-pressures` — pressures (MPa) used during optimization (default:
  `2.85 4.14 5.61 8.27 10.44 13.79`, a 6-point subset spanning the range; all
  16 experimental pressures are used only for the final validation sweep).
- `--pre-exponential0` / `--activation-temperature0` — initial guess
  (defaults `0.0027` 1/Pa/s, `3145` K).
- `--max-nfev` — cap on objective evaluations for the full-pressure-set stage
  (default 30).
- `--no-bracket-stage` / `--bracket-max-nfev` — by default, optimization
  first fits against only the lowest and highest `--fit-pressures` (cheapest
  possible sweep), then refines against the full set from that result;
  `--no-bracket-stage` skips straight to the full-pressure fit.
- `--min-time-low` / `--min-time-high` — seconds of simulated time to
  exclude from the start of each run before measuring the regression rate,
  linearly interpolated per-pressure between the lowest and highest pressure
  in each sweep (default 0.0/0.0, no exclusion). The front position is still
  transient for roughly the first half of the run at low pressure, so e.g.
  `--min-time-low 0.04 --min-time-high 0.01` keeps that startup transient
  out of the steady-state window.
- `--xtol` — `least_squares` convergence tolerance (default 1e-3).
- `--lowmach-bin` / `--template` — override paths if they differ on the
  machine running this (defaults assume
  `/home/mungerct/research/alamo/bin/lowmach-2d-hdf5-clang++` and the
  template next to this file).

### Cost

Each sim is `stop_time = 8.0e-2_s`, `plot_dt = 5.0e-4_s` for enough samples
per run so the steady-state window is well resolved by
`regression_rate.py`'s transient/steady/extinguished classification (after
`--min-time` excludes the startup transient). Each
`least_squares` iteration runs the full fit-pressure subset in parallel, so
the whole optimization is a **multi-hour background job** — run it under
`nohup`/`tmux`/similar. Progress streams to stdout and to
`<workdir>/iterations.jsonl` as it goes, so it's safe to monitor or resume
from partial output.

### Outputs (in `--workdir`)

- `iterations.jsonl` — one JSON record per objective evaluation: iteration
  number, both parameters, per-pressure sim rates, residual norm.
- `best_fit.json` — final `pre_exponential` / `activation_temperature`.
- `input.lm.ap_monopropellant_fullfeedback` — calibrated input rendered at
  3 MPa, ready to run directly.
- `validation/` — the full 16-pressure sweep at the best-fit parameters.
- `fullfeedback_fit.png` — simulated vs. experimental r(P) overlay.

## Verification already done

The full chain (render -> parallel sims -> rate extraction -> `least_squares`
iteration -> validation sweep -> plot) was smoke-tested with an artificially
truncated `stop_time` template to confirm there are no crashes and the CSV/
JSON/PNG outputs are well-formed, including when a run doesn't reach a steady
state (handled as a penalty residual, not a crash). The actual multi-hour
optimization has **not** been run to completion — that's the point of this
runbook.

## If the linear-pressure model can't match the curve shape

`fullfeedback.pressure_dependence` only supports a linear `P/reference_pressure`
mobility term. If the best fit from this pipeline is clearly too far off in
*shape* (not just magnitude), the documented contingency is adding a
`pressure_exponent` power-law term to `src/Model/PhaseField/FullFeedback.H`
and re-fitting `(pre_exponential, activation_temperature, pressure_exponent)`
— out of scope for this pipeline as built.
