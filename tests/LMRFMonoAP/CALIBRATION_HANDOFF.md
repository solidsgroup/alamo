# AP monopropellant regression-rate calibration — handoff

Branch: `eulerian-solids-heat-conduction-laser`. Goal: fit `tests/LMRFMonoAP`
regression rate to `reference.csv` (2-6 MPa) by tuning free params in
`scripts/optimize_ap_regression.py` (`rate_multiplier`, `activation_temperature`,
`w1`, `pressure_exponent`).

## Current best fit (DE iter 127, residual_norm=0.3374)

```
rate_multiplier        = 1545.64
activation_temperature = 610.86 K
w1                      = 1.3797648840875154
pressure_exponent      = 1.9673263657521716
```

| P (MPa) | sim (mm/s) | exp (mm/s) | error |
|---|---|---|---|
| 2 | 2.123 | 2.65 | -19.9% |
| 3 | 3.234 | 3.80 | -14.9% |
| 4 | 4.503 | 5.00 | -9.9% |
| 5 | 5.726 | 6.30 | -9.1% |
| 6 | 6.446 | 7.90 | -18.4% |

Systematic under-prediction at every pressure, U-shaped (worst at low/high P,
best in the middle). Not yet closed.

To resume optimizing: `python3 scripts/optimize_ap_regression.py --workdir <dir>`
(runs stage-1 random search + stage-2 differential evolution; results append
to `<dir>/iterations.jsonl`, one JSON record per eval — best-so-far is
`min(residual_norm)`). Last run was pid 1904695 in
`/home/mungerct/.claude/jobs/05931ebc/tmp/ap_calib_trimmed2` on the original
machine — not portable, restart fresh on the new machine.

## Dead ends already tested (don't re-try without new evidence)

- **`AllenCahn.w12`** (barrier height at η=0.5, fixed at 2.0 in the input,
  not exposed to the optimizer): tested 0.5-8.0 at fixed best-fit chemistry.
  Non-monotonic — peaks near w12≈2.5-3.0 at only ~4.6 mm/s vs 4.5 at w12=2.0
  (~2.5% gain, not enough). w12=4.0 destabilized the solver at 6 MPa (dt
  collapsed to 1e-12, MLMG abort). Leave w12 fixed at 2.0.
- **`amr.max_level` 0→1** (interface resolution, currently fixed at 0 per
  uncommitted `input` change): tested at fixed best-fit chemistry, 4/6 MPa.
  Rate moved *away* from target (4 MPa: 4.516→4.362, 6 MPa: 6.444→6.333), not
  toward it. Doesn't refit the chemistry at the finer resolution though, so
  not a fully fair comparison — but rules out "just bump resolution" as a
  quick win.

## Recommended next step (not yet done)

Check whether the experimental data (`reference.csv`) is even well-described
by a single power law `rate = A * P^n` — curve-fit that directly against the
5 experimental points (no simulation needed). The U-shaped residual pattern
above is a signature of wrong functional form, not just bad parameter values,
so if a single exponent doesn't fit the data's curvature, the model needs a
structural change (e.g. two-regime pressure dependence) rather than more
tuning of the existing 4 parameters.

## Other open items (from earlier in this branch, unrelated to calibration)

- Uncommitted changes: `scripts/optimize_ap_regression.py`,
  `src/Model/Mechanism/PhaseChange.H`, `tests/LMRFMonoAP/input`,
  `tests/LMRFMonoAP/test`. Not staged/committed — do that once fit is settled.
- Once satisfied with the fit, restore `amr.max_level` 0→1 in `input` and
  re-run `./scripts/runtests.py tests/LMRFMonoAP --serial --dim=2` at that
  resolution (see `/home/mungerct/.claude/plans/i-just-merged-in-mutable-lightning.md`
  for the full context on why max_level was dropped to 0 and the conduction-fix
  background).
