# Planar Allen–Cahn rate calibration

This calibration isolates the existing rigid-solid `phase_change` operator in a
narrow periodic-x LowMach strip.  The production input
`tests/LMRFMonoAP/input` is read-only; the commands in
`run_calibration.sh` override the heat flux, conduction, temperature advection,
Rocfire rates, activation/cutoff, and mesh.  The interface is the planar
`eta=1 -> 0` AP-solid profile centered at `y=0`; `eta=0.5` moves toward
decreasing `y`.  `lmrf_interface.py` fits the last half of each history.

## Continuum mapping

For `w0=0`, `w12=2`, `w1=1`, `AllenCahn.H` gives

```text
w'(eta) = 96 eta (1-eta) (9/16 - eta).
```

With constant temperature and Arrhenius factor
`f(T)=exp(-activation_temperature/T)` above cutoff, the exact traveling
profile is

```text
eta(z) = 1/(1 + exp(z/delta)),
delta = sqrt(kappa/(48 lambda)),  z = y-c t,
|c| = 6 M rate_multiplier f lambda delta
     = (sqrt(3)/2) M rate_multiplier f sqrt(lambda kappa).
```

Equivalently, using the tanh width `ell=4 delta`,
`|c| = 1.5 M rate_multiplier f lambda ell`.  For the calibration values

```text
lambda = 8.333333333 Pa, kappa = 1e-8 J/m, M = 0.01 1/(Pa s),
delta = 5 um, ell = 20 um,
```

the coefficient is `2.5e-6 m/s` (`0.0025 mm/s`) per unit
`rate_multiplier` when `f=1`.  Thus

```text
rate_multiplier = v_target / (2.5e-6 m/s * f(T)).
```

For Chen-style Arrhenius terms, use `activation_temperature=7500 K` for the
binder and `11000 K` for AP, with the corresponding `f(T)` in this formula.
The calibration runs below intentionally set activation and cutoff to zero so
that mesh error is measured independently of thermal kinetics.

## Measured speeds

The test uses `rate_multiplier=1e5`, so the exact isothermal prediction is
`250.000 mm/s`.  `dx` is the y-cell size; all runs have one periodic x strip.

| T (K) | Ny | dx (um) | measured (mm/s) | relative to exact | median 10–90 width (um) |
|---:|---:|---:|---:|---:|---:|
| 300 | 40 | 5.0 | 240.687 | -3.73% | 22.22 |
| 300 | 80 | 2.5 | 246.367 | -1.45% | 22.02 |
| 600 | 40 | 5.0 | 241.848 | -3.26% | 22.16 |
| 600 | 80 | 2.5 | 247.423 | -1.03% | 21.98 |

The coarse-to-fine change is 2.36% (300 K) and 2.31% (600 K).  Assuming
second-order spatial error, two-level Richardson estimates are 248.26 and
249.28 mm/s, respectively (within 0.70% and 0.29% of 250 mm/s).  The two
temperatures agree within 0.43% on the fine mesh, as expected when the
Arrhenius factor is set to one.  Width remains close to the analytic 10–90
value `2 ln(9) delta = 21.97 um`.

Measured JSON and histories are in the four `runs/T*_N*_R1e5` directories.
The representative wall times were approximately 0.4 s (300 K, Ny=40),
1.9 s (300 K, Ny=80), 1.5 s (600 K, Ny=40), and 4.5 s (600 K, Ny=80) on the
single-process current binary; dynamic timestepping is much more restrictive
at 600 K.

## Reproduction

From the repository root:

```bash
bash validation_homogeneous_chen2002/calibration/run_calibration.sh
bash validation_homogeneous_chen2002/calibration/analyze_calibration.sh
```

The key per-case command is, for example:

```bash
bin/lowmach-2d-clang++ tests/LMRFMonoAP/input \
  plot_file=validation_homogeneous_chen2002/calibration/runs/T300_N80_R1e5 \
  amr.n_cell="2 80" amr.max_grid_size=80 stop_time=1.5e-4_s \
  amr.plot_dt=1.0e-5_s chemistry.model.rocfire.A="0.0 0.0 0.0 0.0" \
  AP_decomposition.phase_change.rate_multiplier=1.0e5 \
  AP_decomposition.phase_change.activation_temperature=0.0 \
  AP_decomposition.phase_change.temperature_cutoff=0.0 \
  heat_source.ic.expression.constant.qflux=0.0 include_conduction=0 \
  advect_temperature=0
```
