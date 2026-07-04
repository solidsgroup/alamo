# MLMG high-contrast findings

This file is an index for the 2026-07-02 to 2026-07-04 high-contrast elastic
MLMG campaign. The long narrative was split so agents can load only the part
they need.

## Read first

- [Current settings](mlmg_high_contrast_20260702/SETTINGS.md): best measured
  input recipe and caveats.
- [Resync coefficients](mlmg_high_contrast_20260702/RESYNC_COEFFS.md): the
  stale-hierarchy defect and source fix.
- [Psi-weighted Newton convergence](mlmg_high_contrast_20260702/NR_PSI_CONVERGENCE.md):
  the validated fast nonlinear convergence gate.
- [Negative results](mlmg_high_contrast_20260702/NEGATIVE_RESULTS.md):
  refuted paths, including diagonal inflation, residual/stagnation gates,
  lower-floor linear robustness, and quick GMRES wrapping.
- [Evidence bundle](mlmg_high_contrast_20260702/README.md): artifact index and
  run log names.

## Current decision

Commit and use the genuine improvements:

- `elastic.solver.resync_coeffs=1` for high-contrast Newton solves. This
  refreshes coarse MG operators and the smoother/normalize diagonal after each
  Newton relinearization, fixing the stale-hierarchy failure.
- `elastic.solver.nr_convergence=psi_update` for void/solid cases. This gates
  Newton convergence on psi-weighted accepted updates instead of raw max update
  in near-void displacement.
- `elastic.zero_out_displacement=1` for the anchored chamber runs. This avoids
  warm-started void displacement/J-collapse accumulation.
- Keep existing `elastic.solver.line_search=1` enabled.

Do not commit or use the failed experiment knobs:

- Diagonal inflation: refuted and removed from source.
- Residual-only and stagnation-gated Newton exits: not better accuracy/cost
  points than the simple `psi_update` tolerance sweep; removed from source.
- Quick AMReX `GMRES_MLMG` wrapper: crashed before useful iterations/status;
  not committed.

## Best measured point

For the anchor first-elastic-window benchmark:

```text
elastic.psi_floor=0.01
model_void.kappa=8_MPa
model_void.mu=6_MPa
elastic.zero_out_displacement=1
elastic.solver.line_search=1
elastic.solver.resync_coeffs=1
elastic.solver.nr_convergence=psi_update
elastic.solver.nrtolerance=1e-6
elastic.solver.nriters=5
elastic.solver.max_iter=200
```

This stopped at 13 Newton solves in about 27.5 s with
`psi_update ~= 5.7e-7`, versus 178 solves and about 215 s for the raw-update
reference. For a fixed `nriters=3` first-window screen, the same stability recipe
at `psi_floor=0.01` completed in about 8.0 s with MLMG iterations `90/90/113`.

For the known rod-and-tube thin-seam class, add:

```text
elastic.max_coarsening_level=2
elastic.solver.bottom_solver=smoother
```

Do not lower the production floor to `psi_floor=0.005` from this campaign; the
linear robustness screen still fails on the fourth Newton linear solve.
