# Current high-contrast elastic settings

Use this recipe for the best measured performance/stability balance from the
2026-07-04 anchor campaign:

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

Notes:

- `nr_convergence=psi_update` is the performance win. On the anchor first
  elastic window it stopped at 13 Newton solves, about 27.5 s, with
  `psi_update ~= 5.7e-7`. The raw-update reference took 178 Newton solves and
  about 215 s.
- `nrtolerance=1e-6` and `1e-5` stopped at the same point in the measured
  `psi_update` run. Use `1e-6` because it expresses the observed threshold more
  honestly without adding cost in the current case.
- `nriters=5` is the measured natural Newton cap for this recipe. The 8.0 s
  first-window performance screen used `nriters=3`; that was a fixed-budget
  stability screen, not the natural psi-update convergence recipe.
- `nrtolerance=5e-7` is a modest-accuracy point: 22 Newton solves, about
  39.0 s, lower `stress_xx` error than the 13-solve point.
- `nrtolerance=1e-7` is an accuracy point, not a performance point: 128 Newton
  solves and about 164 s.
- `zero_out_displacement=1` should stay on for anchored chamber runs unless the
  experiment is explicitly about warm starts.

For the rod-and-tube thin-seam class, add the validated input-only linear
robustness knobs:

```text
elastic.max_coarsening_level=2
elastic.solver.bottom_solver=smoother
```

Avoid these settings as production recommendations:

- `elastic.psi_floor=0.005` with void `4/3 MPa`: still fails during the fourth
  Newton linear solve.
- Diagonal inflation: refuted and removed.
- Quick AMReX `GMRES_MLMG` wrapping: crashed before useful solver status.
- Residual-only or stagnation-gated Newton exits: not a better cost/accuracy
  point than the simple psi-update sweep.
