# Negative results and removed experiment paths

This file records what not to pursue from the 2026-07-02 to 2026-07-04
campaign.

## Diagonal inflation

The first model predicted that adding first-order interface coupling magnitude
to the smoother diagonal would stabilize high-contrast interfaces. Production
tests refuted that model:

- Tiny 3D CPU zero-floor case: `diag_inflate=0` converged; all positive factors
  tested (`0.25, 0.5, 1, 2, 4, 8, 16`) failed the 80-iteration probe and got
  worse as the factor increased.
- Bounded 3D GPU floor screen: `diag_inflate=0,line_search=1` completed floors
  down to 0; `diag_inflate=1` timed out or failed for most low floors.

Decision: do not commit or use diagonal inflation. The source knob was removed.

## Residual and stagnation Newton exits

Residual-only convergence accepted too early for the anchor case:

- `codex_nrresid_pf010_r008_diagfix_20260704_140419.out` stopped at iteration
  12 with `psi_update ~= 5.1e-5`, before the psi-weighted update had collapsed.

The stagnation gate worked mechanically, but was not a better cost/accuracy
point:

- `codex_nrpsi_stag_pf010_w8_20260704.out`: 20 Newton solves, about 32.7 s,
  `psi_update=5.14982e-7`, solid `stress_xx` rel max/RMS
  `3.895300e-2 / 1.858812e-3`.
- The simple `psi_update` `nrtolerance=5e-7` run took 22 solves, about 39.0 s,
  and produced slightly better stress comparison.

Decision: keep the code surface small. Only `update` and `psi_update` are
committed convergence modes.

## Lower-floor linear robustness

Target case:

```text
elastic.psi_floor=0.005
model_void.kappa=4_MPa
model_void.mu=3_MPa
elastic.solver.nr_convergence=psi_update
```

Results:

| run | result |
|---|---|
| `codex_linrob_pf005_default_20260704.out` | FAIL, fourth Newton linear solve blows up, `resid/resid0 ~ 1.8e20` by iteration 17. |
| `codex_linrob_pf005_smoother_20260704.out` | FAIL, same fourth-solve blow-up. |
| `codex_linrob_pf005_mcl2_smoother_20260704.out` | FAIL, same fourth-solve blow-up. |
| `codex_linrob_pf005_mcl1_smoother_20260704.out` | Stable first solve, but hits 400-iteration cap at `resid/bnorm=2.5810282e-7`. |
| `codex_linrob_pf005_mcl1_smoother_mi900_20260704.out` | Stable but too slow; first two solves take 611/780 iterations and the third hits 900. |
| `codex_linrob_pf005_mcl1_smoother_tol1e6_20260704.out` | First three solves converge in 336/404/502 iterations, then the fourth solve blows up. |
| `codex_linrob_pf005_mcl1_smoother_tol1e7_20260704.out` | First three solves converge in 466/573/739 iterations, then the fourth solve blows up. |

Decision: `psi_floor=0.005` is not a validated production floor from this
campaign. Linear robustness knobs do not solve it.

## GMRES_MLMG wrapper

A quick default-off AMReX `GMRES_MLMG` wrapper was prototyped, but both probes
segfaulted before any useful GMRES iteration/status:

- `codex_linrob_pf005_gmres_default_tol1e6_20260704.out`
- `codex_linrob_pf005_gmres_default_tol1e6_growfix_20260704.out`

Decision: do not commit a crashing opt-in path. GMRES would need a custom
adapter around this nodal ghost-node operator before it is worth testing again.
