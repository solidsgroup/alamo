# Elastic GPU Coarsening Sweep

Date: 2026-06-21

Uniform 2048x2048 GPU case using `bin/alamo_gpu-2d-nofast-cuda86-g++`.

Uniform stiffness overrides:

```text
elastic.use_psi=0
model_prop.kappa=162_MPa
model_prop.mu=113.6_MPa
model_void.kappa=162_MPa
model_void.mu=113.6_MPa
model_casing.kappa=162_MPa
model_casing.mu=113.6_MPa
```

| Label | max_coarsening_level | Verdict | AMR levels | MG levels on coarsest AMR level | Last residual | First 10 resid/bnorm values | Notes |
|---|---:|---|---:|---:|---:|---|---|
| `coarsen0_2048_uniform` | 0 | TIMEOUT at 901 s | 1 | 1 | `3.371003866e-05` | `0.09491393843`, `0.02867154622`, `0.01519545109`, `0.009294161261`, `0.004676572806`, `0.003141141538`, `0.001922058033`, `0.001156408331`, `0.0008388655393`, `0.0006438138162` | Did not diverge. Repeated `MLMG: Bottom solve failed.` No `Final Iter` before harness timeout. |
| `coarsen1_2048_uniform_i60` | 1 | DIVERGED after 12 iterations | 1 | 2 | `8.83843826e+27` | `159.2053539`, `16425.38795`, `1891751.941`, `62060851.29`, `1596084623`, `8.38374002e+10`, `3.585730033e+12`, `1.358588362e+14`, `2.660877641e+15`, `8.927711956e+16` | Added `elastic.solver.max_iter=60`. First coarse transition triggers explosive divergence with default `bicgstab` bottom solver. Failure line: `MLMG: Failing to converge after 12 iterations. resid, resid/bnorm = 8.83843826e+27, 2.624074357e+20`. |
| `coarsen1_2048_uniform_i60_transferfix` | 1 | DIVERGED after 13 iterations | 1 | 2 | `8.921138316e+28` | `131.3041612`, `31231.71438`, `629024.5141`, `26176362.86`, `436020717.6`, `2.220075463e+10`, `1.262692919e+12`, `4.660956993e+13`, `2.691044859e+15`, `4.459360065e+16` | Rebuilt with explicit 2D nodal restriction/interpolation branches in `src/Operator/Operator.cpp`. Trajectory changed but default `bicgstab` still explodes. |
| `coarsen1_2048_uniform_i60_bottom_cg` | 1 | FAILED max_iter 60 | 1 | 2 | `3120866.471` | `1.941338281`, `0.4883254257`, `0.3465672508`, `0.3065399361`, `0.2957172738`, `0.2782667489`, `0.315953119`, `0.2414532012`, `0.2136669055`, `0.2960688403` | `elastic.solver.bottom_solver=cg`. No explosive divergence; failed at cap with `resid/bnorm = 0.09265647885`. |
| `coarsen1_2048_uniform_i60_bottom_smoother` | 1 | FAILED max_iter 60 | 1 | 2 | `24822558.1` | `2.88853321`, `2.616710256`, `2.293487952`, `2.311698363`, `2.045740668`, `2.150623731`, `2.158891384`, `1.926880357`, `1.858839395`, `1.680357561` | `elastic.solver.bottom_solver=smoother elastic.solver.final_smooth=16`. No explosive divergence; failed at cap with `resid/bnorm = 0.7369654712`. |

## Commands

```bash
benchmark/elastic_sensitivity_20260621/run_case.sh \
  coarsen0_2048_uniform \
  bin/alamo_gpu-2d-nofast-cuda86-g++ \
  1 \
  'amr.n_cell=2048 2048' \
  elastic.use_psi=0 \
  model_void.kappa=162_MPa \
  model_void.mu=113.6_MPa \
  model_casing.kappa=162_MPa \
  model_casing.mu=113.6_MPa \
  model_prop.kappa=162_MPa \
  model_prop.mu=113.6_MPa \
  elastic.max_coarsening_level=0
```

## Immediate Interpretation

`elastic.max_coarsening_level=0` restricts the solve to one MG level. This run
was slow and timed out, but it showed monotone residual reduction rather than
the documented 2048x2048 GPU divergence.

`elastic.max_coarsening_level=1` creates two MG levels and immediately triggers
explosive divergence with the default `bicgstab` bottom solver. Switching the
bottom solver to `cg` or `smoother` removes the explosive path, though neither
variant converged within the 60-iteration cap. This narrows the immediate target
to the first coarse-grid transition and especially GPU bottom Krylov behavior /
bottom-solve failure handling.
