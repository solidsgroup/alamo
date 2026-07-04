# Resyncing MLMG coefficients after Newton relinearization

## Problem

`Solver::Nonlocal::Newton` reuses one persistent `amrex::MLMG` object across
Newton iterations. `prepareForSolve()` updates the fine-level elastic
coefficient field through `SetModel()`, but AMReX calls the operator's
`prepareForSolve()` only on the first MLMG solve. That leaves coarse MG-level
operators and the smoother/normalize diagonal stale after Newton
relinearization.

At high contrast this caused the later Newton linear solves to smooth the new
fine operator against an old MG hierarchy, producing explosive divergence.

## Source change

`Operator::Operator<Grid::Node>` now exposes:

```cpp
void SyncCoefficients() { averageDownCoeffs(); Diagonal(true); }
```

`Solver::Nonlocal::Newton` calls it after `prepareForSolve()` on Newton
iterations after the first when `elastic.solver.resync_coeffs=1`.

The option is default-off, so historical input decks keep their old trajectory
unless they opt in.

## Evidence

Anchor, `psi_floor=0`, same binary/config except the gate:

| run | result |
|---|---|
| `validation_20260704/codex_pf0_resync1_20260704_123913.out` | PASS. Newton linear solves converged in 240/248/228 iterations. |
| `validation_20260704/codex_pf0_resync0_20260704_123929.out` | Control FAIL. Newton iteration 2 reproduced stale-hierarchy divergence after 9 MLMG iterations, `resid/resid0 ~= 2.56e22`. |

Additional checks:

- `validation_20260704/codex_gpu3d_resync_parse_smoke_20260704_124158.out`
  parsed and finalized on strict 3D CUDA.
- `validation_20260704/codex_gpu3d_resync_step_smoke_20260704_124217.out`
  entered two Newton/MLMG iterations on strict 3D CUDA and finalized.

## Use

Use:

```text
elastic.solver.resync_coeffs=1
```

for high-contrast Newton solves, especially low `psi_floor` chamber cases.
