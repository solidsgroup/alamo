# Psi-weighted Newton convergence

## Problem

The historical Newton convergence gate uses the maximum displacement correction
over all degrees of freedom. In void/solid problems, low-stiffness void
displacement can keep drifting after the stress-bearing solid region has
settled, so the raw max-update gate can burn many Newton/MLMG solves without a
useful field improvement.

## Source change

`Solver::Nonlocal::Newton` now supports:

```text
elastic.solver.nr_convergence=update      # historical default
elastic.solver.nr_convergence=psi_update  # opt-in
```

`psi_update` computes the maximum accepted displacement update weighted by
local nodal psi. Supporting knobs:

```text
elastic.solver.nr_psi_weight_floor=0
elastic.solver.nr_psi_power=1
elastic.solver.nr_solid_psi=0.5
elastic.solver.nr_diagnostics=0
```

The default mode remains `update`.

## Evidence

Anchor first-elastic-window setup:

```text
elastic.psi_floor=0.01
model_void.kappa=8_MPa
model_void.mu=6_MPa
elastic.zero_out_displacement=1
elastic.solver.line_search=1
elastic.solver.resync_coeffs=1
elastic.solver.nriters=5
```

| run | Newton solves | wall | note |
|---|---:|---:|---|
| `codex_nrdiag_pf010_update_nrtol_20260704_134300.out` | 178 | 215.10 s | Raw-update reference. |
| `codex_nrpsi_pf010_20260704_134719.out` | 13 | 27.42 s | `psi_update ~= 5.7e-7`; best performance point. |
| `codex_nrpsi_pf010_ntol5em7_20260704.out` | 22 | 38.96 s | Modest accuracy improvement. |
| `codex_nrpsi_pf010_ntol1e7_20260704.out` | 128 | 163.67 s | Accuracy point, not performance. |

Lower void modulus at the same floor also passed:

| run | result |
|---|---|
| `codex_nrpsi_pf010_vk4_vm3_20260704_135751.out` | PASS, 13 Newton solves, 27.54 s. |

## Use

For performance, use:

```text
elastic.solver.nr_convergence=psi_update
elastic.solver.nrtolerance=1e-6
```

For a modest accuracy tradeoff:

```text
elastic.solver.nr_convergence=psi_update
elastic.solver.nrtolerance=5e-7
```
