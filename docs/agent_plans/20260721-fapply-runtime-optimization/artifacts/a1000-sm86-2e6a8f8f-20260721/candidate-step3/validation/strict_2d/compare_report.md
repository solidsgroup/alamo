# Physics validation compare report

- reference: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/validation/strict_2d`
- candidate: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/candidate-step3/validation/strict_2d`
- overall verdict: **PASS** (gate verdict: **PASS**)

## canonical_2d_elastic -- PASS (gate: PASS)

### CORRECTNESS

| observable | A | B | abs Δ | rel Δ | tol | status |
| --- | --- | --- | --- | --- | --- | --- |
| eta_field_l2 | 0.171856 | 0.171856 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| eta_field_linf | 1 | 1 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| phi_field_l2 | 0.143237 | 0.143237 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| temp_field_l2 | 53.477 | 53.477 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| disp_field_l2.x | 7.6502e-07 | 7.6502e-07 | 4.235e-22 | 5.536e-16 | 1.000e-06 | PASS |
| disp_field_l2.y | 7.48673e-07 | 7.48673e-07 | 1.059e-22 | 1.414e-16 | 1.000e-06 | PASS |
| disp_field_linf.x | 1.86033e-05 | 1.86033e-05 | 3.388e-21 | 1.821e-16 | 1.000e-06 | PASS |
| disp_field_linf.y | 1.86244e-05 | 1.86244e-05 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| stress_field_l2.xx | 8279.5 | 8279.5 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| stress_field_l2.xy | 2956.11 | 2956.11 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| stress_field_l2.yx | 2956.13 | 2956.13 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| stress_field_l2.yy | 8056.88 | 8056.88 | 3.638e-12 | 4.515e-16 | 1.000e-06 | PASS |
| stress_field_linf.xx | 198551 | 198551 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| stress_field_linf.xy | 85036.3 | 85036.3 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| stress_field_linf.yx | 85025.4 | 85025.4 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| stress_field_linf.yy | 173075 | 173075 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| strain_field_l2.xx | 0.178213 | 0.178213 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| strain_field_l2.xy | 6.07883e-05 | 6.07883e-05 | 6.776e-21 | 1.115e-16 | 1.000e-06 | PASS |
| strain_field_l2.yx | 6.28411e-05 | 6.28411e-05 | 2.711e-20 | 4.313e-16 | 1.000e-06 | PASS |
| strain_field_l2.yy | 0.178214 | 0.178214 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| elastic_residual_at_convergence | 6.89626e-09 | 6.89625e-09 | 4.721e-15 | 6.846e-07 | 1.000e-08 | PASS |

### ENGINEERING-TRAJECTORY

| observable | max rel Δ | max abs Δ | phase-lag | tol | endpoint-only | status |
| --- | --- | --- | --- | --- | --- | --- |
| chamber_pressure | 0.000e+00 | 0.000e+00 | 0 (tol 1) | 0.02 | False | PASS |
| chamber_volume | 0.000e+00 | 0.000e+00 | -- | 0.02 | False | PASS |
| burn_area | 0.000e+00 | 0.000e+00 | -- | 0.02 | False | PASS |
| total_mdot | 0.000e+00 | 0.000e+00 | 0 (tol 1) | 0.02 | False | PASS |
| mdot_max | 0.000e+00 | 0.000e+00 | -- | 0.03 | False | PASS |
| max_temp | 0.000e+00 | 0.000e+00 | -- | 0.02 | False | PASS |
| interface_extent | 0.000e+00 | 0.000e+00 | -- | 0.02 | False | PASS |
| corner_displacement | 0.000e+00 | nan | -- | 0.03 | False | PASS |
| max_von_mises_stress | 0.000e+00 | 0.000e+00 | -- | 0.05 | True | PASS |
| mean_von_mises_stress | 2.662e-16 | 1.455e-11 | -- | 0.05 | True | PASS |
| max_principal_stress | 0.000e+00 | 0.000e+00 | -- | 0.05 | True | PASS |
| elastic_energy | 8.179e-16 | 9.095e-13 | -- | 0.03 | True | PASS |

### SOLVER-HEALTH (non-gating)

| observable | A | B | delta |
| --- | --- | --- | --- |
| newton_iters_per_solve | 2 | 2 | 0 |
| newton_final_relnorm | 4.95084e-07 | 4.95084e-07 | 0.0 |
| mlmg_vcycles_per_solve | 23 | 23 | 0 |
| bottom_solver_iters | {'mean': 1.0, 'max': 1, 'n': 246} | {'mean': 1.0, 'max': 1, 'n': 246} | 0.0 |
| residual_monotonic | False | False | 0 |
