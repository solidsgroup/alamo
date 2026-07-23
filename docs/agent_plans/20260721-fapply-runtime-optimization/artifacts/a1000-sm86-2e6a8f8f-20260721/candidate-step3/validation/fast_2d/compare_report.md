# Physics validation compare report

- reference: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/validation/fast_2d`
- candidate: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/candidate-step3/validation/fast_2d`
- overall verdict: **PASS** (gate verdict: **PASS**)

## canonical_2d_elastic -- PASS (gate: PASS)

### CORRECTNESS

| observable | A | B | abs Δ | rel Δ | tol | status |
| --- | --- | --- | --- | --- | --- | --- |
| eta_field_l2 | 0.171856 | 0.171856 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| eta_field_linf | 1 | 1 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| phi_field_l2 | 0.143237 | 0.143237 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| temp_field_l2 | 53.477 | 53.477 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| disp_field_l2.x | 7.6502e-07 | 7.6502e-07 | 1.874e-20 | 2.450e-14 | 1.000e-06 | PASS |
| disp_field_l2.y | 7.48673e-07 | 7.48673e-07 | 9.211e-21 | 1.230e-14 | 1.000e-06 | PASS |
| disp_field_linf.x | 1.86033e-05 | 1.86033e-05 | 2.575e-19 | 1.384e-14 | 1.000e-06 | PASS |
| disp_field_linf.y | 1.86244e-05 | 1.86244e-05 | 1.355e-20 | 7.277e-16 | 1.000e-06 | PASS |
| stress_field_l2.xx | 8279.5 | 8279.5 | 2.328e-10 | 2.812e-14 | 1.000e-06 | PASS |
| stress_field_l2.xy | 2956.11 | 2956.11 | 6.366e-12 | 2.154e-15 | 1.000e-06 | PASS |
| stress_field_l2.yx | 2956.13 | 2956.13 | 6.366e-12 | 2.154e-15 | 1.000e-06 | PASS |
| stress_field_l2.yy | 8056.88 | 8056.88 | 1.728e-11 | 2.145e-15 | 1.000e-06 | PASS |
| stress_field_linf.xx | 198551 | 198551 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| stress_field_linf.xy | 85036.3 | 85036.3 | 2.357e-09 | 2.772e-14 | 1.000e-06 | PASS |
| stress_field_linf.yx | 85025.4 | 85025.4 | 2.314e-09 | 2.721e-14 | 1.000e-06 | PASS |
| stress_field_linf.yy | 173075 | 173075 | 2.264e-08 | 1.308e-13 | 1.000e-06 | PASS |
| strain_field_l2.xx | 0.178213 | 0.178213 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| strain_field_l2.xy | 6.07883e-05 | 6.07883e-05 | 9.216e-19 | 1.516e-14 | 1.000e-06 | PASS |
| strain_field_l2.yx | 6.28411e-05 | 6.28411e-05 | 2.724e-18 | 4.335e-14 | 1.000e-06 | PASS |
| strain_field_l2.yy | 0.178214 | 0.178214 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| elastic_residual_at_convergence | 6.89591e-09 | 6.89625e-09 | 3.447e-13 | 4.999e-05 | 1.000e-08 | PASS |

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
| max_von_mises_stress | 1.081e-15 | 2.037e-10 | -- | 0.05 | True | PASS |
| mean_von_mises_stress | 1.624e-14 | 8.877e-10 | -- | 0.05 | True | PASS |
| max_principal_stress | 8.785e-16 | 1.746e-10 | -- | 0.05 | True | PASS |
| elastic_energy | 2.045e-14 | 2.274e-11 | -- | 0.03 | True | PASS |

### SOLVER-HEALTH (non-gating)

| observable | A | B | delta |
| --- | --- | --- | --- |
| newton_iters_per_solve | 2 | 2 | 0 |
| newton_final_relnorm | 4.95084e-07 | 4.95084e-07 | 0.0 |
| mlmg_vcycles_per_solve | 23 | 23 | 0 |
| bottom_solver_iters | {'mean': 1.0, 'max': 1, 'n': 246} | {'mean': 1.0, 'max': 1, 'n': 246} | 0.0 |
| residual_monotonic | False | False | 0 |
