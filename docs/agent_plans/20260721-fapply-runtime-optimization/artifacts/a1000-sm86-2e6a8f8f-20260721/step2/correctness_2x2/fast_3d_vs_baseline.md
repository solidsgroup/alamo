# Physics validation compare report

- reference: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/validation/fast_3d`
- candidate: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/step2/correctness_2x2/fast_3d`
- overall verdict: **PASS** (gate verdict: **PASS**)

## centre_bore_3d_128_a2_converged -- PASS (gate: PASS)

### CORRECTNESS

| observable | A | B | abs Δ | rel Δ | tol | status |
| --- | --- | --- | --- | --- | --- | --- |
| eta_field_l2 | 0.0496029 | 0.0496029 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| eta_field_linf | 1 | 1 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| phi_field_l2 | 0.0414434 | 0.0414434 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| temp_field_l2 | 15.7876 | 15.7876 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| disp_field_l2.x | 1.43148e-07 | 1.43148e-07 | 4.095e-17 | 2.861e-10 | 1.000e-06 | PASS |
| disp_field_l2.y | 1.43148e-07 | 1.43148e-07 | 3.006e-17 | 2.100e-10 | 1.000e-06 | PASS |
| disp_field_l2.z | 2.86443e-08 | 2.86443e-08 | 2.647e-16 | 9.242e-09 | 1.000e-06 | PASS |
| disp_field_linf.x | 1.76884e-05 | 1.76884e-05 | 1.693e-15 | 9.573e-11 | 1.000e-06 | PASS |
| disp_field_linf.y | 1.76884e-05 | 1.76884e-05 | 5.345e-16 | 3.022e-11 | 1.000e-06 | PASS |
| disp_field_linf.z | 4.08954e-06 | 4.08954e-06 | 7.295e-14 | 1.784e-08 | 1.000e-06 | PASS |
| stress_field_l2.xx | 962.999 | 962.999 | 1.891e-07 | 1.963e-10 | 1.000e-06 | PASS |
| stress_field_l2.xy | 567.881 | 567.881 | 3.319e-08 | 5.845e-11 | 1.000e-06 | PASS |
| stress_field_l2.xz | 228.282 | 228.282 | 4.208e-07 | 1.843e-09 | 1.000e-06 | PASS |
| stress_field_l2.yx | 567.881 | 567.881 | 3.321e-08 | 5.847e-11 | 1.000e-06 | PASS |
| stress_field_l2.yy | 962.999 | 962.999 | 5.259e-08 | 5.461e-11 | 1.000e-06 | PASS |
| stress_field_l2.yz | 228.282 | 228.282 | 3.913e-07 | 1.714e-09 | 1.000e-06 | PASS |
| stress_field_l2.zx | 228.296 | 228.296 | 4.201e-07 | 1.840e-09 | 1.000e-06 | PASS |
| stress_field_l2.zy | 228.296 | 228.296 | 3.906e-07 | 1.711e-09 | 1.000e-06 | PASS |
| stress_field_l2.zz | 449.323 | 449.323 | 1.634e-07 | 3.636e-10 | 1.000e-06 | PASS |
| stress_field_linf.xx | 65801.1 | 65801.1 | 1.545e-06 | 2.348e-11 | 1.000e-06 | PASS |
| stress_field_linf.xy | 40288 | 40288 | 1.673e-06 | 4.153e-11 | 1.000e-06 | PASS |
| stress_field_linf.xz | 92476.4 | 92476.4 | 2.930e-03 | 3.168e-08 | 1.000e-06 | PASS |
| stress_field_linf.yx | 40288 | 40288 | 1.674e-06 | 4.155e-11 | 1.000e-06 | PASS |
| stress_field_linf.yy | 65801.1 | 65801.1 | 1.032e-06 | 1.569e-11 | 1.000e-06 | PASS |
| stress_field_linf.yz | 92476.4 | 92476.4 | 2.934e-03 | 3.172e-08 | 1.000e-06 | PASS |
| stress_field_linf.zx | 92406.3 | 92406.3 | 2.925e-03 | 3.165e-08 | 1.000e-06 | PASS |
| stress_field_linf.zy | 92406.3 | 92406.3 | 2.928e-03 | 3.169e-08 | 1.000e-06 | PASS |
| stress_field_linf.zz | 75917.8 | 75917.8 | 1.551e-03 | 2.043e-08 | 1.000e-06 | PASS |
| strain_field_l2.xx | 0.0526181 | 0.0526181 | 2.082e-17 | 3.956e-16 | 1.000e-06 | PASS |
| strain_field_l2.xy | 6.73497e-06 | 6.73497e-06 | 2.507e-15 | 3.723e-10 | 1.000e-06 | PASS |
| strain_field_l2.xz | 7.24055e-06 | 7.24055e-06 | 1.899e-14 | 2.623e-09 | 1.000e-06 | PASS |
| strain_field_l2.yx | 6.73497e-06 | 6.73497e-06 | 2.158e-15 | 3.205e-10 | 1.000e-06 | PASS |
| strain_field_l2.yy | 0.0526181 | 0.0526181 | 1.388e-17 | 2.637e-16 | 1.000e-06 | PASS |
| strain_field_l2.yz | 7.24055e-06 | 7.24055e-06 | 1.874e-14 | 2.588e-09 | 1.000e-06 | PASS |
| strain_field_l2.zx | 1.20377e-06 | 1.20377e-06 | 1.091e-14 | 9.064e-09 | 1.000e-06 | PASS |
| strain_field_l2.zy | 1.20377e-06 | 1.20377e-06 | 1.104e-14 | 9.170e-09 | 1.000e-06 | PASS |
| strain_field_l2.zz | 0.0526182 | 0.0526182 | 2.082e-17 | 3.956e-16 | 1.000e-06 | PASS |
| elastic_residual_at_convergence | 4.45907e-09 | 9.72181e-09 | 5.263e-09 | 5.413e-01 | 1.000e-08 | PASS |

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
| max_von_mises_stress | 2.346e-08 | 3.957e-03 | -- | 0.05 | True | PASS |
| mean_von_mises_stress | 1.805e-10 | 4.255e-06 | -- | 0.05 | True | PASS |
| max_principal_stress | 1.328e-08 | 1.897e-03 | -- | 0.05 | True | PASS |
| elastic_energy | 1.427e-10 | 6.490e-09 | -- | 0.03 | True | PASS |

### SOLVER-HEALTH (non-gating)

| observable | A | B | delta |
| --- | --- | --- | --- |
| newton_iters_per_solve | 2 | 2 | 0 |
| newton_final_relnorm | 2.72291e-06 | 2.72291e-06 | 0.0 |
| mlmg_vcycles_per_solve | 24 | 30 | 6 |
| bottom_solver_iters | {'mean': 5.09469696969697, 'max': 6, 'n': 264} | {'mean': 4.8011695906432745, 'max': 6, 'n': 342} | -0.2935273790536952 |
| residual_monotonic | False | False | 0 |
