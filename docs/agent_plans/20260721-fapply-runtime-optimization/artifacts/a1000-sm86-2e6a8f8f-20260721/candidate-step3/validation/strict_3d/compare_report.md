# Physics validation compare report

- reference: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/validation/strict_3d`
- candidate: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/candidate-step3/validation/strict_3d`
- overall verdict: **PASS** (gate verdict: **PASS**)

## centre_bore_3d_128_a2_converged -- PASS (gate: PASS)

### CORRECTNESS

| observable | A | B | abs Δ | rel Δ | tol | status |
| --- | --- | --- | --- | --- | --- | --- |
| eta_field_l2 | 0.0496029 | 0.0496029 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| eta_field_linf | 1 | 1 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| phi_field_l2 | 0.0414434 | 0.0414434 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| temp_field_l2 | 15.7876 | 15.7876 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| disp_field_l2.x | 1.43148e-07 | 1.43148e-07 | 5.327e-18 | 3.721e-11 | 1.000e-06 | PASS |
| disp_field_l2.y | 1.43148e-07 | 1.43148e-07 | 9.745e-19 | 6.808e-12 | 1.000e-06 | PASS |
| disp_field_l2.z | 2.86443e-08 | 2.86443e-08 | 1.172e-18 | 4.092e-11 | 1.000e-06 | PASS |
| disp_field_linf.x | 1.76884e-05 | 1.76884e-05 | 5.086e-16 | 2.875e-11 | 1.000e-06 | PASS |
| disp_field_linf.y | 1.76884e-05 | 1.76884e-05 | 8.414e-16 | 4.757e-11 | 1.000e-06 | PASS |
| disp_field_linf.z | 4.08954e-06 | 4.08954e-06 | 1.086e-16 | 2.656e-11 | 1.000e-06 | PASS |
| stress_field_l2.xx | 962.999 | 962.999 | 9.390e-08 | 9.751e-11 | 1.000e-06 | PASS |
| stress_field_l2.xy | 567.881 | 567.881 | 2.724e-08 | 4.796e-11 | 1.000e-06 | PASS |
| stress_field_l2.xz | 228.282 | 228.282 | 1.356e-08 | 5.938e-11 | 1.000e-06 | PASS |
| stress_field_l2.yx | 567.881 | 567.881 | 2.725e-08 | 4.798e-11 | 1.000e-06 | PASS |
| stress_field_l2.yy | 962.999 | 962.999 | 8.677e-08 | 9.011e-11 | 1.000e-06 | PASS |
| stress_field_l2.yz | 228.282 | 228.282 | 1.382e-08 | 6.056e-11 | 1.000e-06 | PASS |
| stress_field_l2.zx | 228.296 | 228.296 | 1.356e-08 | 5.939e-11 | 1.000e-06 | PASS |
| stress_field_l2.zy | 228.296 | 228.296 | 1.382e-08 | 6.054e-11 | 1.000e-06 | PASS |
| stress_field_l2.zz | 449.323 | 449.323 | 6.461e-08 | 1.438e-10 | 1.000e-06 | PASS |
| stress_field_linf.xx | 65801.1 | 65801.1 | 5.748e-08 | 8.735e-13 | 1.000e-06 | PASS |
| stress_field_linf.xy | 40288 | 40288 | 8.688e-07 | 2.157e-11 | 1.000e-06 | PASS |
| stress_field_linf.xz | 92476.4 | 92476.4 | 1.763e-06 | 1.907e-11 | 1.000e-06 | PASS |
| stress_field_linf.yx | 40288 | 40288 | 8.683e-07 | 2.155e-11 | 1.000e-06 | PASS |
| stress_field_linf.yy | 65801.1 | 65801.1 | 5.308e-07 | 8.067e-12 | 1.000e-06 | PASS |
| stress_field_linf.yz | 92476.4 | 92476.4 | 2.915e-06 | 3.152e-11 | 1.000e-06 | PASS |
| stress_field_linf.zx | 92406.3 | 92406.3 | 1.758e-06 | 1.902e-11 | 1.000e-06 | PASS |
| stress_field_linf.zy | 92406.3 | 92406.3 | 2.907e-06 | 3.145e-11 | 1.000e-06 | PASS |
| stress_field_linf.zz | 75917.8 | 75917.8 | 5.120e-07 | 6.744e-12 | 1.000e-06 | PASS |
| strain_field_l2.xx | 0.0526181 | 0.0526181 | 4.857e-17 | 9.231e-16 | 1.000e-06 | PASS |
| strain_field_l2.xy | 6.73497e-06 | 6.73497e-06 | 2.153e-16 | 3.197e-11 | 1.000e-06 | PASS |
| strain_field_l2.xz | 7.24055e-06 | 7.24055e-06 | 1.216e-16 | 1.679e-11 | 1.000e-06 | PASS |
| strain_field_l2.yx | 6.73497e-06 | 6.73497e-06 | 9.443e-17 | 1.402e-11 | 1.000e-06 | PASS |
| strain_field_l2.yy | 0.0526181 | 0.0526181 | 1.388e-17 | 2.637e-16 | 1.000e-06 | PASS |
| strain_field_l2.yz | 7.24055e-06 | 7.24055e-06 | 3.093e-17 | 4.271e-12 | 1.000e-06 | PASS |
| strain_field_l2.zx | 1.20377e-06 | 1.20377e-06 | 3.490e-17 | 2.900e-11 | 1.000e-06 | PASS |
| strain_field_l2.zy | 1.20377e-06 | 1.20377e-06 | 7.119e-17 | 5.914e-11 | 1.000e-06 | PASS |
| strain_field_l2.zz | 0.0526182 | 0.0526182 | 6.939e-18 | 1.319e-16 | 1.000e-06 | PASS |
| elastic_residual_at_convergence | 4.45907e-09 | 4.45907e-09 | 2.394e-15 | 5.369e-07 | 1.000e-08 | PASS |

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
| max_von_mises_stress | 1.224e-11 | 2.065e-06 | -- | 0.05 | True | PASS |
| mean_von_mises_stress | 9.449e-11 | 2.228e-06 | -- | 0.05 | True | PASS |
| max_principal_stress | 1.698e-11 | 2.426e-06 | -- | 0.05 | True | PASS |
| elastic_energy | 1.507e-10 | 6.855e-09 | -- | 0.03 | True | PASS |

### SOLVER-HEALTH (non-gating)

| observable | A | B | delta |
| --- | --- | --- | --- |
| newton_iters_per_solve | 2 | 2 | 0 |
| newton_final_relnorm | 2.72291e-06 | 2.72291e-06 | 0.0 |
| mlmg_vcycles_per_solve | 24 | 24 | 0 |
| bottom_solver_iters | {'mean': 5.09469696969697, 'max': 6, 'n': 264} | {'mean': 5.09469696969697, 'max': 6, 'n': 264} | 0.0 |
| residual_monotonic | False | False | 0 |
