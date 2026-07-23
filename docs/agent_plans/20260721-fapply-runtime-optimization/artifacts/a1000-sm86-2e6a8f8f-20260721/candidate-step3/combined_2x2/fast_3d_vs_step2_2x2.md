# Physics validation compare report

- reference: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/step2/correctness_2x2/fast_3d`
- candidate: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/candidate-step3/combined_2x2/fast_3d`
- overall verdict: **PASS** (gate verdict: **PASS**)

## centre_bore_3d_128_a2_converged -- PASS (gate: PASS)

### CORRECTNESS

| observable | A | B | abs Δ | rel Δ | tol | status |
| --- | --- | --- | --- | --- | --- | --- |
| eta_field_l2 | 0.0496029 | 0.0496029 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| eta_field_linf | 1 | 1 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| phi_field_l2 | 0.0414434 | 0.0414434 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| temp_field_l2 | 15.7876 | 15.7876 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| disp_field_l2.x | 1.43148e-07 | 1.43148e-07 | 2.027e-18 | 1.416e-11 | 1.000e-06 | PASS |
| disp_field_l2.y | 1.43148e-07 | 1.43148e-07 | 3.188e-18 | 2.227e-11 | 1.000e-06 | PASS |
| disp_field_l2.z | 2.86443e-08 | 2.86443e-08 | 1.356e-18 | 4.734e-11 | 1.000e-06 | PASS |
| disp_field_linf.x | 1.76884e-05 | 1.76884e-05 | 3.260e-17 | 1.843e-12 | 1.000e-06 | PASS |
| disp_field_linf.y | 1.76884e-05 | 1.76884e-05 | 6.093e-16 | 3.445e-11 | 1.000e-06 | PASS |
| disp_field_linf.z | 4.08954e-06 | 4.08954e-06 | 3.616e-17 | 8.842e-12 | 1.000e-06 | PASS |
| stress_field_l2.xx | 962.999 | 962.999 | 1.263e-07 | 1.312e-10 | 1.000e-06 | PASS |
| stress_field_l2.xy | 567.881 | 567.881 | 3.160e-08 | 5.565e-11 | 1.000e-06 | PASS |
| stress_field_l2.xz | 228.282 | 228.282 | 1.780e-08 | 7.799e-11 | 1.000e-06 | PASS |
| stress_field_l2.yx | 567.881 | 567.881 | 3.160e-08 | 5.565e-11 | 1.000e-06 | PASS |
| stress_field_l2.yy | 962.999 | 962.999 | 1.238e-07 | 1.286e-10 | 1.000e-06 | PASS |
| stress_field_l2.yz | 228.282 | 228.282 | 1.520e-08 | 6.661e-11 | 1.000e-06 | PASS |
| stress_field_l2.zx | 228.296 | 228.296 | 1.780e-08 | 7.799e-11 | 1.000e-06 | PASS |
| stress_field_l2.zy | 228.296 | 228.296 | 1.521e-08 | 6.661e-11 | 1.000e-06 | PASS |
| stress_field_l2.zz | 449.323 | 449.323 | 9.044e-08 | 2.013e-10 | 1.000e-06 | PASS |
| stress_field_linf.xx | 65801.1 | 65801.1 | 7.758e-08 | 1.179e-12 | 1.000e-06 | PASS |
| stress_field_linf.xy | 40288 | 40288 | 1.480e-06 | 3.674e-11 | 1.000e-06 | PASS |
| stress_field_linf.xz | 92476.4 | 92476.4 | 3.620e-07 | 3.914e-12 | 1.000e-06 | PASS |
| stress_field_linf.yx | 40288 | 40288 | 1.480e-06 | 3.673e-11 | 1.000e-06 | PASS |
| stress_field_linf.yy | 65801.1 | 65801.1 | 2.682e-07 | 4.076e-12 | 1.000e-06 | PASS |
| stress_field_linf.yz | 92476.4 | 92476.4 | 1.420e-07 | 1.536e-12 | 1.000e-06 | PASS |
| stress_field_linf.zx | 92406.3 | 92406.3 | 3.610e-07 | 3.907e-12 | 1.000e-06 | PASS |
| stress_field_linf.zy | 92406.3 | 92406.3 | 1.411e-07 | 1.527e-12 | 1.000e-06 | PASS |
| stress_field_linf.zz | 75917.8 | 75917.8 | 2.593e-07 | 3.416e-12 | 1.000e-06 | PASS |
| strain_field_l2.xx | 0.0526181 | 0.0526181 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| strain_field_l2.xy | 6.73497e-06 | 6.73497e-06 | 2.863e-17 | 4.252e-12 | 1.000e-06 | PASS |
| strain_field_l2.xz | 7.24055e-06 | 7.24055e-06 | 4.044e-17 | 5.585e-12 | 1.000e-06 | PASS |
| strain_field_l2.yx | 6.73497e-06 | 6.73497e-06 | 9.471e-17 | 1.406e-11 | 1.000e-06 | PASS |
| strain_field_l2.yy | 0.0526181 | 0.0526181 | 1.388e-17 | 2.637e-16 | 1.000e-06 | PASS |
| strain_field_l2.yz | 7.24055e-06 | 7.24055e-06 | 6.847e-17 | 9.456e-12 | 1.000e-06 | PASS |
| strain_field_l2.zx | 1.20377e-06 | 1.20377e-06 | 3.079e-17 | 2.558e-11 | 1.000e-06 | PASS |
| strain_field_l2.zy | 1.20377e-06 | 1.20377e-06 | 6.210e-18 | 5.158e-12 | 1.000e-06 | PASS |
| strain_field_l2.zz | 0.0526182 | 0.0526182 | 1.388e-17 | 2.637e-16 | 1.000e-06 | PASS |
| elastic_residual_at_convergence | 9.72181e-09 | 9.72181e-09 | 7.400e-17 | 7.612e-09 | 1.000e-08 | PASS |

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
| max_von_mises_stress | 8.784e-12 | 1.482e-06 | -- | 0.05 | True | PASS |
| mean_von_mises_stress | 1.254e-10 | 2.956e-06 | -- | 0.05 | True | PASS |
| max_principal_stress | 1.136e-11 | 1.622e-06 | -- | 0.05 | True | PASS |
| elastic_energy | 2.128e-10 | 9.677e-09 | -- | 0.03 | True | PASS |

### SOLVER-HEALTH (non-gating)

| observable | A | B | delta |
| --- | --- | --- | --- |
| newton_iters_per_solve | 2 | 2 | 0 |
| newton_final_relnorm | 2.72291e-06 | 2.72291e-06 | 0.0 |
| mlmg_vcycles_per_solve | 30 | 30 | 0 |
| bottom_solver_iters | {'mean': 4.8011695906432745, 'max': 6, 'n': 342} | {'mean': 4.8011695906432745, 'max': 6, 'n': 342} | 0.0 |
| residual_monotonic | False | False | 0 |
