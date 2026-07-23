# Physics validation compare report

- reference: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/validation/fast_3d`
- candidate: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/candidate-step3/validation/fast_3d`
- overall verdict: **PASS** (gate verdict: **PASS**)

## centre_bore_3d_128_a2_converged -- PASS (gate: PASS)

### CORRECTNESS

| observable | A | B | abs Δ | rel Δ | tol | status |
| --- | --- | --- | --- | --- | --- | --- |
| eta_field_l2 | 0.0496029 | 0.0496029 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| eta_field_linf | 1 | 1 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| phi_field_l2 | 0.0414434 | 0.0414434 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| temp_field_l2 | 15.7876 | 15.7876 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| disp_field_l2.x | 1.43148e-07 | 1.43148e-07 | 7.406e-18 | 5.174e-11 | 1.000e-06 | PASS |
| disp_field_l2.y | 1.43148e-07 | 1.43148e-07 | 2.363e-18 | 1.651e-11 | 1.000e-06 | PASS |
| disp_field_l2.z | 2.86443e-08 | 2.86443e-08 | 1.332e-18 | 4.650e-11 | 1.000e-06 | PASS |
| disp_field_linf.x | 1.76884e-05 | 1.76884e-05 | 8.666e-16 | 4.899e-11 | 1.000e-06 | PASS |
| disp_field_linf.y | 1.76884e-05 | 1.76884e-05 | 6.687e-17 | 3.780e-12 | 1.000e-06 | PASS |
| disp_field_linf.z | 4.08954e-06 | 4.08954e-06 | 1.077e-16 | 2.634e-11 | 1.000e-06 | PASS |
| stress_field_l2.xx | 962.999 | 962.999 | 1.989e-07 | 2.066e-10 | 1.000e-06 | PASS |
| stress_field_l2.xy | 567.881 | 567.881 | 2.707e-08 | 4.767e-11 | 1.000e-06 | PASS |
| stress_field_l2.xz | 228.282 | 228.282 | 3.139e-08 | 1.375e-10 | 1.000e-06 | PASS |
| stress_field_l2.yx | 567.881 | 567.881 | 2.708e-08 | 4.769e-11 | 1.000e-06 | PASS |
| stress_field_l2.yy | 962.999 | 962.999 | 7.081e-08 | 7.353e-11 | 1.000e-06 | PASS |
| stress_field_l2.yz | 228.282 | 228.282 | 5.059e-09 | 2.216e-11 | 1.000e-06 | PASS |
| stress_field_l2.zx | 228.296 | 228.296 | 3.140e-08 | 1.375e-10 | 1.000e-06 | PASS |
| stress_field_l2.zy | 228.296 | 228.296 | 5.054e-09 | 2.214e-11 | 1.000e-06 | PASS |
| stress_field_l2.zz | 449.323 | 449.323 | 9.426e-08 | 2.098e-10 | 1.000e-06 | PASS |
| stress_field_linf.xx | 65801.1 | 65801.1 | 3.224e-07 | 4.900e-12 | 1.000e-06 | PASS |
| stress_field_linf.xy | 40288 | 40288 | 2.371e-06 | 5.884e-11 | 1.000e-06 | PASS |
| stress_field_linf.xz | 92476.4 | 92476.4 | 2.672e-06 | 2.890e-11 | 1.000e-06 | PASS |
| stress_field_linf.yx | 40288 | 40288 | 2.372e-06 | 5.888e-11 | 1.000e-06 | PASS |
| stress_field_linf.yy | 65801.1 | 65801.1 | 5.894e-08 | 8.957e-13 | 1.000e-06 | PASS |
| stress_field_linf.yz | 92476.4 | 92476.4 | 9.769e-07 | 1.056e-11 | 1.000e-06 | PASS |
| stress_field_linf.zx | 92406.3 | 92406.3 | 2.664e-06 | 2.883e-11 | 1.000e-06 | PASS |
| stress_field_linf.zy | 92406.3 | 92406.3 | 9.734e-07 | 1.053e-11 | 1.000e-06 | PASS |
| stress_field_linf.zz | 75917.8 | 75917.8 | 5.218e-06 | 6.873e-11 | 1.000e-06 | PASS |
| strain_field_l2.xx | 0.0526181 | 0.0526181 | 6.939e-18 | 1.319e-16 | 1.000e-06 | PASS |
| strain_field_l2.xy | 6.73497e-06 | 6.73497e-06 | 1.898e-16 | 2.819e-11 | 1.000e-06 | PASS |
| strain_field_l2.xz | 7.24055e-06 | 7.24055e-06 | 1.641e-16 | 2.266e-11 | 1.000e-06 | PASS |
| strain_field_l2.yx | 6.73497e-06 | 6.73497e-06 | 1.190e-16 | 1.767e-11 | 1.000e-06 | PASS |
| strain_field_l2.yy | 0.0526181 | 0.0526181 | 4.857e-17 | 9.231e-16 | 1.000e-06 | PASS |
| strain_field_l2.yz | 7.24055e-06 | 7.24055e-06 | 5.909e-17 | 8.160e-12 | 1.000e-06 | PASS |
| strain_field_l2.zx | 1.20377e-06 | 1.20377e-06 | 3.949e-17 | 3.280e-11 | 1.000e-06 | PASS |
| strain_field_l2.zy | 1.20377e-06 | 1.20377e-06 | 7.631e-17 | 6.339e-11 | 1.000e-06 | PASS |
| strain_field_l2.zz | 0.0526182 | 0.0526182 | 6.939e-18 | 1.319e-16 | 1.000e-06 | PASS |
| elastic_residual_at_convergence | 4.45907e-09 | 4.45907e-09 | 1.225e-15 | 2.747e-07 | 1.000e-08 | PASS |

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
| max_von_mises_stress | 1.369e-11 | 2.308e-06 | -- | 0.05 | True | PASS |
| mean_von_mises_stress | 1.263e-10 | 2.979e-06 | -- | 0.05 | True | PASS |
| max_principal_stress | 2.183e-11 | 3.118e-06 | -- | 0.05 | True | PASS |
| elastic_energy | 2.173e-10 | 9.883e-09 | -- | 0.03 | True | PASS |

### SOLVER-HEALTH (non-gating)

| observable | A | B | delta |
| --- | --- | --- | --- |
| newton_iters_per_solve | 2 | 2 | 0 |
| newton_final_relnorm | 2.72291e-06 | 2.72291e-06 | 0.0 |
| mlmg_vcycles_per_solve | 24 | 24 | 0 |
| bottom_solver_iters | {'mean': 5.09469696969697, 'max': 6, 'n': 264} | {'mean': 5.09469696969697, 'max': 6, 'n': 264} | 0.0 |
| residual_monotonic | False | False | 0 |
