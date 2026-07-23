# Physics validation compare report

- reference: `/work/brunnels/jackplum/alamo-fapply-20260722/nova-results-11751561/validation_4x4`
- candidate: `/work/brunnels/jackplum/alamo-fapply-20260722/nova-results-11751561/validation_2x2`
- overall verdict: **PASS** (gate verdict: **PASS**)

## centre_bore_3d_128_a2_converged -- PASS (gate: PASS)

### CORRECTNESS

| observable | A | B | abs Δ | rel Δ | tol | status |
| --- | --- | --- | --- | --- | --- | --- |
| eta_field_l2 | 0.0496029 | 0.0496029 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| eta_field_linf | 1 | 1 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| phi_field_l2 | 0.0414434 | 0.0414434 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| temp_field_l2 | 15.7876 | 15.7876 | 0.000e+00 | 0.000e+00 | 1.000e-06 | PASS |
| disp_field_l2.x | 1.43148e-07 | 1.43148e-07 | 3.645e-17 | 2.546e-10 | 1.000e-06 | PASS |
| disp_field_l2.y | 1.43148e-07 | 1.43148e-07 | 2.761e-17 | 1.929e-10 | 1.000e-06 | PASS |
| disp_field_l2.z | 2.86443e-08 | 2.86443e-08 | 2.629e-16 | 9.179e-09 | 1.000e-06 | PASS |
| disp_field_linf.x | 1.76884e-05 | 1.76884e-05 | 1.426e-15 | 8.064e-11 | 1.000e-06 | PASS |
| disp_field_linf.y | 1.76884e-05 | 1.76884e-05 | 3.305e-16 | 1.868e-11 | 1.000e-06 | PASS |
| disp_field_linf.z | 4.08954e-06 | 4.08954e-06 | 7.290e-14 | 1.783e-08 | 1.000e-06 | PASS |
| stress_field_l2.xx | 962.999 | 962.999 | 3.974e-08 | 4.127e-11 | 1.000e-06 | PASS |
| stress_field_l2.xy | 567.881 | 567.881 | 5.983e-09 | 1.054e-11 | 1.000e-06 | PASS |
| stress_field_l2.xz | 228.282 | 228.282 | 4.021e-07 | 1.762e-09 | 1.000e-06 | PASS |
| stress_field_l2.yx | 567.881 | 567.881 | 5.990e-09 | 1.055e-11 | 1.000e-06 | PASS |
| stress_field_l2.yy | 962.999 | 962.999 | 1.389e-07 | 1.442e-10 | 1.000e-06 | PASS |
| stress_field_l2.yz | 228.282 | 228.282 | 3.687e-07 | 1.615e-09 | 1.000e-06 | PASS |
| stress_field_l2.zx | 228.296 | 228.296 | 4.014e-07 | 1.758e-09 | 1.000e-06 | PASS |
| stress_field_l2.zy | 228.296 | 228.296 | 3.680e-07 | 1.612e-09 | 1.000e-06 | PASS |
| stress_field_l2.zz | 449.323 | 449.323 | 4.200e-08 | 9.347e-11 | 1.000e-06 | PASS |
| stress_field_linf.xx | 65801.1 | 65801.1 | 1.570e-06 | 2.386e-11 | 1.000e-06 | PASS |
| stress_field_linf.xy | 40288 | 40288 | 1.716e-07 | 4.258e-12 | 1.000e-06 | PASS |
| stress_field_linf.xz | 92476.4 | 92476.4 | 2.931e-03 | 3.170e-08 | 1.000e-06 | PASS |
| stress_field_linf.yx | 40288 | 40288 | 1.729e-07 | 4.292e-12 | 1.000e-06 | PASS |
| stress_field_linf.yy | 65801.1 | 65801.1 | 1.416e-06 | 2.152e-11 | 1.000e-06 | PASS |
| stress_field_linf.yz | 92476.4 | 92476.4 | 2.934e-03 | 3.173e-08 | 1.000e-06 | PASS |
| stress_field_linf.zx | 92406.3 | 92406.3 | 2.926e-03 | 3.166e-08 | 1.000e-06 | PASS |
| stress_field_linf.zy | 92406.3 | 92406.3 | 2.929e-03 | 3.169e-08 | 1.000e-06 | PASS |
| stress_field_linf.zz | 75917.8 | 75917.8 | 1.553e-03 | 2.046e-08 | 1.000e-06 | PASS |
| strain_field_l2.xx | 0.0526181 | 0.0526181 | 4.857e-17 | 9.231e-16 | 1.000e-06 | PASS |
| strain_field_l2.xy | 6.73497e-06 | 6.73497e-06 | 2.429e-15 | 3.606e-10 | 1.000e-06 | PASS |
| strain_field_l2.xz | 7.24055e-06 | 7.24055e-06 | 1.889e-14 | 2.609e-09 | 1.000e-06 | PASS |
| strain_field_l2.yx | 6.73497e-06 | 6.73497e-06 | 2.186e-15 | 3.245e-10 | 1.000e-06 | PASS |
| strain_field_l2.yy | 0.0526181 | 0.0526181 | 1.388e-16 | 2.637e-15 | 1.000e-06 | PASS |
| strain_field_l2.yz | 7.24055e-06 | 7.24055e-06 | 1.869e-14 | 2.582e-09 | 1.000e-06 | PASS |
| strain_field_l2.zx | 1.20377e-06 | 1.20377e-06 | 1.090e-14 | 9.058e-09 | 1.000e-06 | PASS |
| strain_field_l2.zy | 1.20377e-06 | 1.20377e-06 | 1.100e-14 | 9.138e-09 | 1.000e-06 | PASS |
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
| max_von_mises_stress | 2.350e-08 | 3.963e-03 | -- | 0.05 | True | PASS |
| mean_von_mises_stress | 3.152e-11 | 7.433e-07 | -- | 0.05 | True | PASS |
| max_principal_stress | 1.333e-08 | 1.903e-03 | -- | 0.05 | True | PASS |
| elastic_energy | 1.132e-10 | 5.150e-09 | -- | 0.03 | True | PASS |

### SOLVER-HEALTH (non-gating)

| observable | A | B | delta |
| --- | --- | --- | --- |
| newton_iters_per_solve | 2 | 2 | 0 |
| newton_final_relnorm | 2.72291e-06 | 2.72291e-06 | 0.0 |
| mlmg_vcycles_per_solve | 24 | 30 | 6 |
| bottom_solver_iters | {'mean': 5.09469696969697, 'max': 6, 'n': 264} | {'mean': 4.8011695906432745, 'max': 6, 'n': 342} | -0.2935273790536952 |
| residual_monotonic | False | False | 0 |
