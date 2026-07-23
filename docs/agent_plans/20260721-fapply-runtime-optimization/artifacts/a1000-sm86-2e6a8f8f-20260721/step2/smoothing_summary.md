# Step 2 paired smoothing summary

All values below exclude the warmup and use five measured repetitions.

| regime | setting | external wall median ± MAD (s) | MLMG solve median ± MAD (s) | FApply calls | FApply wall median ± MAD (s) |
|---|---:|---:|---:|---:|---:|
| 2d_conservative | drift_4x4 | 1.72 ± 0.01 | 0.9444 ± 0.004 | 11808 | 0.463 ± 0.002 |
| 2d_conservative | 2x2 | 1.5 ± 0 | 0.7181 ± 0.0073 | 8904 | 0.3496 ± 0.0026 |

- 2d_conservative matched 2x2 gains: external 12.791%, MLMG solve 23.962%, FApply wall 24.492%, calls 24.593%.

| 3d_psi | drift_4x4 | 288.58 ± 0.15 | 230.5 ± 4.1 | 15477 | 211.9 ± 3.5 |
| 3d_psi | 2x2 | 227.84 ± 2.02 | 176.5 ± 1.4 | 12660 | 161.6 ± 1.1 |

- 3d_psi matched 2x2 gains: external 21.048%, MLMG solve 23.427%, FApply wall 23.738%, calls 18.201%.

