# Step 3 isolated timing summary

Warmups are excluded; each entry is the median ± MAD of five runs.

| regime | arm | external wall (s) | MLMG solve (s) | FApply calls | FApply wall (s) | FApply per call (µs) |
|---|---|---:|---:|---:|---:|---:|
| 2d_conservative_2x2 | baseline | 1.49 ± 0 | 0.7008 ± 0.0113 | 8904 | 0.3433 ± 0.0017 | 38.5557 ± 0.191 |
| 2d_conservative_2x2 | candidate | 1.46 ± 0.01 | 0.7577 ± 0.0164 | 8904 | 0.3414 ± 0.0013 | 38.3423 ± 0.146 |

- 2d_conservative_2x2 candidate gains: external 2.013%, MLMG solve -8.119%, FApply wall 0.553%, per-call 0.553%.

| 3d_psi_2x2 | baseline | 224.93 ± 0.64 | not run (fail-fast) | — | — | — |
| 3d_psi_2x2 | candidate | 228.35 ± 5.2 | not run (fail-fast) | — | — | — |

- 3d_psi_2x2 candidate gains: external -1.520%.

