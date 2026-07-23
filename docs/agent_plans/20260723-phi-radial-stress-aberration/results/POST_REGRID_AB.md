# Controlled post-regrid mechanics A/B

Verdict: **solution 1 is rejected as a fix for the visible radial-stress
aberration**.

The controlled pair uses the same frozen restart at
`1.0003999471664429 s`. Arm A is the restored `05002node` state; arm B is
`05003node` after one mechanics solve on the unchanged hierarchy. The solve
took two Newton iterations and accepted an update of `4.96028e-6`, with the
logged nonlinear residual `2.32284e6 Pa/m`.

Pair validity:

- all three hierarchy coverage hashes match;
- all 33 `phi`, `eta`, temperature, material, and RHS hashes match;
- all 18 displacement and stress hashes change, proving that the mechanics
  solve occurred;
- no physics-time advance or regrid occurs between arms.

The radial metrics do not meet the required 80% reduction:

| Axis / representation | Arm A Nyquist | Arm B Nyquist | Reduction |
|---|---:|---:|---:|
| bottom saved polar | `121.5718 Pa` | `121.5586 Pa` | `0.0108%` |
| bottom native face | `49.2044 Pa` | `49.1933 Pa` | `0.0227%` |
| right saved polar | `282.5042 Pa` | `282.4803 Pa` | `0.00845%` |
| right native face | `90.2255 Pa` | `90.2398 Pa` | `-0.0159%` |

All four monotone-envelope overshoot values are exactly zero in both arms.
The detrended maxima are also unchanged: approximately `15.97 kPa` on the
bottom native faces and `3.70 kPa` on the right native faces.

The synchronized arm already agrees with the exact face-average reconstruction
to `1.88e-5 Pa` in the finest-level `phi` band. The offline band residual in
this oracle is not the same norm/region as the Newton acceptance residual and
must not be used as a solver-convergence claim.

Conclusion: solving after regrid remains independently justified to prevent
stale derived mechanics output, but it does not remove this visible band. No
production lifecycle change is retained by this task because the requested
visual-aberration threshold was missed by more than three orders of magnitude.

Machine-readable comparison:
[`controlled-comparison/comparison.json`](controlled-comparison/comparison.json).

