# Existing-pair verdict

Verdict: **invalid as a controlled A/B; historical evidence only**.

The original `05000node` and the later final-hierarchy `05002node` do not
represent the same frozen state:

- output times are `0.9999999999999225 s` and `1.0003999999999225 s`;
- level-2 AMR coverage differs;
- `eta`, temperature, material, RHS, and level-2 `phi` hashes differ.

The pair therefore cannot select a production fix. It does reproduce the
earlier lifecycle finding: in the current `phi=0.1-0.9` band, the saved polar
stress versus exact face-average reconstruction mismatch falls from
`8.232 kPa` to `1.93e-5 Pa`. The earlier broader-band diagnostic measured the
same defect as `48.60 kPa` and showed the interior residual falling from
`815.7 MPa/m` to `2.32284 MPa/m` after solving on the final hierarchy.

The native radial feature itself does not materially change in this historical
comparison:

| Axis | Existing native Nyquist | Final-hierarchy native Nyquist | Change |
|---|---:|---:|---:|
| bottom | `49.2145 Pa` | `49.2044 Pa` | `-0.0205%` |
| right | `90.2413 Pa` | `90.2255 Pa` | `-0.0175%` |

This pair supports post-regrid solving as a hierarchy-consistency correction,
but it does not show that the correction removes the visible `phi`-following
radial-stress band.

Machine-readable comparison:
[`existing-comparison/comparison.json`](existing-comparison/comparison.json).

