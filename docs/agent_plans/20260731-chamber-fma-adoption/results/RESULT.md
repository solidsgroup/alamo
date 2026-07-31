# RESULT — chamber-fma-adoption

Status: **PAUSED AT STEP 3**

## Adopted source state

The current `chamber-gpu` target before adoption was
`508a8785d972191235a80394614fc2e33a7ef57c`. The adoption plan was committed as
`c1f49caba6f415dd301d189f973cff993482d3c5`, followed by the reviewed source
pair:

- `8110c5e5c` — `perf: skip interior physical boundary fills`
- `f73154f7e` — `fix: preserve index type in BC interior skip`

The resulting source diff is exactly 13 insertions in `src/BC/BC.H` and
`src/Integrator/BaseField.H`. Its source-patch SHA-256 is
`6369167a6b823b60760d6a258cdbc71bc6d2dd72a29f96418e2614183f634585`,
identical to the frozen campaign patch.

## Completed target checks

- `GOLDEN_MODE=cpu BUILD_JOBS=8 benchmark/ci_golden_compare.sh`: PASS.
  All four CPU golden cases and the NaN-flag smoke passed.
- Strict GPU build with CUDA 12.6.85: build PASS. The installed CUDA 12.0
  `nvcc` wrapper was unusable because it delegates to missing
  `/usr/lib/nvidia-cuda-toolkit/bin/nvcc`; the self-contained
  `/home/jackplum/Projects/alamo/.local/cuda-12.6.3-redist` toolchain was used.
- Three strict GPU golden cases passed:
  `canonical_step1`, `canonical_step2`, and `eta_expression_step1`.

## Blocking target discrepancy

`rod_and_tube_step2/gpu_strict` fails four traction observables against the
stored reference:

| Observable | Absolute delta | Relative delta |
| --- | ---: | ---: |
| `trac_xhi_x` | `1.384000e+02` | `2.265232e-03` |
| `trac_xhi_y` | `4.585000e-01` | `2.349704e-02` |
| `trac_yhi_x` | `6.512000e-01` | `2.642621e-01` |
| `trac_yhi_y` | `1.256000e+02` | `2.100370e-03` |

An independent detached-worktree diagnostic ran the identical strict GPU gate
at pre-adoption target `508a8785d`. It reproduced all four failures with the
same deltas. The temporary worktree was removed afterward. This proves the
failure predates OPT-10, but it leaves the required target oracle red.

Per the task operating rules, sanitizer and target-bound 800-step production
validation have not advanced. Updating the scientific reference or waiving the
gate is outside this adoption task and requires explicit adjudication.
