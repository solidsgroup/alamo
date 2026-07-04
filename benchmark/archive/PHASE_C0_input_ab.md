# Phase C0 Input A/B

## 2026-06-30 - First Phase C NOVA sweep (`gpu_3d_a2_logs.tgz`)

Source bundle: `~/Downloads/gpu_3d_a2_logs.tgz`

Scope: analyze the first Phase C NOVA A100 sweep for the cheap input levers in
`input_3d_centre_bore_256_a2_tuned` versus the untuned `..._256_a2` deck. This
bundle contains log output only; no plotfiles, stress compares, or solver-side
CSV exports are included.

### Jobs in bundle

All jobs used the same binary and device shape:

- Binary: `bin/alamo_gpu-3d-profile-cuda80-g++`
- Device: A100 / `sm_80`
- MPI shape: 1 rank
- Mode: `bench`

Exploratory scavenger runs on 2026-06-28:

| Job | Partition | Input | Outcome |
|---|---|---|---|
| `11307829` | `scavenger` | `input_3d_centre_bore_128_a2` | preempted after 50 steps |
| `11307831` | `scavenger` | `input_3d_centre_bore_256_a2` | preempted after 15,150 steps |

First proper NOVA Phase C sweep on 2026-06-29:

| Job | Partition | Input | Outcome | Last completed step | Sim time reached |
|---|---|---|---|---:|---:|
| `11333287` | `nova` | `input_3d_centre_bore_256_a2` | time limit | 5,800 | 0.5800 |
| `11333288` | `nova` | `input_3d_centre_bore_256_a2` | time limit | 5,600 | 0.5600 |
| `11333289` | `nova` | `input_3d_centre_bore_256_a2` | time limit | 5,800 | 0.5800 |
| `11333290` | `nova` | `input_3d_centre_bore_256_a2` | time limit | 5,600 | 0.5600 |
| `11333291` | `nova` | `input_3d_centre_bore_256_a2_tuned` | time limit | 14,500 | 1.4500 |

### Main result

Under the same one-hour `nova` partition limit, the tuned 256^3 deck advanced
substantially farther before timeout than the untuned baseline:

- Baseline `256_a2`: 5,600-5,800 completed steps before timeout
- Tuned `256_a2_tuned`: 14,500 completed steps before timeout
- Ratio versus baseline: `14500 / 5800 = 2.50x` at the best baseline row,
  `14500 / 5600 = 2.59x` at the worst baseline row
- Ratio versus the four-run baseline mean (5,700 steps): `2.54x`

This is the useful signal in the bundle. The 2026-06-28 scavenger runs were
interrupted by preemption and should be treated as exploratory only.

### Interpretation

The C0 input levers are directionally successful on NOVA A100:

- The tuned deck makes about 2.5x more wall-clock progress than the untuned
  baseline under the same partition and binary.
- No job in the bundle shows a CUDA fault, AMReX abort, or numerical blow-up.
- The gain is large enough to justify carrying C0 forward before source-level
  kernel work.

What this bundle does **not** prove yet:

- No run completed to the intended endpoint, so this is not yet a full A/B
  wall-time result.
- No stress-field parity artifact is present, so correctness of the relaxed
  tolerances is not closed from this bundle alone.
- No per-solve wall, V-cycle count, or bottom-solver iteration table is present
  in stdout, so attribution across `tol`, `bottom_*`, and `interval` is still
  missing.
- Only one tuned run is present.

### Disposition

Status: **promising but incomplete**.

The first Phase C NOVA sweep supports the C0 hypothesis strongly enough to keep
the tuned deck in play, but it does **not** satisfy the roadmap done-when yet.
The missing pieces are:

1. Complete the explicit A/B matrix `{baseline, tol, bottom, interval, all}`.
2. Record stress-field parity with `benchmark/compare_thermo.py`.
3. Capture per-solve wall and V-cycle counts so the win can be attributed, not
   just observed as more completed timesteps before the wall limit.

### Recommended next NOVA run

- Use the dedicated `benchmark/phase_c_elastic_ab.sh` harness rather than raw
  one-off slurm logs.
- Give the tuned case enough wall time to finish the target stop condition, or
  reduce the stop condition so all A/B rows complete within one allocation.
- Preserve the current binary/device shape for comparability:
  `bin/alamo_gpu-3d-profile-cuda80-g++`, single-rank A100, `nova` partition.
