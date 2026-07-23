# FApply runtime optimization result

## Outcome

Retain configuration-only `elastic.solver.pre_smooth=2` and
`elastic.solver.post_smooth=2`. No FApply source optimization is retained;
protected Elastic source matches HEAD `2e6a8f8f5430a58dd6a85b07b08d472e2bd1d5cd`.

The accepted evidence is limited to the frozen two-step cases. It passes both
the local A1000 gates and a matched NOVA/A100 confirmation; it is not a
longer-evolution stability claim.

## Retained result

Warmups are excluded; values are medians of five measurements per arm/mode.

| regime | setting | external wall (s) | MLMG solve (s) | FApply calls | FApply wall (s) |
|---|---|---:|---:|---:|---:|
| 2D conservative | 4/4 drift | 1.72 | 0.9444 | 11,808 | 0.4630 |
| 2D conservative | 2/2 | 1.50 | 0.7181 | 8,904 | 0.3496 |
| 3D psi | 4/4 drift | 288.58 | 230.5 | 15,477 | 211.9 |
| 3D psi | 2/2 | 227.84 | 176.5 | 12,660 | 161.6 |

Matched 2/2 gains are 12.791% external / 23.962% solve in 2D and
21.048% external / 23.427% solve in 3D. Correctness passed for conservative
single-/multi-box and 3D psi physics-budget cases. Newton iterations remain two;
the 3D final residual is `9.72181e-09`, inside but close to the `1e-08` gate.
This is not a longer-evolution stability claim.

Primary evidence:

- `step2/smoothing_summary.md` and `step2/smoothing_summary.json`
- `step2/raw/correctness_2x2.log`
- `step2/raw/timing_pair_2d.log`
- `step2/raw/timing_pair_3d_rep5_recovery.log`

All paths above are under
`artifacts/a1000-sm86-2e6a8f8f-20260721/` in this task directory.

## NOVA/A100 confirmation

An isolated sm_80 profile-fast build at the same HEAD used binary SHA-256
`036c7ec520c6ae301b31d1e75a107b5dff2e343f02729c880f7be3ceeda98c84` on an
NVIDIA A100 80GB PCIe. Warmups are excluded; values are median ± MAD of five
alternating measurements per arm.

| metric | 4/4 | 2/2 | reduction |
|---|---:|---:|---:|
| external wall | 15.58 ± 0.02 s | 13.00 ± 0.01 s | 16.560% |
| MLMG solve | 11.42 ± 0.01 s | 8.846 ± 0.002 s | 22.539% |
| FApply wall | 8.388 ± 0.008 s | 6.434 ± 0.001 s | 23.295% |
| FApply calls | 15,477 | 12,660 | 18.201% |

The two A100 physics bundles passed the overall and gating comparisons. The
same residual caveat remains: 2/2 ends at `9.72181e-09`, inside but close to the
`1e-08` gate. A one-launch Nsight Compute capture reports 254 registers/thread,
one register-limited block/SM, and 12.50% theoretical / 12.02% achieved
occupancy for the hot FApply kernel.

Primary A100 evidence is under
`artifacts/nova-a100-sm80-2e6a8f8f-20260722/`; its `SUMMARY.md` records the job
chain, profiler recovery, hashes, and compact-copy scope.

## Rejected experiments

| experiment | correctness | decisive performance result | disposition |
|---|---|---|---|
| Step 3 conservative specialization | pass | 2D external +2.013% gain, below 3%; 3D external 1.520% regression | reverted |
| Step 4 sequential Cgrad contraction | quick runtime completed | 3D external +0.979%, MLMG solve +0.519%, FApply +0.581% | reverted |

Step 4 reduced the static hot-kernel resource footprint from 254 to 250
registers and from 440 to 288 stack bytes, but that did not translate into a
material runtime gain. Its abbreviated rejected-candidate timing is not used as
accepted-result evidence.

## Conditional branches

- Psi caching was skipped: the 3D case uses a 26,059 MB managed-arena high-water
  mark on a 7,829 MB physical GPU and fails the mandatory 1 GiB-headroom gate.
- Multi-box launch fusion was skipped: the retained 2/2 trace contains exactly
  8,904 FApply kernels for 8,904 FApply calls, already one launch per invocation
  in the target case. There is no per-box launch multiplicity to eliminate.
- Fresh antagonistic review found no remaining blocker to local-finalist
  acceptance. Findings and adjudications are in `results/REVIEW.md`.

## Remaining risk and ship status

The retained 2/2 configuration is cleared by both the A1000 performance contract
and the NOVA/A100 pre-ship confirmation. No source optimization is retained.
The stability claim remains restricted to `max_step=2`; a longer production
evolution would be a separate qualification task.
