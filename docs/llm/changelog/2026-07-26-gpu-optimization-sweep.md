# GPU optimization sweep — smoother fusion, read-only coefficients, Fapply launch width — 2026-07-26

Branch `chamber-gpu`. Task folder:
`docs/agent_plans/20260726-gpu-optimization-sweep/` (evidence in
`results/RESULT.md`, traces in `artifacts/nsys/`).

## TL;DR

| Change | Effect (2D conservative / 3D section, medians of 5) | Status |
| --- | --- | --- |
| `Operator::Fsmooth` elementwise fusion — drop `Dx`/`Rx` MultiFabs and 4 Copy/Multiply/Subtract kernels per Jacobi half-sweep | -18.0% / -3.0% wall | **retained** |
| `Fapply` reads `DDW`/`psi` as `const_array`; `ddw` bound by reference instead of copying 45 doubles | -7.1% / -2.0% wall | **retained** |
| `Fapply` launch width 256 -> 128 (`ALAMO_ELASTIC_FAPPLY_MT`) | -9.8% / -3.4% wall | **retained** |
| `Set::MulCol` in the conservative face-flux branch | +0.6%, noise | reverted |
| `Fapply` launch width 64 | +0.7% vs 128 | not adopted |

Cumulative, paired: **-24.2%** wall 2D, **-9.8%** wall 3D section;
GPU kernel time -28.4%; kernel launches -56.2% (74,756 -> 32,772);
`Fapply` kernel time -22.8% at an unchanged 11,808 calls.

## Why the launch width mattered

Grouping 2D `Fapply` launches by grid size showed 3,936 launches with a grid of
**2 blocks** on a 16-SM device, each still ~30 us, and 55% of `Fapply` time in
launches of <= 17 blocks. Coarse MG levels, not the fine level, dominate the
launch population. Occupancy per SM is unchanged by the block width (254
regs/thread caps an SM near 256 resident threads either way); what changes is
how many SMs are engaged at all.

## Correctness

Strict (`--fmad=false`) CUDA builds of the retained set are **bit-identical** to
the unmodified source on `input` and `input_rod_and_tube_2d`. Device lint, CPU
golden compare (4/4) and `local_a100_gate` TIERS=1 all PASS.

Fast-math builds differ at ~1e-9 relative because the fused `Ax - x*diag` can
contract to an FMA; MLMG iteration counts are unchanged.

## Pre-existing red gate confirmed, not fixed

`rod_and_tube_step2/gpu_strict` and `/gpu_fast` fail with **identical deltas**
using the untouched 2026-07-24 binary: the GPU references date from 2026-07-05
and were never regenerated after the 2026-07-18 CPU recalibration (already noted
in `2026-07-18-rodtube-gate-recalibration.md`). Fresh GPU output now agrees with
the CPU reference. Re-recording the GPU references is left as a separate call.
`benchmark/status.sh` does not exercise this leg, which is why it stayed unseen.

## Notes

- The 3D-psi case `input_3d_centre_bore_128_a2` is not locally measurable:
  ~26 GB managed high-water on a 7.8 GB card, >18 min/rep. The `use_psi=1`
  (non-conservative) `Fapply` branch is unmeasured here and belongs on NOVA.
- Nsight Compute is unavailable unprivileged on this driver
  (`ERR_NVGPUCTRPERM`); `scripts/run_root_ncu.sh` is ready for a `sudo` run.
- Next candidates, ranked, are in the task `results/RESULT.md`: face-flux kernel
  fission, shared-memory tiling of `U`, coarse-level launch consolidation,
  `m_ddw_mf` AoS->SoA, and dropping the unread face components of `m_ddw_mf` on
  the psi path (75% off the largest array).
