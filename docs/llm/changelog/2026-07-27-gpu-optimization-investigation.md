# GPU optimization investigation — smoothing config win, coarsening rejected — 2026-07-27

Branch `chamber-gpu`, HEAD `19132c02c`. Task folder
`docs/agent_plans/20260727-gpu-optimization-investigation/`.
**No source change retained**; `src/` clean at HEAD. Gates green before and
after.

## TL;DR

Broad sweep of the GPU optimization space on `Operator::Elastic`. One
configuration win, four measured rejections, several candidates closed on
existing evidence.

| finding | result | disposition |
| --- | --- | --- |
| pre/post smoothing 4/4 -> 2/2 | -10.5% (2D), **-25.3%** (3D section), -6.9% (rod_and_tube), -5.9% (production `input`); 0 failures | **retain — recommended, not applied** |
| `max_coarsening_level=1` | -9.4% in 2D but **aborts** `input_rod_and_tube_2d` and is +11.6% in 3D | rejected |
| deeper coarsening (mcl 3/4/6) | +12.0% / +27.5% / +46.7% for <= 1 fewer iteration | rejected |
| `final_smooth` 8 -> 4 / 2 | bottom time -72% but wall +1.3% / +1.7% | rejected |
| level-dependent `Fapply` block width | -0.7% 2D (noise), +2.6% 3D (regression) | rejected, reverted |

## The actionable gap

Every deck in the repo sets `pre_smooth=4 post_smooth=4`
(`tests/ElasticSoftVoid/input:168-169`, `input_rod_and_tube_2d:188-189`,
`input_rod_and_tube_3d:157-158`, `input:139-140`). Task 3.3 measured 2/2 on
A100 in July (external wall -16.56%, MLMG solve -22.54%) and recommended it;
the recommendation never reached an input file. This task reproduces it locally
and extends validation from `max_step=2` to the hard rod-and-tube deck and the
production deck, with zero convergence failures on any of them.

Not applied here on purpose: 2/2 needs 27-41% more MLMG iterations, the longest
run in this task was 125 steps, and production runs to `stop_time = 6.5_s`. A
long-horizon run would close that gap. One-line change per deck when it does.

## Why coarsening reduction was tempting and is wrong

`mcl=1` was the best 2D number of the session (-9.4% at a *better* final
residual). It fails outright on `input_rod_and_tube_2d` —
`amrex::Abort::0::MLMG failed` at the first elastic solve — which is the
high-contrast failure mode `MLMG_HIGH_CONTRAST_FINDINGS.md` documents and which
that deck's `max_coarsening_level=2` stability recipe exists to prevent. It is
also +11.6% in 3D. `mcl=0` does not converge at all. mcl=2 is already tuned and
sits one step from a cliff.

Conversely, *raising* `mcl` is monotonically worse (+12% to +47%) for at most
one saved iteration — the clearest confirmation yet that added MG depth is a net
loss on GPU here, because each level's tiny-grid launches cost more than the
algorithmic benefit.

## Method

Nsight Compute is unusable on this workstation (`ERR_NVGPUCTRPERM`, needs root).
Register/spill/occupancy data therefore come from `cuobjdump -res-usage`, which
needs no privileges — `scripts/res_usage.sh`. Production instantiation is
`Elastic<Sym::Major>` (`Sym::Major == 1`, `src/Set/Matrix4.H:12`):

| kernel | 2D regs / stack | 3D regs / stack |
| --- | --- | --- |
| `Fapply` | 94 / 0 B | 252 / 288 B |
| `Diagonal` | 126 / 240 B | 255 / 1408 B |

2D `Fapply` has no register problem at all, which means historical tuning
measured on the 2D case could never have shown a register win.

## Closed without experiment

Host-device transfers and async copies (already nil: 268 transfers, ~0 GB,
0.12 ms); `Matrix4` loop unrolling and branch divergence (already hand-unrolled,
`Matrix4_Major.H:541-620`); occupancy/register tuning (measured dead end, twice);
`__launch_bounds__` (swept in task 3.2); fast math (already default); FP64
tensor cores (absent on sm_86); multi-stream overlap (one box per level, so
AMReX's stream rotation never advances); sparse formats (operator is
matrix-free).

## Reversal of a prior ranking

The 20260726 addendum ranked face-flux kernel fission second. This task argues
against it on its own data: splitting `Fapply` doubles 11,808 launches to
~23,616, and at the measured ~12 us mean pre-launch gap that adds ~140 ms to a
776 ms window — plausibly more than the `Matrix4` traffic it saves. On a
launch-latency-bound problem, fission is the wrong direction. It may still pay
on A100 where the launch/compute ratio differs; measure there first.

Remaining scoped work, unchanged in value: dropping the unread `m_ddw_mf` face
components on the psi path (a **memory-footprint** fix — measured `ddw`
ghost-exchange is 0.0% of kernel time, but the trim is 75% off the largest
array), shared-memory tiling of `U`, and CUDA graphs over the V-cycle — the only
candidate that attacks the 354 ms of idle rather than the 421.8 ms of busy.
