# RESULT — GPU optimization investigation (2026-07-27)

Branch `chamber-gpu`, HEAD `19132c02c`. RTX A1000 (sm_86, 16 SM, 8 GB).
Gates green before and after. **No source change is retained**; `src/` is clean
at HEAD. The one retained finding is a configuration change to input decks,
described below and deliberately *not* applied — see "Recommendation".

## Headline

The elastic solve is **launch-latency bound, not compute bound**, and the single
largest available win on this hardware is a solver configuration that four
separate decks agree on:

| deck | pre/post 4/4 | pre/post 2/2 | delta | MLMG iters | failures |
|---|---|---|---|---|---|
| ElasticSoftVoid 2D conservative | 1.367 s | 1.223 s | **-10.5%** | 82 -> 86 | 0 |
| ElasticSoftVoid 3D section | 7.044 s | 5.265 s | **-25.3%** | 72 -> 76 | 0 |
| `input_rod_and_tube_2d` (hard) | 2.18 s | 2.03 s | -6.9% | 64 -> 90 | 0 |
| production `input` | 12.11 s | 11.39 s | -5.9% | 79 -> 100 | 0 |

Every deck in the repo currently sets 4/4 (`tests/ElasticSoftVoid/input:168-169`,
`input_rod_and_tube_2d:188-189`, `input_rod_and_tube_3d:157-158`,
`input:139-140`). Task 3.3 measured this same configuration on A100 in July and
recommended retaining it (external wall -16.56%, MLMG solve -22.54%); that
recommendation was never propagated into any input file. This task independently
reproduces it locally and extends validation from `max_step=2` to the hard
rod-and-tube deck and the production deck.

## Method note

Nsight Compute is unusable here (`ERR_NVGPUCTRPERM`, needs root). Register,
spill and occupancy-ceiling data therefore come from `cuobjdump -res-usage`,
which needs no privileges — see `scripts/res_usage.sh`. Runtime evidence is
paired wall-clock A/B plus the MLMG iteration counts and relative residual from
`elastic.solver.verbose=2`, so that a wall win bought by weaker convergence is
visible rather than hidden (`scripts/mg_sweep.sh`).

## Baseline device resources, production `Elastic<Sym::Major>`

`Sym::Major == 1` (`src/Set/Matrix4.H:12`); both Flame and the benchmark use
`NeoHookeanPredeformed : NeoHookean : Solid<Set::Sym::Major>`.

| kernel | 2D regs / stack | 2D occ ceiling | 3D regs / stack | 3D occ ceiling |
|---|---|---|---|---|
| `Fapply`   |  94 / 0 B   | 45.4% | 252 / 288 B  | 16.9% |
| `Diagonal` | 126 / 240 B | 33.9% | 255 / 1408 B | 16.7% |
| `Energy`   |  52 / 0 B   | 82.0% | 120 / 0 B    | 35.5% |

2D `Fapply` has no register problem at all (94 regs, zero spill). The register
story is 3D-only — which matters because most historical tuning was measured on
the 2D case, where a register win cannot express itself.

## Measured and retained

**R1 — pre/post smoothing 4/4 -> 2/2.** Configuration only. Costs 5-41% more
MLMG iterations and returns 6-25% of wall on every deck tested, with zero
convergence failures. The gain is largest in 3D because each smoothing sweep is
a full `Fapply` + fused-`Fsmooth` launch pair, and 3D pays more per launch.

## Measured and rejected

**J1 — `max_coarsening_level` reduction (mcl=1).** Looked like the best result of
the session in 2D: **-9.4%** wall (1.239 s vs 1.367 s) at a *better* final
residual. It is nonetheless rejected on two independent grounds:

- It **fails outright** on `input_rod_and_tube_2d`: `amrex::Abort::0::MLMG
  failed` at the first elastic solve, Newton iteration 2. This is the exact
  high-contrast failure `MLMG_HIGH_CONTRAST_FINDINGS.md` warns about, and the
  deck documents `max_coarsening_level=2` as part of its validated stability
  recipe.
- It is **+11.6% slower in 3D** (7.864 s vs 7.044 s), with bottom-solve time
  doubled. The 2D win does not generalise.

`mcl=0` does not converge at all (200 iterations, `resid/bnorm = 3.40e-05`,
MPI_ABORT), so mcl=2 sits one step from a cliff and is already the tuned value.

**J2 — deeper coarsening.** Monotonically worse for zero convergence benefit:

| mcl | 2D wall | vs base | iters |
|---|---|---|---|
| 2 (current) | 1.367 s | — | 82 |
| 3 | 1.531 s | +12.0% | 81 |
| 4 | 1.743 s | +27.5% | 81 |
| 6 | 2.006 s | +46.7% | 81 |

Each added MG level buys at most one iteration and costs 12-16% of wall in
tiny-grid launch overhead. This is the clearest direct confirmation of the
launch-bound thesis in the whole investigation: **more algorithmic depth is a
net loss on GPU here.** Nobody should raise `mcl` on this hardware.

**J3 — `final_smooth` reduction.** Cuts bottom-solve time by up to 4x
(0.40 s -> 0.11 s at `final_smooth=2`) and still makes wall *worse*
(+1.3% at 4, +1.7% at 2): the saving is repaid in outer iterations. Reject.

**J4 — level-dependent `Fapply` block width.** The one occupancy item left
unmeasured by the 20260726 task. Implemented as a runtime dispatch (wide 256 for
boxes >= 8192 nodes, narrow 128 otherwise), verified live against the measured
grid histogram (fine tier 12,800 nodes triggers wide; the 4,352-node tier
correctly stays narrow). Result: **-0.7% in 2D (noise), +2.6% in 3D
(regression)**. Mechanism is visible in the static data: at 252 regs a
256-thread block fits exactly one block/SM (64,512 of 65,536 registers) while
128 threads fit two — identical occupancy, coarser scheduling granularity,
worse tail. 128 was already the right choice. Reverted.

**J5 — `bottom_solver`.** No-op: the decks already set
`bottom_solver = smoother` (`tests/ElasticSoftVoid/input:172`). Measured
bit-identical, as expected. Recorded because it also means the safety
precondition for touching `mcl` was already satisfied — mcl=1 failed for
genuinely algorithmic reasons, not for want of a smoother bottom.

## Closed by existing evidence — no experiment run

| candidate | why closed |
|---|---|
| Minimizing host-device transfers; async transfers | Already nil: 268 transfers, ~0 GB, 0.12 ms for the whole run. Nothing to recover. |
| Loop unrolling / branch divergence in `Matrix4` | Already done. `operator*(Matrix4,Matrix)`, `(Matrix4,Matrix3)` and `MulCol` are hand-unrolled to literal `data[]` indices (`Matrix4_Major.H:541-620`). The `uid` if-chain in `operator()` is not on the hot path. |
| Occupancy / register-pressure tuning | Measured dead end. Task 3.3 Step 4 cut `Fapply` 254->250 regs and 440->288 stack bytes for no material gain; J4 above reconfirms. At <= 9 blocks on 16 SMs, per-SM occupancy is irrelevant — most SMs hold no block at all. |
| `__launch_bounds__` | Already swept (task 3.2). `Fapply` spills hard at `min_blocks >= 2` (128-reg cap vs 252 live). Knob exists, default off. |
| Fast math | Already the default build; `--cuda-fp strict` is the opt-out used for bit-identity checks. |
| FP64 tensor cores (DMMA) | Absent on sm_86 consumer silicon; the 9x9 contraction is too small and irregular to map onto A100 DMMA. |
| Multi-stream overlap | No concurrency available: one box per level, so AMReX's per-MFIter stream rotation (`AMReX_MFIter.cpp:378`, `max_gpu_streams` default 4) never advances — 32,484 of 32,772 launches on one stream. |
| Sparse-matrix formats | Not applicable: the operator is matrix-free, assembled nowhere. |

## Argued against, deliberately not built

**Kernel fission of the face-flux branch.** This was ranked #2 by the 20260726
task, and this investigation's own data now argues against it. `Fapply` is
11,808 launches; splitting flux computation into a separate kernel doubles that
to ~23,616. At the measured ~12 us mean gap before a `Fapply` launch, the added
launches cost roughly 140 ms against a 776 ms window — plausibly more than the
`Matrix4` traffic the split saves. **On a launch-latency-bound problem, fission
is the wrong direction.** It may still pay on A100, where the launch/compute
ratio differs; it should be measured there first, not here.

## Scoped, not executed

**S1 — drop unread `m_ddw_mf` face components on the psi path.** `m_ddw_mf` is
allocated with `AMREX_SPACEDIM + 1` `Matrix4` components and 2 ghost nodes
(`Elastic.cpp:73-75`). The conservative branch reads components `1..dim`
(`Elastic.cpp:285-287`); the psi branch reads only component 0. Production runs
the psi branch — `Newton.H:377` gates conservative flux on `m_psi == nullptr`,
and chamber sims always set psi. In 3D that is 3 of 4 components, i.e. 1080 of
1440 bytes per node, written, coarsened and ghost-exchanged for nothing.

Correct framing: **this is a memory-footprint fix, not a speed fix.** Measured
`ddw` ghost-exchange is 4 launches and 0.1 ms — 0.0% of 2D kernel time. Its
value is the 75% cut to the largest array, which is what forces managed memory
and makes the 3D psi case (26 GB high-water on an 8 GB card) unmeasurable
locally. Implementation complication: `SetConservativeFaceFlux` is called after
`define()` (`Newton.H:353`), so allocation currently precedes the decision;
either defer allocation or reallocate on the flag.

**S2 — shared-memory tiling of `U` in `Fapply`.** Unchanged launch count, so
unlike fission it does not fight the launch-bound regime. Needs a custom launch;
AMReX `ParallelFor` gives no tiled block shape.

**S3 — CUDA graphs over the MLMG V-cycle.** The only candidate that attacks the
354 ms of idle directly rather than shaving the 421.8 ms of busy. AMReX carries
`AMReX_CudaGraph.H`. Largest scope of anything here; correspondingly the largest
ceiling.

## Recommendation

Apply 2/2 smoothing to the input decks — a one-line change per deck in
`tests/ElasticSoftVoid/input`, `input_rod_and_tube_2d`, `input_rod_and_tube_3d`
and `input`. **Not applied by this task**, for one reason: these are production
physics-configuration files, the longest validation run here was 125 steps, and
production runs to `stop_time = 6.5_s`. 2/2 consistently needs 27-41% more MLMG
iterations to reach the same tolerance; on a long horizon with a harder
conditioning state that margin is exactly what gets consumed. Task 3.3 hit the
same boundary and recorded the same limit. A long-horizon run on one geometry
would close it.

## Caveats

- All timings are RTX A1000 (sm_86). FP64 throughput there is ~1/64 of FP32
  versus ~1/2 on A100, so any FP64-throughput-bound conclusion may not transfer.
  The retained finding (R1) is a launch-count effect, which should transfer, and
  already agrees with task 3.3's A100 measurement.
- The production `input` validation was run with the user's uncommitted
  `Flame.{cpp,H}` chi-rename WIP compiled in. It does not touch the elastic
  solver, but the run is not from a clean tree.
- The `use_psi=1` 3D branch remains unmeasured locally (26 GB managed
  high-water on an 8 GB card); it belongs on NOVA.

## Evidence

- `artifacts/e1_mg_sweep/<arm>/` — per-rep wall, stdout, iteration counts
- `scripts/res_usage.sh` — static regs/stack/occupancy (no root required)
- `scripts/mg_sweep.sh` — paired config A/B with convergence reporting
