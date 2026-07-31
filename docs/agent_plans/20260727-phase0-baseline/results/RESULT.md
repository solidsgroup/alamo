# RESULT — phase0-baseline

Campaign: `docs/agent_plans/20260727-gpu-memory-strategy/PLAN.md` §6.
Branch: `chamber-gpu-mem`. Started 2026-07-27.

Status: **IN PROGRESS** — local legs L1-L4 done; L5 (NOVA batch) composed and
awaiting the pre-submit checkpoint; N1-N5 and the decision gate pending.

---

## §L1 NVTX verification — DONE 2026-07-27

**No source work is needed. The campaign's "highest-leverage QOL investment" is
already paid for.**

Vendored AMReX 26.06 `Src/Base/AMReX_TinyProfiler.cpp` wraps every `BL_PROFILE`
region in an NVTX range whenever `AMREX_USE_CUDA` is defined:

- `:20-27` — includes `nvtx3/nvtx3.hpp`, `nvtx3/nvToolsExt.h`, or `nvToolsExt.h`
- `:134` — `nvtxRangePush(fname.c_str())` in `TinyProfiler::start()`
- `:211` — `nvtxRangePop()` in `TinyProfiler::stop()`

There is **no NVTX-specific build flag**. The requirement is `TINY_PROFILE=TRUE`,
which `configure --profile` supplies (`configure:533` →
`--enable-tiny-profile=yes`). `bin/alamo_gpu-2d-profile-cuda86-g++` already has
it.

Verified empirically, not just by source read: a 2-step 2D capture
(`nsys profile -t cuda,nvtx`, Nsight Systems 2026.1.3 from
`.local/nsight/opt/nvidia/nsight-systems/2026.1.3/bin/nsys`) yields a populated
`nvtx_sum` report. Top ranges by total time:

| Range | Instances | Total (ns) |
|---|---:|---:|
| `Integrator::Evolve` | 1 | 3,078,961,763 |
| `Integrator::Flame::TimeStepBegin` | 2 | 2,816,400,109 |
| `Integrator::Base::Mechanics::TimeStepBegin` | 2 | 2,613,574,528 |
| `MLMG::solve()` | 2 | 2,574,584,540 |
| `MLMG::oneIter()` | 48 | 2,515,630,123 |
| `MLMG::mgVcycle()` | 528 | 2,141,060,303 |
| `Operator::Fsmooth()` | 6,912 | 1,902,881,128 |
| `Operator::Elastic::Fapply()` | 16,304 | 1,440,144,876 |
| `Integrator::TimeStep` | 30 | 811,125,981 |
| `FillBoundary_nowait()` | 43,682 | 189,002,630 |

**Consequence for F10.** The profile build is the *diagnostic* binary, not the
timing binary: 16,304 NVTX push/pops for `Fapply` alone in two steps, and
`nova_flame_gpu.slurm`'s `MODE=bench` additionally sets
`tiny_profiler.device_synchronize_around_region=1`, which serializes exactly the
asynchrony this campaign is trying to measure. **F10 wall-time must come from the
non-profile binary**, or the headline number measures the profiler. Encoded in
the L5 batch as two separate runs.

## §L2 Step-loop region coverage — DONE 2026-07-27

`BL_PROFILE` counts: `Integrator.cpp` 19, `Operator/Elastic.cpp` 14,
`Flame.cpp` 8, `Newton.H` 2.

Named and usable for F1 today: `Integrator::Evolve` → `Integrator::TimeStep` →
`Flame::TimeStepBegin` / `Flame::Advance` / `Flame::TimeStepComplete` /
`Integrator::IntegrateVariables` / `Flame::Integrate` / `Flame::Regrid`, plus
the whole MLMG stack (`solve` → `oneIter` → `mgVcycle` → `mgVcycle_bottom` /
`miniCycle` / `mgFcycle`) and the elastic kernels (`Fapply`, `Fsmooth`), plus
`FillPatch` / `FillBoundary`.

**Gap, but not a blocking one:** `Solver/Nonlocal/Newton.H` carries only 2
regions, so Newton iteration structure — relinearization, line-search
backtracks, the per-iteration `norm0` reductions of campaign §9 / backlog 3.I —
is not separately attributable in F1; it aggregates into `MLMG::solve()` and its
children. Adding regions there is a **tier-2 source edit** and is therefore *not*
done under this tier-1 task. Deferred to Phase 2, which is where Newton sync is
the subject; it is not needed to read the Phase 0 baseline.

No proposal to add `BL_PROFILE` is being made now, so the L2 checkpoint does not
fire.

## §L3 Footprint budget — DONE 2026-07-27 (local anchor; NOVA refines)

### Measured, local A1000, 2D deck `input`, elastic every step

Method: three runs at `amr.n_cell = 64 / 128 / 256`, `max_grid_size=32`,
`max_step=2`, `amrex.the_arena_init_size=64 MiB` (small on purpose, so the arena
grows on demand instead of masking the high-water). Peak sampled at 5 Hz from
`nvidia-smi --query-compute-apps=pid,used_memory`.

| base `n_cell` | base cells | peak used (MiB) |
|---:|---:|---:|
| 64 | 4,096 | 298 |
| 128 | 16,384 | 390 |
| 256 | 65,536 | 730 |

Two-point slopes: 64→128 gives 7.85 KB per base-level cell; 128→256 gives
7.25 KB. Consistent. Linear fit intercept ≈ **267 MiB fixed overhead** — CUDA
context plus the arena floor, not field data.

Converting to resident cells: the verbose run of the same deck reports per-level
cell counts 4,096 / 16,384 / 34,240 / 31,936 = 86,656 total, i.e. an AMR
multiplier of ≈21× over base at `max_level=3` in 2D. So

> **≈ 350-380 B per resident cell**, elastic enabled, FP64, this deck.

That independently corroborates the 512 B/node elastic-enabled figure hard-coded
in `benchmark/phase3_memory_budget.py` (whose breakdown attributes ~360 B/node to
the `Matrix4<3>` elastic model fab alone), and it means the script's model is
conservative rather than optimistic.

### Budget against device memory

Analytic, `phase3_memory_budget.py --bytes-per-node 512 --max-level 3
--refine-fill 0.1` (AMR multiplier 59.4× in 3D, ghost factor 1.3, 85% headroom):

| device | mem | largest wide-shallow base grid | est. |
|---|---|---|---|
| A1000 | 8 GB | 64 × 64 × 32 | 4.83 GB |
| A100-40 | 40 GB | 64 × 64 × 32 | 4.83 GB |
| A100-80 / H100-80 | 80 GB | 128 × 128 × 64 | 38.61 GB |
| H200 | 141 GB | 128 × 128 × 64 | 38.61 GB |

**Verdicts:**

1. **FP32 storage is optional, not structural** — for chamber. The 2D campaign
   deck peaks at 298 MiB on an 8 GB card, ~90% of which is CUDA context. Nothing
   in Phases 1-3 has to be built around halving the footprint. It stays a
   bandwidth lever for Phase 3 (§10: FP32 halves traffic), not a survival
   requirement. Campaign §8's "if the budget clears device memory, FP32 becomes
   structural" branch is **not taken**.
2. **Local correctness work is not memory-limited.** Every 2D campaign deck, and
   3D up to roughly 64 × 64 × 32 base with 3 levels, fits in 8 GB. The A1000
   constraint on this campaign is HMM and the 50 W cap, not capacity.
3. **`multicomponent/FMA` is where this arithmetic bites**, exactly as campaign
   §8 predicts: species count multiplies field count against the same 512 B/node
   baseline. The measured 360 B/resident-cell anchor is the reusable number.

Caveats: all local, FP64, 2D-anchored; the 3D column is the script's first-order
model, not measurement. The authoritative 3D bytes/node comes from N1.

## §L4 `rod_and_tube_step2` GPU golden reference — RESOLVED (re-recorded)

Campaign §6 lists this as a Phase 1 prerequisite. It was real, it is now fixed,
and the fix **tightens** the gate rather than loosening it.

### The failure

`baseline_suite.py check --profiles=gpu_strict` against
`bin/alamo_gpu-2d-nofast-cuda86-g++`: three cases ok, `rod_and_tube_step2` FAIL
on all four traction columns.

```
trac_xhi_x: abs=1.384000e+02 rel=2.265232e-03 FAIL
trac_xhi_y: abs=4.585000e-01 rel=2.349704e-02 FAIL
trac_yhi_x: abs=6.512000e-01 rel=2.642621e-01 FAIL
trac_yhi_y: abs=1.256000e+02 rel=2.100370e-03 FAIL
```
`gpu_fast` failed with byte-identical deltas.

### Why re-recording is correct here, and not test-weakening

For every one of the four columns, `|fresh_gpu − gpu_ref|` equals
`|cpu_ref − gpu_ref|` **exactly**:

| column | cpu ref | old gpu ref | \|cpu−gpuref\| | reported \|fresh−ref\| |
|---|---:|---:|---:|---:|
| `trac_xhi_x` | −61097.5000 | −60959.1000 | 138.4000 | 138.4000 |
| `trac_xhi_y` | 19.5131 | 19.0546 | 0.4585 | 0.4585 |
| `trac_yhi_x` | 2.4642 | 1.8130 | 0.6512 | 0.6512 |
| `trac_yhi_y` | −59799.0000 | −59673.4000 | 125.6000 | 125.6000 |

The fresh GPU result **has moved onto the CPU reference**. The stored GPU
references were recorded before the GPU elastic path was repaired (the
elastic-void recovery / `Matrix4` device-annotation lineage); they pin a
disagreement that no longer exists. Keeping them would gate against known-wrong
values.

Independent confirmation that the CPU side is the fixed point:
`check --profiles=cpu --case rod_and_tube_step2` passes today.

### Action taken

- Old references preserved as evidence in
  `results/stale_refs_backup/{gpu_strict,gpu_fast}.json`.
- Re-recorded with `baseline_suite.py record --profiles=gpu_strict,gpu_fast
  --case rod_and_tube_step2`.
- Post-record GPU-vs-CPU agreement: **worst relative deviation 4.06e-6**
  (`trac_yhi_x`), down from 2.64e-1. No tolerance was touched.
- Full re-check: **8/8 green** — all four cases × `gpu_strict` + `gpu_fast`.

**The gpu_strict leg is now trustworthy. Campaign §6's prerequisite and NOTES.md
N4 are both discharged.** Phase 1 may rely on it.

## §L5 NOVA batch — SUBMITTED 2026-07-27

Pre-submit checkpoint cleared by the user. Build job `11772112` COMPLETED
(00:15:37), producing all four binaries on sm_80:
`alamo_gpu-{2d,3d}-cuda80-g++` (timing) and `alamo_gpu-{2d,3d}-profile-cuda80-g++`
(diagnostic). These are the first non-profile NOVA binaries this repo has had —
`build_alamo_nova.sh` previously built only `--profile`.

### Capture jobs

| Job | Deck | Legs | Purpose |
|---|---|---|---|
| 11772151 | `input_copy` | env timing flip nsys ncu | **Primary.** Production-condition deck |
| 11772152 | `input` | env timing flip nsys ncu | Prior 2D baseline, kept for continuity |
| 11772153 | `input_3d_centre_bore_128_a2` | env timing flip nsys ncu | 3D (Fapply-compute-bound regime) |
| 11772154 | `input_copy` | smooth | 4/4 vs 2/2 A/B, 800 steps |

`11772151` and its retry `11772159` both landed on `nova21-gpu-10` and died at
CUDA init; `11772161` is the `input_copy` capture that ran. Outcomes in §N1.

### Two design corrections made before submitting

**1. The original horizons measured nothing.** Elastic is ~95% of GPU wall, and
every production deck fires it on an interval after a `tstart` delay:

| deck | dt | tstart | interval | 1st solve | 2nd solve |
|---|---|---|---|---|---|
| `input_copy` | 2.5e-4 | 0.01 s | 40 | step ~40 | step ~80 |
| `input` | 1.0e-4 | 0.5e-4 | 50 | step ~50 | step ~100 |
| `input_3d_centre_bore_128_a2` | 1.0e-4 | 0.5e-4 | 50 | step ~50 | step ~100 |

The first draft used `MAX_STEP=10` and `TRACE_STEP=3`, which contain **zero
elastic solves** on every one of these decks. That baseline would have been
confidently wrong — it would have profiled the phase-field/thermal steps and
missed the entire dominant cost. Horizons are now per-deck (90/110 steps) and
sized to span two solves at **production cadence**. Deliberately *not* fixed by
overriding `elastic.interval=1`, which measures a configuration nobody runs.

**2. `input_copy` would have aborted at `InitData` on NOVA.** Its eta IC is
`blur5_rod_and_tube.bmp`, untracked locally and absent on the remote. `push` now
ships `*.bmp`.

## §L6 Smoothing A/B — the experiment two prior tasks declined to run

User report 2026-07-27: 2/2 pre/post smoothing causes MLMG divergence in
production runs, against perf studies recommending it. Both studies located, and
**both predict this**:

- `docs/agent_plans/20260721-fapply-runtime-optimization/results/RESULT.md`
  retains 2/2 on measured gains (2D wall −12.79%, 3D −21.05%; A100 wall −16.56%,
  MLMG solve −22.54%, FApply −23.30%) and states **twice** that it is "not a
  longer-evolution stability claim". The tell: the 3D final residual is
  `9.72181e-09` against a `1e-08` gate — clearing convergence by 3% on a
  **two-step** deck.
- `docs/agent_plans/20260727-gpu-optimization-investigation/results/RESULT.md`
  R1 measures 2/2 costing **5-41% more MLMG iterations** and explicitly refuses
  to apply it: *"the longest validation run here was 125 steps, and production
  runs to `stop_time = 6.5_s`… that margin is exactly what gets consumed."*

Mechanism: 2/2 buys wall by doing less work per V-cycle and repays it in
iteration count. On short horizons the saving exceeds the repayment. `input_copy`
runs ~6,000 steps (`stop_time=1.5_s` at `dt=2.5e-4`) with ~150 elastic solves at
`interval=40`, and sets `tol_abs=1e-8` — the exact gate 2/2 cleared by 3%. The
divergence is the predicted outcome of an unvalidated extrapolation, not a
contradiction of the measurements.

Job `11772154` ran both arms at 800 steps (vs the 125-step prior maximum),
reporting **MLMG iteration counts and survival**, wall secondary.

### Verdict — 2/2 REFUTED at production horizon

Node `nova22-amp-5`, deck `input_copy`, A100-SXM4-80GB. The two arms differ in
`elastic.solver.pre_smooth` / `post_smooth` and nothing else
(`benchmark/phase0_capture.slurm`, leg `smooth`).

| arm | rc | wall_s | elastic solves | MLMG iters mean | max |
|---|---:|---:|---:|---:|---:|
| 4/4 | 0 | 431.264 | 37 | 157.41 | 289 |
| 2/2 | 6 | 6.482 | 0 | — | — |

2/2 does not cost iterations. **It diverges on the first elastic solve.** The
residual grows monotonically by ~4.4× per iteration from the outset:

```
MLMG: Iteration   7 Fine resid/bnorm = 11.48471813
MLMG: Iteration  20 Fine resid/bnorm = 2394726743
MLMG: Iteration  37 Fine resid/bnorm = 1.735000728e+20
MLMG: Failing to converge after 37 iterations. resid, resid/bnorm = 5.388825455e+27, 1.735000728e+20
amrex::Abort::0::MLMG failing so lets stop here !!!
```

Dead at 6.5 s of an expected 431 s, before completing one solve. Evidence:
`smooth/{4x4,2x2}/summary.txt` and `smooth/2x2/run.log`.

**This is stronger than the mechanism proposed above.** The pre-run reasoning —
2/2 repays wall in iteration count and eventually consumes a 3% convergence
margin — predicts gradual degradation. What happens is immediate divergence: the
first production-cadence solve never converges at all. The 20260727 study's
measured "5-41% more MLMG iterations" was taken on solves that still converged;
that regime does not extend to `input_copy` at `interval=40`.

**Consequences.**

1. Keep 4/4. No deck should be changed; none currently is (all decks already set
   4/4, so nothing shipped and nothing needs reverting).
2. `docs/llm/PLAN.md:26-36` next-task item 1 (3.3) reads "Retain
   configuration-only 2/2 pre/post smoothing". That recommendation is now
   refuted at production horizon and must be annotated — not deleted, since its
   own two-step measurements stand and it explicitly bounded itself to that
   horizon. **Out of scope for this tier-1 task; raised for the campaign.**
3. Both prior tasks were right to refuse to apply it. The open question was
   never whether 2/2 is faster on a short deck; it is, measurably. It is that
   the short deck cannot see the failure mode.
4. The 4/4 arm's own profile is worth carrying into Phase 1: solves 1-3 cost
   124/289/273 iterations before settling into a stable ~138-162 band. The
   startup transient is ~2× the steady state.

### Unrelated defect found in the same deck

`input_copy:175` sets `elastic.solver.bottom_solver = bigcstab` — transposed
letters. `src/Solver/Nonlocal/Linear.H:286-288` matches only `cg`, `bicgstab`,
`smoother`, with **no `else` branch**, so an unrecognized value is silently
discarded. Inert here by luck: `MLLinOp::getDefaultBottomSolver()`
(`AMReX_MLLinOp.H:264`) returns `bicgstab` and Elastic does not override it, so
the deck gets what it meant. The hazard is the next one — `smoother` is the
documented fix for the Mode-B high-contrast failure, and a typo there would
silently no-op while the deck reads as though it were set. Logged, not fixed.

## Decision gate

Decision gate **answered 2026-07-27** (campaign PLAN §17): proceed, targeting
Flame + Elastic. Chamber/Ballistic needs no work; Hydro deferred to a follow-on
port (campaign PLAN §18).

## §N1 Baseline capture — PARTIAL

All four jobs reached `COMPLETED` at the job level; two legs failed inside them.

| Job | Deck | env | timing | flip | nsys | ncu |
|---|---|---|---|---|---|---|
| 11772152 | `input` (2D) | ok | ok | ok | **DEAD** | 3/3 |
| 11772153 | `input_3d_centre_bore_128_a2` | ok | ok | ok | **DEAD** | **0/3 OOM** |
| 11772161 | `input_copy` (2D) | ok | ok | ok | **DEAD** | 2/3 |

Job `11772151` and `11772159` are the `nova21-gpu-10` casualties (see §L5); both
re-ran elsewhere. **F1, F2, F3, F8, F9 are not captured for any deck** — those
five come from nsys, and nsys produced no `trace.nsys-rep` anywhere. F4-F7 exist
for 2D only. Root causes and fixes in §H below.

### F10 — wall/step, non-profile binary (NOT yet admissible)

| deck | steps | managed | device | Δ |
|---|---:|---:|---:|---:|
| `input` (2D) | 110 | 12.803 s / 0.1164 per step | 13.278 s / 0.1207 per step | device **+3.7%** |
| `input_copy` (2D) | 90 | 33.581 s / 0.3731 per step | 31.726 s / 0.3525 per step | device **−5.5%** |
| `input_3d_centre_bore_128_a2` | 110 | 12.772 s / 0.1161 per step | 12.453 s / 0.1132 per step | device **−2.5%** |

**These are not yet a result.** Campaign PLAN and `docs/llm/PLAN.md:18-22` both
require a startup calibration alongside wall/step; none was run. Every number
above includes CUDA context creation and `InitData`, which at 12-13 s totals is a
large fraction. Run `max_step=1` per deck before quoting any of this.

Sign of the arena effect flips by deck (+3.7% / −5.5% / −2.5%), which is what an
uncalibrated ~2-6% measurement looks like. Do not conclude "device arena wins" or
"loses" from this table.

**Open question for the re-capture:** the 3D deck costs the same per step as 2D
`input` (0.1132 vs 0.1207 s). If elastic is ~95% of GPU wall, a 3D case should
not match a 2D one. Either the 3D deck is far smaller than assumed or its elastic
solves are not landing inside the 110-step horizon. Verify before the gap table
uses any 3D number.

## §N2 Device-arena flip inventory — DONE. The crash set is empty.

`amrex.the_arena_is_managed=0` plus `amrex.abort_on_out_of_gpu_memory=1`,
`the_arena_init_size` pinned to 0.5 × 81920 MiB = 42949672960 B per rank.

| deck | rc | failures.txt |
|---|---:|---|
| `input` (2D) | 0 | empty (0 lines) |
| `input_copy` (2D) | 0 | empty (0 lines) |
| `input_3d_centre_bore_128_a2` | 0 | empty (0 lines) |

Grep pattern was `illegal|CUDA error|Abort|out of memory|segmentation` over the
full run log. Zero hits on any deck, all three running their full horizon.

**The finding is that there is nothing to triage.** The crash-triage
classification table (timestep loop / diagnostic / checkpoint / debug leftover)
has no rows to fill, because no host-pointer dereference or arena exhaustion
occurred. Phase 1 was budgeted around triaging a device-arena crash set that does
not exist on these decks at these horizons.

> **T1 claim WITHDRAWN 2026-07-28.** This section originally read "T1 passes."
> It does not follow. `ext/AMReX-Codes/amrex/Src/Base/AMReX_Arena.cpp:59` sets
> `the_arena_is_managed = false` by default, so the device arm is what the code
> does anyway and the *managed* arm was the synthetic one. Per campaign PLAN
> v2 §3.3, a clean run cannot distinguish "no managed allocations" from
> "managed allocations that happen to work", and `The_Managed_Arena()` persists
> as a separate pool regardless of this flag. T1 under v2 requires a static
> call-site inventory plus an allocation trace attributed by pool. The
> measurement below stands as a measurement; the conclusion drawn from it does
> not. First cut at the inventory: NOTES.md N9.

Scope of the claim, stated honestly:

- Three decks, 90-110 steps, **single rank**, one GPU, 80 GB. Multi-rank
  device-arena behaviour is untested here.
- 0.5 × 80 GB per rank is a generous arena. It does not probe the capacity edge,
  which is the other thing `abort_on_out_of_gpu_memory` is for.
- Horizons span two elastic solves but no checkpoint/restart. The known restart
  temp-fab crash (`elastic-void-recovery` memory) is outside this window, so
  "checkpoint" as a crash class remains genuinely unprobed.

## §N3 Nsight Compute — PARTIAL

Counter permission **granted** on NOVA (`ncu counters permitted.` from a real
capture, not a `nvidia-smi` inference). This retires the `ERR_NVGPUCTRPERM`
risk that Step N3's VERIFY was gating on, and means `benchmark/res_usage.sh` is
not needed as a substitute here.

Captured, 2D only: `Operator::Elastic::Fapply()` and `Operator::Fsmooth()` on
both `input` and `input_copy`, with `.ncu-rep` + details CSV. `MLMG::mgVcycle()`
succeeded on `input`, OOM-killed on `input_copy`. All three 3D ranges OOM-killed.
F5/F6/F7 therefore exist for 2D and not for 3D — and 3D is the
Fapply-compute-bound regime, so the roofline that sets the T7 band is exactly the
one still missing.

## §N4 GPU-aware MPI — DONE. Inactive on NOVA. Campaign §9 is debt.

Compute node `nova21-gpu-12`, inside the allocation:

```
mca:mpi:base:param:mpi_built_with_cuda_support:value:false
mca:mpi:base:param:mpi_built_with_cuda_support:source:default
```

Combined with Phase 0.5 §2 (same verdict on kermit), GPU-aware MPI is inactive on
**both** machines. Campaign §9's "device-buffer Allreduce" cannot be implemented
as written and must ship as acknowledged debt, not as a planned optimization.
Any device buffer handed to `MPI_Allreduce` will be staged through the host
regardless of what the call site looks like.

### The probe was reading the wrong MPI

The env leg loaded `MPI_MOD=openmpi/4.1.8-lctrfpx` and reported *that*
`ompi_info`. `ldd bin/alamo_gpu-2d-cuda80-g++` resolves `libmpi.so.40` to
**openmpi-5.0.8** (`.../openmpi-5.0.8-mybw4ysihjo77xmjcknse5my7koichix/lib`) via
RPATH — a different implementation behind the same soname. The nsys crash
backtrace independently confirms 5.0.8 and UCX 1.18.1 are what actually load at
runtime.

Re-checked the linked build directly, 2026-07-28:

```
mca:opal:base:param:opal_built_with_cuda_support:value:false
mca:opal:base:param:opal_cuda_support:value:false
mca:accelerator:null:version:"component:5.0.8"
```

**Verdict unchanged** — 5.0.8 is also built without CUDA support, and its
`accelerator` framework resolves to the `null` component. But the evidence in
`11772152/env/inventory.txt` cites a library the run never used, and a future
capture would have re-recorded the same wrong provenance. Probe fixed in §H.

Both builds export `libmpi.so.40`, so the 4.1.8 module could in principle shadow
the 5.0.8 the binary was built against. Tested on NOVA 2026-07-28 with
`MPI_MOD` loaded exactly as the job does: `ldd` still resolves 5.0.8, so RPATH
wins and the mismatch is inert. `MPI_MOD` is therefore left at 4.1.8 — changing
it without evidence risks breaking a working capture — but the probe now prints
a `MISMATCH` line so the next reader is not misled the way this one was.

## §N5 Gap table and T7 threshold — BLOCKED

Cannot be populated. N5 needs F2/F3 (bus traffic, nsys) and the 3D roofline
(F5, ncu), and neither exists. Step N5's CHECK — "no gap-table cell is empty;
T7 has a number and a derivation" — is unreachable until the re-capture lands.
No partial gap table is written here, because a half-populated one invites
exactly the assumption-driven scoping the task exists to prevent.

## §H Harness defects found and fixed (2026-07-28)

Three, all in `benchmark/phase0_capture.slurm`. No `src/` change; tier 1 holds.

### H1 — nsys traced MPI and died before the first timestep

`-t cuda,nvtx,osrt,mpi`. nsys 2024.6.2's `libToolsInjection64.so` segfaults
inside `ompi_mpi_init` → `opal_common_ucx_mca_register` against Open MPI 5.0.8 /
UCX 1.18.1:

```
libToolsInjection64.so(+0x2b445a)
libopen-pal.so.80(opal_common_ucx_mca_register+0x4f)
libmpi.so.40(ompi_mpi_init+0x98)
```

`nsys rc=137` on all three jobs, no `trace.nsys-rep` written. **Fix:** dropped
`mpi` from `-t`. Every run in this harness is `-n1`, so MPI tracing was never
load-bearing.

### H2 — the harness then misreported H1 as a version quirk

With no trace file, all eight `nsys stats` calls failed, and the loop printed
`stats MISSING in this nsys version: <report>` eight times — which reads as
"this nsys is old" rather than "the capture is dead". `stats.log` held the real
answer (`ERROR: Specified input file ... does not exist`) but nothing surfaced
it. **Fix:** guard on `trace.nsys-rep` existence; on absence, print an explicit
`CAPTURE FAILED`, name the affected figures, and drop a `CAPTURE_FAILED`
marker file instead of running the loop.

### H3 — ncu profiled 3 launches, then let the app run 80 more steps

`ncu --kill` defaults to `0` (verified against NOVA's ncu 2025.1.1:
`--kill arg (=0)`). Without it, ncu takes its `--launch-count 3` matches at the
first elastic solve and then follows the app to the end of the horizon with NVTX
collection still attached, retaining range state for every `BL_PROFILE` push/pop
(~16,304 `Fapply` ranges per 2D step, per §L1). `MLMGmgVcycle.log` shows the app
reaching **step 81** before the SLURM OOM killer fired:

```
==ERROR== The application returned an error code (9).
error: Detected 1 oom_kill event in StepId=11772161.7. Some of the step tasks have been OOM Killed.
```

The two ranges that survived still wrote 173 MB reports for 3 profiled launches.
**Fix:** `--kill yes` (spelling verified on the NOVA module: valid choices are
`on|off`, `yes|no`, `1|0`, `true|false`), plus `--mem` raised 64G → 128G for
headroom. Note rc is now expected to be nonzero on this leg — ncu terminates the
target by design — so the leg is judged by the `.ncu-rep`, not by rc. The echo
line comments record this so the next reader does not "fix" it back.

**Not yet re-run.** These fixes are unvalidated on hardware; the next submission
is the test.

### H4 — OSRT, not MPI or `srun`, is the injection trigger

The H1 fix was incomplete: the 2026-07-28 re-capture still aborted on all three
decks after `mpi` was removed from the trace list. Two focused NOVA jobs on
2026-07-31 isolated the trigger:

| Job | Arm | Result |
|---|---|---|
| 11825224 | `/bin/true`, direct `nsys` | PASS |
| 11825224 | ALAMO, direct `nsys` | FAIL, same MPI-init exception |
| 11825224 | ALAMO, `srun nsys` | FAIL, same MPI-init exception |
| 11825250 | `-t cuda` | PASS |
| 11825250 | `-t nvtx` | PASS |
| 11825250 | `-t osrt` | **FAIL** |
| 11825250 | `-t cuda,nvtx` | PASS |
| 11825250 | `-t cuda,nvtx,osrt` | **FAIL** |

The failing arms throw `Expected shared object name, found a path delimiter`
from `libToolsInjection64.so` while Open MPI 5.0.8 registers UCX. The passing
`cuda,nvtx` arm initializes CUDA, completes one timestep, and emits a nonempty
report. Older Nsight Systems 2024.5, 2023.4, and 2023.3 also pass with
`cuda,nvtx`, confirming that changing toolkits is unnecessary.

**Fix:** the Phase 0 nsys leg now traces `cuda,nvtx` only. CUDA tracing retains
the CUDA API calls needed for transfer and synchronization accounting; OSRT is
not required by F1/F2/F3/F8/F9. Hardware validation still requires a full-horizon
re-capture.

---

## Local preview (indicative only, not admissible evidence)

From the same 2-step L1 trace, `cuda_gpu_mem_size_sum`:

| Direction | Total (MB) | Count | Avg (MB) |
|---|---:|---:|---:|
| Device-to-Device | 65,235.4 | 5,291 | 12.33 |
| Host-to-Device | 61.3 | 8,328 | 0.007 |
| Device-to-Host | 22.8 | 382 | 0.060 |

Two leads for N1 to confirm or kill on real hardware:

1. **65 GB of device-to-device copies in two steps**, 5,291 copies averaging
   12.3 MB. This is not bus traffic, so it does not threaten T3 — but it is
   memory traffic on a bandwidth-bound code, and campaign §10's framing
   ("optimize bytes moved") applies to it directly. Worth attributing to MLMG
   level transfers vs `FillPatch` in F2.
2. **8,328 host-to-device copies averaging 7 KB.** Small but numerous, and H2D
   inside a step loop is what T3 is about. Needs per-step attribution against
   the NVTX ranges before it means anything — the count includes `InitData` and
   the plotfile write.

Both figures are from a shared, power-capped A1000 under a profile build. They
set the questions for N1; they answer nothing.
