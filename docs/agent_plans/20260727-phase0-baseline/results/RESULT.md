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

Job `11772154` runs both arms at 800 steps (vs the 125-step prior maximum),
reporting **MLMG iteration counts and survival**, wall secondary.

### Unrelated defect found in the same deck

`input_copy:175` sets `elastic.solver.bottom_solver = bigcstab` — transposed
letters. `src/Solver/Nonlocal/Linear.H:286-288` matches only `cg`, `bicgstab`,
`smoother`, with **no `else` branch**, so an unrecognized value is silently
discarded. Inert here by luck: `MLLinOp::getDefaultBottomSolver()`
(`AMReX_MLLinOp.H:264`) returns `bicgstab` and Elastic does not override it, so
the deck gets what it meant. The hazard is the next one — `smoother` is the
documented fix for the Mode-B high-contrast failure, and a typo there would
silently no-op while the deck reads as though it were set. Logged, not fixed.

## §N1-N5, decision gate

Decision gate **answered 2026-07-27** (campaign PLAN §17): proceed, targeting
Flame + Elastic. Chamber/Ballistic needs no work; Hydro deferred to a follow-on
port (campaign PLAN §18).

N1-N5 pending job completion.

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
