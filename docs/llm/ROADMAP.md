# chamber-gpu -> Optimized GPU Build: Transition Plan

**Branch:** `chamber-gpu` (ALAMO / AMReX `Integrator::Flame`)
**Branch policy:** This GPU work lives permanently on `chamber-gpu` and is **never merged to `master`** — the branch itself is the deliverable. "Done" means the branch meets its Phase 5 definition-of-done, not that it lands on master.
**Scope:** Take the current GPU port from "runs, but loses to one CPU core everywhere and has its safety nets off" to a build that is provably correct, launch/sync-optimized, and demonstrates a real win in the regime where one is achievable - or documents, with evidence, that no win exists for the target problem class.

---

## Guiding principles

These shape every phase. If a step violates one of these, stop.

1. **Correctness and honest measurement come before optimization.** You cannot optimize what you cannot trust or measure. The current build has disabled tripwires and a single-core baseline; both must be fixed first.
2. **The two paths are separate problems.** The phase-field / AMR path is launch-latency-bound and recoverable by implementation work. The elastic MLMG solve is structurally GPU-hostile (sequential V-cycle, tiny nodal boxes, large per-thread tensor state) and requires a decision, not a port. Do not let elastic effort block phase-field progress, and do not assume the same fix applies to both.
3. **The lever is cells-per-kernel, not resolution.** Deep AMR raises launch count and loses; wide grids raise work-per-launch and win. 3D helps only because per-box work scales as b^3 vs b^2, and only the phase-field path benefits - the elastic tensor state (Matrix4<3> = 45 doubles/thread vs 10 in 2D) gets worse.
4. **Every gate has explicit numeric exit criteria.** "Seems faster" is not a gate. Each gate below states the measured condition required to proceed.
5. **Every optimization change re-runs the correctness gate.** Optimizations must be numerically transparent - verified by golden compare, not assumed exact. (The current report claims "Option B is exact" while deferring the only test that proves it; do not repeat that.)
6. **"No win exists for the target problem" is a valid, fundable outcome.** If the analysis lands there, the deliverable is the evidence and the CPU-elastic/GPU-phasefield (or CPU-only) recommendation, not a forced GPU result.

---

## Phase map (dependencies)

```
P0 Stabilize + Instrument ──► P1 Elastic divergence + disposition ──┐
        │                                                           │
        └──────────────────► P2 Phase-field launch/sync opt ────────┼──► P3 Regime scaling
                                                                    │       (2D->3D, multi-GPU)
                                                                    │            │
                                                                    └────────────┴──► P4 Framework
                                                                                       dispatch decision
                                                                                            │
                                                                                            ▼
                                                                                       P5 Harden + lock
```

P1 and P2 can run in parallel once P0 passes (different code paths). P3 needs both. P4 needs P3's data to know whether generalization is even worth it.

---

## Phase 0 - Stabilize and Instrument

**Purpose:** Establish a trustworthy, measurable base. Nothing downstream is meaningful until this passes.

### Steps

- **0.1 Fix the known correctness defects.**
  - `Util/BMP.H`: `m_max[k] = std::min(...)` -> `std::max`; it currently pins `m_max` at `INT_MIN`.
  - `Util/Util.{H,cpp}`: device `Util::Abort` is a `// todo!` no-op, so base `Solid::DW/DDW` "not implemented" aborts silently return `Zero()` on device. Route device aborts to `amrex::Abort()`/`__trap()` (the `Assert` path already does this - make all of them consistent).
  - Restore NaN detection in a device-safe form: write a per-cell error flag into a scratch fab inside the kernel, reduce it once after the kernel, abort on the host. This replaces the deleted in-kernel `Util::Abort` NaN guards in `Flame.cpp` / `FullFeedback` / `Homogenize` without host-only calls on device.
- **0.2 Audit and fix the host-loop IC/BC landmines.** `BC/Operator/Elastic/Expression.H` and `IC/{PNG,PSRead,Trig,Laminate,Expression(vector)}.H` were demoted to `LoopConcurrentOnCpu`, which writes into device-arena fabs from the host - safe only under the managed arena. For each: either device it properly (as `Constant` BC was) or guard it so a pure-device-arena run cannot silently reach it. Document which ICs/BCs are currently GPU-safe (today: BMP IC + Constant BC only).
- **0.3 Build the golden-compare harness.** Generate a CPU reference plotfile for the chamber 20 MPa case (no-fast-math, IEEE) and a tolerance comparison tool (HDF5/plotfile diff). Define pass tolerance: <= 1e-12 relative for no-fast-math (FP reassociation only); a separate, looser, documented tolerance for fast-math.
- **0.4 Define the build matrix.** Two canonical builds, clearly named: `no-fast-math` (correctness / bit-compare gate, `--fmad=false` or conservative mode) and `fast-math` (performance). All correctness claims use the former; all perf headline numbers state which build produced them.
- **0.5 Replace the rigged baseline.** Stand up a **fully-subscribed CPU baseline** (16+ MPI ranks on a NOVA node), not 1 core. The current 1-GPU-vs-1-core comparison flatters the elastic "win" and understates the phase-field loss. Record CPU wall/step at the target resolutions as the real number to beat.
- **0.6 Extend the analysis suite with occupancy metrics.** Add Nsight Compute (`ncu`) capture for the hot kernels: achieved occupancy, registers/thread, SM efficiency, achieved DRAM bandwidth. The existing suite measures launches and sync fraction (good) but never occupancy - that gap hid the register-pressure ceiling.

### Information gathered
Launch count/step, kernel-duration histogram, `cudaStreamSynchronize` fraction (already have); **plus** occupancy + registers/thread for `Flame::Advance`, `Fapply`, `Diagonal`, `Newton::prepareForSolve`; CPU node baseline wall/step; golden-compare residuals.

### Report R0 - "Baseline of record"
A single page: corrected-build status, GPU-safe IC/BC matrix, golden-compare result at coarse res, multi-core CPU baseline numbers, and the occupancy/register table for the four hot kernels. This is the reference all later reports compare against.

### Gate G0 (must pass to proceed)
- [ ] `no-fast-math` build links and runs; golden compare passes at coarse res (<= 1e-12).
- [ ] BMP, device-abort, and NaN-detection fixes committed; no host-loop write into device-arena fabs reachable in the canonical input.
- [ ] Multi-core CPU baseline recorded.
- [ ] Occupancy/registers captured for the four hot kernels.

If G0 fails, do not optimize - the numbers will be untrustworthy.

---

## Phase 1 - Elastic divergence + path disposition

**REOPENED 2026-06-22.** The original D1 verdict below (CPU-resident) was
overridden by an explicit user directive: elastic must run on GPU. Treat this
phase as back in progress, not gated-closed — see `docs/llm/CURRENT.md` and
`VERSIONS.md` for the current investigation status. The steps/decision-tree
below are preserved as the original methodology; the new investigation is
narrower (a specific bottom-MG-level operator-evaluation defect, not the
general fast-math/conditioning question 1.1-1.2 already answered) and is
tracked outside this file for now (Claude memory: `gpu-elastic-d1-reversed`,
`gpu-elastic-nondeterminism`) pending a formal task-folder writeup.

**Purpose:** Resolve the high-res MLMG divergence and decide, with data, whether the elastic solve lives on the GPU, the CPU, or a hybrid. This is a decision phase.

### Steps

- **1.1 Run the no-fast-math GPU elastic at 2048**^**2** - the exact case that SIGABRTs under fast-math. Record convergence and MLMG iteration count vs the CPU on the identical system.
- **1.2 Branch on the result** (see D1). If it converges, the fast-math x conditioning hypothesis is confirmed and the fix is the build-flag split from 0.4. If it still diverges, it is **not** the flag - open an operator/conditioning investigation (psi formulation, `psi_floor`, masked-operator conditioning at thin features) before any further device-elastic perf work.
- **1.3 Occupancy/register profile** of `Fapply` / `Diagonal` / `Newton` with `maxrregcount` removed. Quantify the Eigen + `Matrix4` per-thread state cost.
- **1.4 Fair elastic benchmark:** GPU elastic vs the fully-subscribed CPU MLMG from 0.5, at the resolution where the elastic solve actually matters (>= 1024^2 finest), counting only the solve regions.

### Decision tree D1 - Elastic path disposition

```
D1
├─ no-fast-math GPU elastic converges at 2048^2?
│   ├─ NO  -> root-cause the operator/conditioning (1.2b).
│   │         Freeze device-elastic perf work. Elastic stays CPU-RESIDENT
│   │         until convergence is restored. Re-enter D1 after fix.
│   └─ YES -> continue
├─ achieved occupancy of Fapply/Diagonal (ncu)?
│   ├─ < 25%  -> register-bound, weak GPU fit; bias toward CPU/hybrid
│   └─ >= 25% -> viable on GPU
├─ GPU elastic wall vs N-rank CPU MLMG (fair baseline, >=1024^2)?
│   ├─ GPU slower than CPU node          -> CPU-RESIDENT elastic
│   │                                       (GPU phase-field + explicit batched
│   │                                        H2D/D2H around the every-50-steps solve)
│   ├─ within ~1.5x of CPU node          -> HYBRID: keep CPU-resident now,
│   │                                       re-evaluate after P3 (3D may shift it)
│   └─ GPU faster than CPU node          -> GPU-RESIDENT elastic
└─ Output: elastic disposition = {GPU | CPU | HYBRID}, with the numbers behind it.
```

### Report R1 - "Elastic path disposition"
Convergence result, root-cause (if applicable), occupancy/register findings, fair benchmark table, and the D1 verdict with confidence. State explicitly whether the "1.65x faster" coarse result survives a multi-core baseline.

### Gate G1
- [ ] Divergence resolved (converges) or root-caused (real bug identified).
- [ ] Elastic disposition decided and recorded.

---

## Phase 2 - Phase-field launch/sync optimization

**Purpose:** Recover the phase-field path. Apply the levers in order of leverage; measure after each; re-run the golden compare after each.

### Steps (apply one at a time, measure between each)

- **2.1 Box-sizing sweep (config, highest leverage, zero code).** Sweep `amr.blocking_factor`, `amr.max_grid_size`, `amr.grid_eff` (toward 0.9), `amr.regrid_int`, `amr.nsubsteps`, and AMR depth. Since GPU tiling is off by default and AMReX launches at the Box level, launch count ~ sum over levels of (boxes x 2^level). This step sets the ceiling for everything after it.
- **2.2 Reduction/sync removal.** Restructure `Flame::Integrate` to accumulate one `ReduceData` across all boxes (ideally all levels) and pull `.value()` once per step, not per box. Move the T/mdot/L min/max/sum diagnostics in `TimeStepComplete` to every `thermo.int` and batch them.
- **2.3 Pass fusion (vertical).** Fuse the three `TagCellsForRefinement` `ParallelFor`s into one. Fuse `UpdateModel`/`Advance` passes where a field is produced then immediately consumed, so intermediates stay in registers.
- **2.4 Level fusion for residual small-box launches.** For levels that still have many small boxes after 2.1, switch the hot MFIter+`ParallelFor` loops to AMReX's MultiFab-level `ParallelFor` (one launch per level when safe) or `ParallelFor(Tag)` (`AMReX_TagParallelFor.H`, one kernel over an irregular box set), gated on `MultiFab::isFusingCandidate()`. **Verify the signature works with node-centered `Set::Vector`/`Set::Scalar` fabs and the `Set::Field` wrapper against the AMReX source before committing.**
- **2.5 Field packing.** Combine co-evolving fields (eta, temp, ...) into one multi-component MultiFab so a single `FillBoundary` covers them, cutting ghost-exchange launches proportionally.
- **2.6 (optional) CUDA Graphs** on the fixed-topology per-step chain between regrids; re-capture on regrid; do not graph the MLMG cycle. Verify AMReX's current graph entry points before relying on this.

### Per-lever measurement protocol
After each sub-step, capture: launches/step, kernel-duration histogram, `cudaStreamSynchronize` fraction, wall/step vs the CPU node, **and** the golden compare. Record the delta attributable to that lever. Stop pursuing a lever when its marginal wall/step improvement is < 5 percent.

### Decision tree D2 - Optimization stop/continue

```
D2 (evaluated after each lever)
├─ golden compare still passes?
│   ├─ NO  -> the lever broke numerics. Revert, isolate, fix, re-test.
│   └─ YES -> continue
├─ launches/step still >> CPU-competitive AND kernel avg still ~ launch latency (~3-5 us)?
│   ├─ YES -> more launch-bound headroom: proceed to next lever
│   └─ NO  -> launch-bound largely cured; check the real target:
│       ├─ wall/step competitive with multi-core CPU node on a saturating config?
│       │   ├─ YES -> Phase-field path DONE for this regime -> G2
│       │   └─ NO  -> remaining loss is compute/occupancy or memory-bandwidth,
│       │             not launch overhead -> this is a REGIME problem -> go to P3,
│       │             do not keep grinding launch levers.
```

### Report R2 - "Phase-field optimization attribution"
A table: lever -> (launches/step before/after, avg kernel duration before/after, sync fraction before/after, wall/step vs CPU node, golden-compare residual). Plus the D2 verdict: is the path competitive, or is the residual gap a regime problem for P3?

### Gate G2
- [ ] Phase-field path competitive with (or beating) the multi-core CPU node on a saturating config, OR D2 has concluded the residual is a regime problem.
- [ ] Golden compare preserved through every lever.

---

## Phase 3 - Regime scaling (2D -> 3D, single -> multi-GPU)

**Purpose:** Find where a real win lives, if it exists. This is where the b^3 scaling and multi-GPU come in.

### Steps

- **3.1 3D readiness audit.** Confirm the 3D Eigen device paths (`.determinant()`, `.inverse()`, quaternion ops) compile and run; measure their occupancy/register cost. Quantify the `Matrix4<3>` (45 doubles/thread) elastic penalty - this is where the elastic disposition from D1 may need revisiting.
- **3.2 Memory budgeting.** Compute bytes/node (model fab at `Matrix4<3>` ~ 360 bytes/node, plus eta/phi/temp/displacement-vector fields) and the largest grid that fits on the target device. The A1000's 8 GB will not reach the saturating regime; budget for NOVA A100/H100.
- **3.3 3D box strategy.** Big base grid + shallow AMR + large `blocking_factor` (e.g. 32 -> 32768 cells/box) to land in the saturating regime that 2D could only approximate. Do **not** carry the deep-AMR/tiny-box pattern into 3D - it is as bad or worse there.
- **3.4 Multi-GPU.** 1 rank/GPU; characterize halo-exchange overhead; run strong- and weak-scaling studies on NOVA (A100 and H200 binaries already exist in the build scripts).
- **3.5 Crossover hunt.** Sweep problem size and GPU count to locate the point where the GPU build beats the CPU node, and by how much.

### Decision tree D3 - Regime / win-existence

```
D3
├─ Does a saturating 3D problem fit in target-GPU memory?
│   ├─ NO on A1000 -> move to NOVA A100/H100 (expected) -> continue
│   └─ YES -> continue
├─ At the largest fitting problem, phase-field kernels saturate (occupancy high,
│   avg kernel duration >> launch latency)?
│   ├─ NO  -> still latency/occupancy-bound even in 3D -> revisit P2 levers
│   │         and elastic disposition; the problem may be intrinsically too small
│   └─ YES -> continue
├─ GPU (single) wall vs CPU node at this size?
│   ├─ GPU faster      -> WIN at single-GPU; quantify, then test multi-GPU scaling
│   ├─ GPU within ~1.5x -> test multi-GPU; win likely emerges with scale
│   └─ GPU slower      -> single-GPU loses; only multi-GPU weak-scaling can justify
├─ Multi-GPU weak scaling efficient (>~70%) up to target node count?
│   ├─ YES -> WIN at scale -> G3 (win documented)
│   └─ NO  -> communication-bound; document the ceiling
└─ Output: {WIN @ single | WIN @ scale | NO WIN for target class}
            -- all three are valid, reported outcomes.
```

### Report R3 - "Crossover analysis"
The problem-size x hardware grid, with wall/step GPU vs CPU node at each point, occupancy at each point, and the crossover (or its absence). State the regime in which the GPU build is the right tool, and the regime in which CPU remains correct. This report drives the project's go/no-go on full GPU investment.

### Gate G3
- [ ] A net win demonstrated in a defined regime, **or** a documented, evidence-backed conclusion that no win exists for the target problem class (with the CPU-elastic/GPU-phasefield or CPU-only recommendation that follows).

---

## Phase 4 - Framework dispatch decision and generalization

**Purpose:** Resolve the whole-codebase implication. The de-virtualization of the shared `Solid`/`BC`/`Operator` hierarchies and the public-ification of integrator internals are global changes validated on one integrator via a de-facto fork (`alamo_gpu.cc` + hand-curated Makefile source list). Decide deliberately.

### Steps

- **4.1 CPU-semantics regression.** Run the full CPU test suite for **all** integrators built from `chamber-gpu` to confirm removing `virtual` from `Solid` did not change any integrator that relied on runtime `Solid*` dispatch. (The current report asserts this is "absent from the hot path" but shows no audit.)
- **4.2 Dispatch decision** (see D4): commit the hierarchies to static dispatch framework-wide, or isolate the device-model behind `ALAMO_GPU` / a separate device-model type so the shared CPU build is untouched for other integrators.
- **4.3 De-fork.** Replace the hand-curated CUDA source list and `alamo_gpu.cc` with a principled per-integrator GPU-clean policy (which integrators are nvcc-clean, tracked in the build system, not by a static file list).
- **4.4 Encapsulation convention.** Standardize the "copy POD params into a local, capture by value" pattern (currently ad hoc shadowing in `Flame.cpp`) into a documented convention, instead of making integrator internals `public` case by case.

### Decision tree D4 - Framework dispatch

```
D4
├─ Did P3 conclude there is a usable GPU win for the target class?
│   ├─ NO  -> Do NOT generalize. Isolate the de-virtualization behind ALAMO_GPU /
│   │         device-model type so the CPU build is unaffected. Ship Flame-GPU
│   │         (or CPU-elastic/GPU-phasefield) as a contained capability. STOP.
│   └─ YES -> continue
├─ Do other integrators need GPU within the planning horizon?
│   ├─ NO  -> isolate now, generalize later when needed
│   └─ YES -> commit to framework-wide static dispatch (CRTP), budget the full
│             port: de-virtualize + test every integrator together, with the
│             CPU regression suite (4.1) as the gate.
└─ Output: {ISOLATE device-model} or {COMMIT framework-wide static dispatch}
```

### Report R4 - "Framework dispatch strategy and generalization roadmap"
The CPU-regression result, the D4 verdict, and - if COMMIT - the per-integrator port order, effort estimate, and risk.

### Gate G4
- [ ] CPU regression suite green across all integrators on the branch.
- [ ] Dispatch decision ratified with the group (Runnels).

---

## Phase 5 - Hardening

**Purpose:** Lock in correctness and prevent regressions. This work stays on
`chamber-gpu` permanently; the goal is a trustworthy, self-contained branch, not
a merge to master.

### Steps
- **5.1** Re-arm correctness gates in CI: golden bit/tolerance compare on every push to `chamber-gpu` (no-fast-math), and a NaN-flag assertion build.
- **5.2** Performance-regression tracking: launches/step, avg kernel duration, sync fraction, wall/step vs CPU node, recorded per commit on the canonical case.
- **5.3** Documentation: the GPU-safe IC/BC matrix, the build matrix (no-fast-math vs fast-math and when to use each), and the elastic disposition decision.
- **5.4** Branch definition-of-done (below).

### Branch definition-of-done (stays on `chamber-gpu`; never merged to master)
- [ ] Golden compare passes (no-fast-math) on coarse and at least one saturating config.
- [ ] No host-loop writes into device-arena fabs reachable from any documented-supported input.
- [ ] Device aborts/NaN detection active.
- [ ] CPU regression suite green for all integrators.
- [ ] R3 win (or documented no-win + recommendation) recorded on the branch.

---

## Master pivot tree (when to change strategy)

```
At any gate:
├─ Correctness gate fails and cannot be fixed within the phase budget
│      -> halt optimization; correctness is non-negotiable.
├─ D1 = CPU-RESIDENT elastic (HISTORICAL — overridden 2026-06-22 by explicit
│      user directive; elastic device-residency is back in scope. Do not
│      route effort off device-elastic on the basis of this branch anymore.)
│      -> [superseded] redirect P2/P3 effort to the GPU phase-field + CPU
│         elastic split; this is the most likely "best realistic" architecture.
├─ D2 concludes residual gap is a regime problem
│      -> stop grinding launch levers; jump to P3.
├─ D3 = NO WIN for target class
│      -> deliverable becomes the evidence + CPU-only or hybrid recommendation;
│         D4 -> ISOLATE; do not generalize. This is a successful, fundable result,
│         not a failure.
└─ D4 = ISOLATE
       -> contain the device-model; ship Flame-GPU as a capability without
          touching the shared CPU build.
```

---

## Benchmark strategy addendum (2026-06-22 external review)

Findings from an external review of `benchmark/PHASE3_R3_crossover.md` and
this roadmap. These amend, not replace, the standing metric set and report
template below.

**Multi-GPU is on the back burner — not a near-term goal (explicit user
directive, 2026-06-24).** 2 GPUs is 3.5x-13.3x *slower* than 1 GPU at every
tested size (128^3/256^3/512^3, both NOVA sweeps) — a confirmed red flag, not
"poor scaling." This was previously frozen pending an `nsys` trace
(breaking down kernel time, CUDA API time, MPI time, `FillBoundary`/halo
cost, H2D/D2H/D2D copy volume — GPU oversubscription is ruled out;
non-GPU-aware MPI, bad box distribution, excessive halo-exchange, and
host-staged copies remain open hypotheses), but that diagnosis work is now
deprioritized rather than the next task. Revisit only if multi-GPU becomes a
near-term goal again.

**Three named benchmark tracks, not one blended number.** Future benchmark
reports should separate:
- `flame_only_3d` — elastic disabled, phase-field ceiling (what
  `PHASE3_R3_crossover.md` currently measures).
- `elastic_solve_only` — MLMG/Newton solve timing and convergence in
  isolation (what `PHASE1_ELASTIC_DISPOSITION.md` and
  `docs/gpu_elastic_device_port_plan.md` measure).
- `full_chamber` — end-to-end coupled run, elastic frequency stated
  explicitly (`elastic.interval`). This is the number that matters for the
  actual deliverable and is not yet measured anywhere.

**Amdahl's-law caveat — do not conflate flame-only speedup with full-chamber
speedup.** If elastic stays CPU-resident and is ~27% of CPU wall time (the
device-port plan's own estimate), a 13.1x flame-only speedup caps the
*full-chamber* speedup at roughly:
```
S = 1 / (0.73/13.1 + 0.27) ≈ 3.1x
```
If elastic's share grows to 50% of wall time, the same flame-only speedup
caps full-chamber speedup at roughly 1.9x. Every report that quotes a
flame-only number must state this cap alongside it, not let the flame-only
number stand in for the deliverable's actual speedup.

**Headline comparison numbers must come from fixed-step runs.** Time-limit
estimates (`1800s / steps_completed`) are acceptable for scouting but
understate steady-state throughput, especially at low step counts. Final
comparison tables use: a controlled `max_step`, output/plotting disabled
unless I/O is specifically under test, a discarded warmup window, and the
median of repeated runs — not a single time-boxed SLURM job.

**Profiling tool caveats:**
- `nvidia-smi` utilization/power timelines are coarse evidence only. They
  cannot substitute for `ncu` achieved occupancy, registers/thread, memory
  throughput, or warp-stall-reason data.
- `strace -f -c` aggregates syscall time across threads, so its percentages
  are not a clean wall-clock fraction. `perf stat` on the GPU binary mostly
  profiles host-side activity, not the device. Both are supporting evidence,
  not decisive.

**Benchmark data requirements (for any run intended to produce a headline
number):**
| Category | Data needed |
|---|---|
| Reproducibility | git commit + dirty-diff state, input file hash, exact command, module/AMReX/CUDA/compiler versions, GPU model, CPU node spec, MPI rank count, memory request |
| Wall timing | TinyProfiler region times, steady-state s/step, startup excluded, output on/off stated |
| GPU behavior | `nsys` CUDA API summary, kernel summary, memcpy summary, idle gaps, launches/step, sync fraction |
| Kernel quality | `ncu` achieved occupancy, registers/thread, SM efficiency, DRAM bandwidth, warp-stall reasons for `Flame::Advance`, `Fapply`, `Diagonal`, `Newton::prepareForSolve` |
| CPU baseline | fully-subscribed CPU-node timing, rank sweep, memory usage |
| Correctness | no-fast-math field diffs (not just step 1), thermo time-series diffs, pressure/mdot/volume/burn-area histories, NaN flags, elastic residual/iteration histories when elastic is enabled |

---

## Testing and reporting cadence (cross-cutting)

**Every phase produces:** (a) an information-gathering step with a fixed metric set, (b) a report, (c) a gate with numeric exit criteria, (d) a decision tree where the path forks.

**Standing metric set** (capture identically every time, for comparability):
- launches/step (`nsys` CUDA API summary)
- kernel-duration histogram + average (`nsys` kernel summary)
- `cudaStreamSynchronize` fraction of CUDA-API time
- achieved occupancy + registers/thread for the hot kernels (`ncu`)
- wall/step: GPU vs **multi-core CPU node**
- golden-compare residual (no-fast-math)

**Report template** (use for R0-R4):
1. Objective (one line)
2. Method (configs, builds, hardware)
3. Data (the standing metric set + anything phase-specific)
4. Decision-tree verdict
5. Recommendation + confidence (state it as a percentage; flag anything inferred rather than measured)
6. Open risks / what would change the verdict

---

## Risk register (initial)

| Risk | Phase | Mitigation |
|---|---|---|
| Divergence is a real operator bug, not fast-math | P1 | 1.2 root-cause before any device-elastic perf work |
| Elastic never competitive on GPU even in 3D (45-double/thread state) | P1/P3 | D1/D3 default to CPU-resident elastic; that split may be the best architecture anyway |
| Single A1000 too small to reach saturating regime | P3 | Budget NOVA A100/H100 early; do not draw regime conclusions from the A1000 |
| De-virtualization silently changed a CPU integrator | P4 | 4.1 CPU regression across all integrators before relying on the de-virtualized build |
| Optimization silently breaks numerics | all | golden compare after every lever (principle 5) |
| "No win" outcome treated as failure | P3 | principle 6: it is a valid, evidence-backed deliverable |

---

## One-line summary

Re-establish correctness and an honest multi-core baseline (P0); resolve the elastic divergence and decide GPU-vs-CPU elastic on data (P1); recover the phase-field path by cutting launches and syncs in leverage order, re-verifying numerics each step (P2); find the regime where a win exists - large 3D, shallow AMR, multi-GPU - or prove it does not (P3); then make the framework dispatch decision deliberately rather than by fork (P4) and harden (P5). At every fork, the decision trees route on measured thresholds, and "no win for this problem class" is an acceptable, useful result.
