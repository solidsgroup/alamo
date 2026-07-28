# chamber-gpu-mem: Memory Strategy Migration Plan

**Version:** 2.0
**Branch:** `chamber-gpu-mem` (from `chamber-gpu`; Flame solver, ALAMO/AMReX)
**Status:** Phase 0 blocked. Target set revised. Decision gate reopened.
**Character:** Strategy and scope. Execution detail deferred to per-phase work items.

---

## 0. Changelog from v1.0

v1.0 was rejected by hostile review. The rejection was substantially correct. Changes:

| Change | Origin |
|---|---|
| **T0 wall-time non-regression added as blocking gate** | Omission: campaign could pass every target and ship slower production |
| **T1, T3, T4, T5, T6 redefined** | Metrics were structurally incapable of detecting their own pathologies |
| **False-pass analysis now mandatory per target** | Generalized fix for the above; the meta-defect |
| **Phase 1 split into 1a (mechanical) and 1b (lifetime redesign)** | Review correctly rejected "Phase 1 is mechanical" as a whole-phase claim |
| **T6 split; idle-fraction gate moved to Phase 3** | Phase 2 could not exit on a gate whose fix was scheduled for Phase 3 |
| **Ballistic scope contradiction escalated to explicit decision** | Plan declared Ballistic out of scope while requiring device-side pressure |
| **Sync inventory now mechanically derived, not authored** | Authored inventory missed the largest contributors |
| **`AbortIfDeviceError` measurement promoted to early Phase 0** | Cheap, removable, may reshape the idle profile everything else is scoped against |
| **Benchmark methodology section rewritten** | No repetitions, startup included, horizon shorter than the startup transient |
| **Reproducibility requirements added to capture** | Baseline captured from a 574-file dirty tree with no manifest or hashes |
| **Oracle integrity requirements added** | Main gate ran a managed-memory smoke while reporting sanitizer coverage |
| **Decision gate reopened** | Cost estimate moved materially; prior answer was against the old cost |

### The v1 root defect

v1 stated that arena pools recycle after warmup, then selected metrics that pooling defeats. The same error repeated across four targets:

| Target | v1 measured | Pathology lives in | Defeated by |
|---|---|---|---|
| T1 | "runs clean, managed off" | managed arena as a separate pool | pool separation |
| T3 | bytes per step | transfer count, blocking calls | small payloads |
| T4 | `cudaMalloc` calls | arena request count | pool recycling |
| T5 | high-water mark | per-step churn | pool recycling |

System-level symptoms were used to certify application-level properties. Pooling, batching, and caching sit between the layers and hide exactly what was claimed to be detected.

**Structural remedy, Section 3.1:** no target enters the set without a written false-pass scenario and a demonstration that the metric catches it.

---

## 1. Scope

### In scope

- Memory residency and arena discipline for the chamber timestep loop
- Allocation lifetime redesign where Phase 0 evidence shows per-step churn
- Elimination and reduction of per-step host/device synchronization
- Launch configuration and kernel shape, to the extent memory traffic drives it
- Multi-GPU *structural* decisions only (Section 11)
- A repeatable, reproducible profiling harness with a fixed figure set

### Out of scope

- Algorithmic changes to the Flame solver
- Elastic/MLMG solver algorithm work (separate thread)
- `multicomponent/FMA` porting (this plan generates patterns for it)
- Multi-GPU performance tuning, halo overlap, load balance

### Scope items pending decision

- **Ballistic** — declared out of scope in v1 while the design required device-side pressure. Unresolved contradiction. See Section 9.2.
- **Newton convergence policy** — check-every-N is not an algorithmic change; fixed-iteration-count is. See Section 9.3.

### Governing rule

> Device owns field data permanently. Host copies are transient, explicit, created for I/O, and destroyed immediately. No field crosses the bus inside a timestep.

Resolve ambiguity against this sentence. Note that it is a *goal*, not a gate — the gates are Section 3, and the v1 failure was assuming a plausible gate followed from a clear goal.

---

## 2. Blocking preconditions

None of Phases 1-3 may begin until all four clear.

| # | Precondition | Why blocking |
|---|---|---|
| P1 | Target set v2 complete with false-pass analysis (Section 3) | Better measurement of wrong targets is wasted work |
| P2 | Capture reproducibility fixed (Section 5.1) | Ledger keyed on commit SHA cannot identify a dirty-tree measurement |
| P3 | Oracle integrity fixed (Section 5.2) | Agent self-certification against an inadequate oracle already occurred |
| P4 | Decision gate re-run at revised cost (Section 6) | Prior answer was against a materially lower cost estimate |

---

## 3. Target definition, v2

### 3.1 Meta-gate: false-pass analysis

Every target carries a written false-pass scenario: a concrete way the pathology exists while the metric reports PASS. The target is admitted only when the metric demonstrably catches that scenario.

Where a known instance of the false-pass already exists in the tree, the metric is validated against it directly before the target is used for certification.

This gate applies to the target set itself and is checked before Phase 1.

### 3.2 Targets

| # | Target | Measurement | Gate at | Threshold |
|---|---|---|---|---|
| **T0** | **Production wall time does not regress** | Steady-state step time, production deck, repetitions with uncertainty | **Every phase exit** | No regression outside uncertainty |
| T1 | No managed allocations | Static call-site inventory + allocation trace by arena | 1a | Zero bytes from managed pool in steady state |
| T2 | No page-fault migration | nsys unified-memory counters | 1a | 0 events |
| T3a | Field data does not cross bus per step | Bytes per step, by direction | 2 | Scalars only |
| T3b | Transfer count bounded | Transfer count per step | 2 | Bounded, each justified |
| T3c | Blocking transfer calls bounded | Synchronizing API call count per step | 2 | Bounded, each justified |
| T4 | No unjustified per-step allocation | **Arena-level** alloc/free request count per step | 1b | Bounded, each individually justified |
| T5a | Device footprint stable | Arena high-water between regrids | 1a | Flat |
| T5b | Live allocation count stable | Arena live-allocation count per step | 1b | Flat |
| T5c | Host footprint stable | Host-side container growth per step | 1b | Bounded |
| T6a | Synchronizing operations bounded | Count per step vs mechanical inventory | 2 | Bounded, each justified |
| T6b | GPU idle fraction | Timeline idle inside step loop, coarse-NVTX build | **3** | Threshold set at Phase 0 |
| T7 | Kernels near their applicable roof | Per-kernel, classified by limiter | 3 | Per-class, set at Phase 0 |

### 3.3 False-pass analysis per target

**T0** — Passes on a small case while production regresses. *Catch:* production deck, production-length horizon, startup excluded, uncertainty estimated. Small-case timing is not admissible for T0.

**T1** — v1 form ("runs clean with device arena") fails three ways: `The_Managed_Arena()` persists as a separate pool when the default arena is device-only; the vendored AMReX already defaults `the_arena_is_managed=false`, so the v1 "flip" proves nothing about the code; and a successful run does not distinguish "no managed allocations" from "managed allocations that happen to work." *Catch:* static inventory of managed-arena call sites plus an allocation trace attributing bytes to pools.

*Definitional prerequisite:* "current state" must be pinned to one of AMReX defaults, benchmark `MODE=bench`, or production launch practice. v1 left this undefined and Phase 0 consequently manufactured a managed baseline that does not correspond to how the code is run.

**T2** — Meaningless if managed was never enabled. Report alongside the T1 definitional answer or not at all.

**T3a alone** — Thousands of tiny scalar transfers pass a bytes-per-step threshold. This pathology is already visible in the local trace. *Catch:* T3b and T3c. Count and blocking behavior carry their own thresholds.

*Multi-rank caveat:* T3 cannot pass while GPU-aware MPI is inactive, since ghost exchange and device-buffer collectives may stage field data through host. Single-rank measurement does not license a device-residency claim.

**T4** — v1 measured post-warmup `cudaMalloc`. AMReX arenas cache and reuse backing allocations, so unbounded application-level churn produces zero `cudaMalloc` calls. *Known instances in tree:* MLMG operator and solver constructed and destroyed per elastic solve; `FieldNorm0` allocating composite MultiFabs per call including line-search backtracks; per-step device error flag allocating a `DeviceScalar` with H2D and D2H round trip. *Catch:* instrument at the arena request layer, not the CUDA API. Validate the instrumentation against these three before use.

**T5a** — Same defeat as T4: pooling holds high-water flat under arbitrary churn. Retained as necessary-not-sufficient; a *growing* high-water is still a real signal. *Catch:* T5b.

**T5c** — Device-side metrics watch the wrong side of the machine for host-side growth. `Ballistic::Advance` appending to a host `std::vector` every step is unbounded growth invisible to every other target here.

**T6a vs T6b** — v1 gated Phase 2 on under-5% idle while stating that 2D is 45.6% idle and launch-latency bound, with launch work deferred to Phase 3. Unsatisfiable by construction. Split: T6a counts synchronizing operations, which Phase 2 controls; T6b measures idle fraction, gated at Phase 3, threshold derived from the measured launch-bound floor rather than assumed.

*Measurement integrity:* a build emitting 16,000+ NVTX ranges in a two-step run perturbs the launch gaps T6b measures. Two-build policy, Section 12.2.

**T7** — v1 assumed the HBM roof applies. Fapply history points at register pressure, spill, and occupancy limits. A kernel can sit far from the HBM roof and be correctly optimized. *Catch:* classify each top kernel by actual limiter first; assign the target appropriate to its class. "Near the roof" is not a universal target.

### 3.4 Explicitly not targets

Zero bus traffic is wrong. Three transfers are legitimate: plotfile and checkpoint writes (bulk, amortized, async, pinned); regrid redistribution (infrequent, device-to-device); convergence and diagnostic scalars (few, batched, counted under T3b/T3c).

The distinction is not whether data moves but whether it is scheduled, explicit, visible in source, counted, and outside the inner loop.

---

## 4. Phase overview

| Phase | Name | Character | Exit gates |
|---|---|---|---|
| 0.5 | Two-rank correctness probe | Correctness | Multi-rank matches single-rank |
| 0 | Baseline, inventory, target validation | Measurement | Gap table; T6b/T7 thresholds set; P1-P4 clear |
| 1a | Arena hygiene | Mechanical | T0, T1, T2, T5a |
| 1b | Allocation lifetime redesign | Refactor, evidence-scoped | T0, T4, T5b, T5c |
| 2 | Sync elimination | Design | T0, T3a-c, T6a |
| 3 | Launch config and kernel shape | Tuning | T0, T6b, T7, **re-check T4/T5** |
| 4 | Multi-GPU | Deferred | Out of scope here |

Ordering 1a to 1b to 2 to 3 is strict. T0 gates every phase exit, not just the last.

T4 and T5 are re-checked at Phase 3 exit: fusion work can reintroduce a step-loop allocation and silently undo Phase 1b.

---

## 5. Phase 0 — Baseline, inventory, target validation

### 5.1 Capture reproducibility (P2)

The existing capture records HEAD plus a dirty-file count. With 574 dirty files including `Flame.*`, a ledger row keyed on commit SHA cannot identify the measured source, and the branch cannot bisect or reproduce.

Requirements:

- Capture refuses to run on a dirty tree, **or** records a full diff plus per-file content hashes and a manifest of every copied artifact including untracked inputs
- Ledger key becomes `(commit_sha, tree_hash, case, n_ranks, device, build_config)`
- Results absent locally (`results/figures/`, metrics ledger) are produced or the run is not admissible

No measurement taken before this is fixed enters the ledger.

### 5.2 Oracle integrity (P3)

The main gate reported three green labels while its sanitizer leg forced `TIERS=1` — a managed-memory pre-elastic smoke, neither compute-sanitizer nor elastic. The load-bearing sanitizer is Tier 2. The golden leg defaults to CPU rather than `gpu_strict`.

This is the predicted failure mode, confirmed: an agent self-certified correctly against an inadequate oracle.

Requirements:

- Sanitizer leg runs Tier 2
- Golden leg defaults `gpu_strict`
- Label text matches what the leg actually executes
- **Coverage requirements are explicit:** pressure history over enough steps to diverge, at least one regrid, two ranks, at least one checkpoint/restart cycle

Harness coverage defines what agent-produced code can be trusted to be. It is designed deliberately, not inherited from whatever small case existed.

### 5.3 Early cheap measurement

**Before scoping Phase 2 or 3.**

`AbortIfDeviceError` costs two unconditional full-stream synchronizations per Flame level per step and is a debug facility. If it compiles out, the effort is near zero and the idle profile may change substantially.

Measure with it disabled first. Scoping Phase 2 against an idle profile dominated by a removable debug sync wastes the scoping.

### 5.4 Mechanical sync inventory

v1 described its authored inventory as confirmed. It missed the two largest contributors.

The Phase 0 deliverable is a **mechanically derived** inventory: enumerate synchronizing call sites from source, tag each with call frequency per step, and rank by measured cost. Authored lists are not accepted.

Known omissions from the v1 list, to be included and not treated as exhaustive:

- `AbortIfDeviceError` full-stream syncs, per level per step
- Explicit `streamSynchronizeAll()` in Newton and its line search
- Per-box reduction landings in `HasNonFinite`
- Elastic-operator lifetime synchronization
- Full-MultiFab copies and refluxing inside every `FieldNorm0`

These may outweigh the chamber scalar path. Phase 2 scope follows the inventory, not the reverse.

### 5.5 Footprint budget

```
resident = cells * components * (1 + ghost_overhead) * sizeof(Real) * n_multifab_copies
```

Summed across levels and solver working copies, compared against device capacity divided by ranks-per-GPU.

Device-arena success established on single-rank 80 GB A100 does not test the 40 GB budget. Test it explicitly.

If the budget clears capacity, FP32 storage becomes structural rather than optional, and that changes what Phases 1-3 build. On a bandwidth-bound code, storage precision is a memory-strategy decision: it halves footprint and halves traffic, the same lever twice.

### 5.6 Threshold derivation

Phase 0 sets, from measurement rather than assumption:

- **T0** baseline and uncertainty band
- **T6b** idle threshold, from the measured launch-bound floor
- **T7** per-class targets, after classifying top kernels by actual limiter

### 5.7 Deliverables

1. Reproducible baseline capture (5.1)
2. Fixed oracle with stated coverage (5.2)
3. `AbortIfDeviceError` disabled comparison (5.3)
4. Mechanical sync inventory (5.4)
5. Footprint budget including 40 GB case (5.5)
6. T0/T6b/T7 thresholds (5.6)
7. Target set v2 with false-pass analysis validated against known in-tree instances (3.1, 3.3)
8. Gap table
9. Revised cost estimate for the decision gate (Section 6)

---

## 6. Decision gate — reopened

v1's gate was marked answered by a scope clarification while its evidence items remained incomplete. That is a process failure. The substantive failure matters more:

**Phase 1 moved from "mechanical, low-risk" to "mechanical plus lifetime redesign."** Reaching genuinely allocation-free steady state requires operator lifetime restructuring and `FieldNorm0` rework — a real refactor with correctness risk. The campaign cost estimate moved materially.

The gate asked whether chamber optimization is worth doing *at the v1 cost*. That answer does not carry to the v2 cost.

Re-run it with:

- Revised Phase 1a + 1b + 2 estimate
- Ballistic scope decision resolved (9.2), since it changes Phase 2 cost
- Honest statement of whether chamber currently gates SRM paper throughput

If chamber-gpu already turns around the runs the paper needs, `multicomponent/FMA` is the better claim on the time. Readiness is not a reason.

---

## 7. Phase 1a — Arena hygiene

**Character:** Mechanical, low-risk. Unconditional.

| Data class | Arena | Rationale |
|---|---|---|
| MultiFab field data | `The_Arena()`, device, non-managed | Governing rule |
| Kernel temporaries | `The_Async_Arena()` | Stream-ordered free, no sync on destruction |
| I/O staging, reduction landing | `The_Pinned_Arena()` | DMA-direct; pageable D2H forces a synchronous staging copy |
| MPI buffers (non-GPU-aware path) | `The_Comms_Arena()` | Isolation from field pool |
| Managed | Nothing | Each exception justified in writing and recorded under T1 |

Configuration: device arena permanent; `the_arena_init_size` as a **fraction of device memory divided by ranks-per-GPU**, never a hardcoded byte count; `abort_on_out_of_gpu_memory=1`.

Parameter names and defaults move between AMReX versions. Verify against the vendored version. Note in particular that the vendored version already defaults `the_arena_is_managed=false`, which changes what this phase is actually doing relative to v1's description.

**Method:** flip per subsystem, one commit per subsystem. Porting everything then flipping yields simultaneous segfaults with no attribution.

**Crash triage — classify before fixing:**

| Crash site | Fix |
|---|---|
| Inside the timestep loop | Port the kernel |
| Diagnostic / screen output | Explicit staged copy through pinned |
| Checkpoint / plotfile | Explicit staged copy, async |
| Debug leftovers | Delete |

Classification is a required step, not advice. The default agent response to a segfault is to port the kernel, which is wrong for three of the four rows.

**Exit:** T0, T1, T2, T5a. Harness green including the 5.2 coverage requirements.

---

## 8. Phase 1b — Allocation lifetime redesign

**Character:** Refactor. Correctness risk. Evidence-scoped — content determined by the Phase 0 gap table and arena-level instrumentation.

This phase exists because T4 was redefined. Under the v1 metric it was invisible.

### Known candidates

| Site | Behavior | Direction |
|---|---|---|
| MLMG operator and solver | Constructed and destroyed per elastic solve | Hoist lifetime above the solve loop |
| `FieldNorm0` | Full composite MultiFabs per call, including line-search backtracks | Persistent scratch; remove copies and refluxing from the norm path |
| Per-step device error flag | `DeviceScalar` alloc, H2D, D2H, free each step | Persistent device-resident flag |
| `Ballistic` history | Host `std::vector` append per step, unbounded | Bounded or externalized (T5c) |

Each candidate is admitted on measured per-step request count, not on the list above.

**Risk:** lifetime changes are correctness-sensitive in ways arena hygiene is not. Hoisted operator lifetime interacts with regrid. Persistent scratch interacts with level count changes. Both need explicit regrid coverage in the harness.

**Exit:** T0, T4, T5b, T5c. Per-step allocation set bounded and each entry individually justified in writing.

---

## 9. Phase 2 — Sync elimination

**Character:** Design. Scoped by the Section 5.4 mechanical inventory, not by the chamber scalar path alone.

### 9.1 Scope follows the inventory

v1 scoped this phase against the chamber feedback loop. The mechanical inventory may rank `AbortIfDeviceError`, Newton line-search syncs, `HasNonFinite` reductions, and `FieldNorm0` above it. Scope after the inventory exists and after the 5.3 measurement.

### 9.2 Chamber feedback loop — scope contradiction, requires decision

The chamber model is a per-step scalar round trip: reduce regression rate over burning surface, mass flux, scalar ODE for pressure, pressure feeds Arrhenius mobility and burn rate. Host-side implementation drains the pipeline every step and appears as dead timeline, not as a slow kernel.

The v1 target — pressure resident on device, updated by a one-thread kernel — **contradicts the v1 scope declaration that Ballistic is out of scope.** `Ballistic::Advance` is host-only, appends to a host `std::vector`, and its pressure feeds host-side Propellant state and elastic traction.

Three options, all with costs:

| Option | Cost |
|---|---|
| Port Ballistic | In-scope expansion; changes its history and diagnostic contract |
| Duplicate the formula in Flame | Maintenance trap; two sources of truth for the pressure law |
| Keep host update | Fails the target; Phase 2 exits with a known per-step drain |

**Additional obstacle:** `ReduceData::devicePtr()` holds block partials. The finalized value comes from `.value()`, which lands on host. A true device-result reduction needs additional reduction machinery beyond what v1 assumed. This raises the cost of options 1 and 2.

**Decision required before Phase 2 scoping, and it feeds Section 6.**

### 9.3 Global reduction shape

The burning-surface reduction is global. Multi-rank, the scalar update is an `MPI_Allreduce`, not a local reduce.

Build as `device-local reduce -> device-buffer Allreduce -> device scalar update`. Single-rank the Allreduce is a no-op. Cost now is nil; correctness later is by construction.

This must be specified explicitly in any implementation spec. Device-local reduce is the dominant pattern in training data, works single-rank, and passes the harness. It is the predicted agent failure for this phase.

### 9.4 Newton convergence

If the Newton convergence test is a host-side comparison on a reduced residual, that is a sync per Newton iteration stacked on the per-step drain. Confirm whether the earlier CPU stall fix left a host-side branch in place — correct CPU design, serious GPU defect, survives review unremarked.

**Scope clarification:** check-every-N is *not* an algorithmic change provided it only ever adds iterations past convergence — the converged answer is identical to tolerance. Fixed-iteration-count *is* an algorithmic change and remains out of scope.

**Honest limit:** device-side convergence testing does not remove the host control dependency on whether to continue or terminate. Some sync survives. T6a is a bounded-and-justified target, not a zero target, for this reason.

**Exit:** T0, T3a, T3b, T3c, T6a. Idle fraction is *not* a Phase 2 gate.

---

## 10. Phase 3 — Launch configuration and kernel shape

**Character:** Tuning. Last, because the bottleneck relocates after Phase 2.

**Box size.** `max_grid_size` up from CPU-tuned values; 128+ typical. Confirm `TilingIfNotGPU()` throughout. *Constraint:* box count is decomposition granularity — choose against a target rank count or accept a Phase 4 retune.

**Fusion.** Fuse when it removes a field read/write cycle; split when registers spill. Verify with `-Xptxas -v` or `--resource-usage`.

**Occupancy limiter.** Identify the actual cap — registers, shared memory, block count — before adjusting.

**Framing.** Phase-field and stencil physics are memory-bound: optimize bytes moved; FP32 gains come from bandwidth; hardware selection is an HBM bandwidth question. But per T7, verify this per kernel rather than assuming it. Fapply history suggests some kernels are register- and occupancy-limited, and those need different treatment.

**Exit:** T0, T6b, T7, **plus re-check of T4, T5a, T5b**. Fusion can reintroduce step-loop allocation and silently undo Phase 1b.

---

## 11. Phase 4 — Multi-GPU (deferred)

Deferred: halo overlap and comm hiding; load balance and regrid distribution; scaling studies. All require an optimized single-GPU baseline to be measurable.

Carried forward: Allreduce-shaped reduction (9.3); rank-aware `init_size` (Section 7); box size against target rank count (Section 10); GPU-aware MPI confirmed active (Section 13).

Scope depends on hardware. Two cards in a workstation and a cluster allocation are different projects.

---

## 12. Benchmark methodology

v1 had no methodology section. The Phase 0 results correctly recorded themselves as inadmissible.

### 12.1 Timing

- **Repetitions with uncertainty estimate.** One managed run and one device run, unrepeated, is not a measurement.
- **Interleaved and randomized order** across configurations to absorb machine drift.
- **Startup excluded.** Calibration derived up front, not deferred to later analysis.
- **Horizon exceeds the startup transient.** Two elastic solves cannot characterize steady state when the first three are transient. Set the horizon from the measured transient length.
- **Production-length stability run** for T5, separate from timing runs.

### 12.2 Build policy

Two builds, both in the ledger:

| Build | NVTX | Used for |
|---|---|---|
| Coarse | Stage-level ranges only | T0, T6b, all timing |
| Fine | Full annotation | Diagnosis, F1 reading |

Never measure T0 or T6b on the fine build. 16,000+ ranges in a two-step run perturbs exactly the launch gaps being measured.

### 12.3 Deck matrix

Present: thermal-on, variable-pressure-on, static elastic. Missing, and required:

- Thermal-off
- Constant-pressure
- Checkpoint / restart cycle
- Per-field multi-rank comparison
- Multi-GPU device arena
- 40 GB capacity case
- Production-length memory stability

---

## 13. Phase 0.5 — Two-rank correctness probe

**Runs first. Not a performance activity.**

The burning-surface reduction is inherently global. If implemented rank-locally, multi-rank chamber results are already wrong — each rank integrates its own pressure from a partial surface. Physics error, not numerics.

Questions:

1. Is the reduction global or rank-local?
2. Does the port survive domain decomposition — ghost fill and BC at rank boundaries?
3. Is GPU-aware MPI actually active, or silently staging halos through host?

Question 3 matters now, not at Phase 4: silent host staging violates the governing rule, pollutes the baseline, and makes multi-rank T3 unpassable.

**Note:** the claim that local hardware cannot run device arena is contradicted by this probe — with a smaller `init_size` it omits the managed override and exercises the AMReX device default successfully. Local HMM still limits some defect detection, but device-arena testing is not categorically unavailable locally.

**Exit:** two-rank run reproduces single-rank pressure history to solver tolerance. Failure is fixed before Phase 1a, not deferred.

---

## 14. Instrumentation and figure set

### 14.1 Figure set

| # | Figure | Reads on | Status |
|---|---|---|---|
| F1 | Annotated timeline, one full timestep | T6a, T6b | Missing |
| F2 | Memcpy summary: bytes, **count**, **blocking calls**, by direction | T3a-c | Missing; v1 lacked count and blocking columns |
| F3 | Unified-memory page-fault count | T2 | Missing |
| F4 | Kernel time Pareto, **discovered** top 10 | Targeting | Present but hardcoded to three ranges; must be discovered |
| F5 | Roofline, top kernels, **annotated with actual limiter** | T7 | 3D missing |
| F6 | Achieved bandwidth vs applicable roof, per kernel class | T7 | Pending T7 redefinition |
| F7 | Occupancy with limiting resource identified | T7 | |
| F8 | **Arena live allocations and high-water** over steps | T4, T5a, T5b | Missing; `--cuda-memory-usage` sees CUDA backing allocations, not live arena allocations — needs arena-level instrumentation |
| F9 | GPU idle fraction inside step loop | T6b | Not computed by the script |
| F10 | Wall time per step, phase over phase | **T0** | Promoted from figure to gate |

F1 is read first at every phase boundary. Human pattern recognition on a timeline outperforms aggregate metrics at "why is there a gap here."

### 14.2 Ledger

One CSV, appended at every phase boundary, keyed `(phase, commit_sha, tree_hash, case, n_ranks, device, build_config)`. Columns: T0 through T7 plus wall time.

Purpose is regression detection across phases. Reuse the existing path-keyed ledger machinery.

### 14.3 Automation

Target: one command emits the full figure set plus the ledger row.

```
make profile CASE=<deck> PHASE=<n> BUILD=<coarse|fine>
```

This is not QOL. It gates whether Phases 1-3 can be agent-driven at all: if T3, T4, and T6 emerge as ledger numbers, agents self-certify and humans review at boundaries. If not, humans are in the loop on every check.

Build it during Phase 0.

---

## 15. Execution model

### 15.1 Phase tractability

| Phase | Agent fit | Human load | Note |
|---|---|---|---|
| 0.5 | Low | Low | Run, look, verdict |
| 0 | High to build, low to interpret | Medium | Automation is ideal agent work; reading F1 is not |
| 1a | Highest in plan | Medium volume, low stakes | Mechanical, pattern-referenced, hard oracle |
| 1b | Medium | Medium | Lifetime changes are correctness-sensitive; regrid interaction needs human review |
| 2 | **Lowest in plan** | Low volume, highest stakes | Coupled design; not decomposable |
| 3 | High to sweep, low to restructure | Low | Parameter search is agent work; fusion calls are not |

### 15.2 Review posture

Phase 1a is skim-for-pattern-violation on high volume. Phase 2 is line-by-line on a small diff. Different postures; do not conflate.

### 15.3 Predicted agent failures

| Failure | Phase | Countermeasure |
|---|---|---|
| Device-local reduce instead of Allreduce-shaped | 2 | Specify shape explicitly in the spec, not the goal |
| Over-porting during crash triage | 1a | Triage table as required classification step |
| Newton host branch judged correct | 2 | Explicit audit item, not general instruction |
| Fusion undoes Phase 1b | 3 | T4/T5 re-check as Phase 3 exit gate |
| Self-certification against inadequate oracle | All | Section 5.2; **already occurred** |

### 15.4 Human decisions

1. Chamber vs `multicomponent/FMA` at revised cost (Section 6)
2. Ballistic scope resolution (9.2)
3. Newton convergence policy (9.3)
4. T6b and T7 thresholds (5.6)
5. Multi-rank correctness verdict (Section 13)
6. Box size vs target rank count (Section 10)
7. Definition of "current state" for T1 (3.3)

Everything else is supervision.

### 15.5 Branch mechanics

- One commit per subsystem flip in Phase 1a, so bisect works when the harness reddens
- Tag at each phase exit; the ledger keys on `commit_sha` and tags make phase-over-phase comparison navigable
- Phase 2 on a child branch — the one redesign that may need abandoning
- **Decide now:** does `chamber-gpu` survive as fallback, or does `chamber-gpu-mem` become trunk? Determines whether Phase 1a edits must stay cherry-pickable
- Captures refuse to run on a dirty tree (5.1)

---

## 16. Pattern capture

Chamber work is also template generation for `multicomponent/FMA`. Capture as work proceeds:

| Pattern | Destination |
|---|---|
| Arena assignment by data class | `gpu_manual` |
| Async arena discipline for temporaries | `gpu_manual` |
| Allocation lifetime patterns (operator hoisting, persistent scratch) | `gpu_manual` |
| Device-resident scalar with Allreduce-shaped reduction | `gpu_manual` |
| Arena-level allocation instrumentation | `gpu_manual` |
| False-pass analysis method for target sets | `gpu_manual` — generalizes beyond this campaign |
| Mechanical sync inventory procedure | `gpu_manual` |
| Gap table, figure set, revised targets | `CHAMBER_GPU_MASTER_REPORT.md`, new section |

---

## 17. Risk register

| Risk | Phase | Response |
|---|---|---|
| **Metrics structurally blind to their pathology** | All | Section 3.1 meta-gate; validate against known in-tree instances |
| Campaign passes all targets, production slower | All | T0 blocking at every phase exit |
| Reduction is rank-local; existing multi-rank results invalid | 0.5 | Probe first; physics bug; fix before proceeding |
| Baseline not reproducible from dirty tree | 0 | 5.1; no dirty-tree measurement enters the ledger |
| Agent self-certifies against weak oracle | All | 5.2; already occurred once |
| Ballistic scope contradiction unresolved | 2 | Forced decision, 9.2; feeds Section 6 |
| Phase 1b lifetime changes break regrid | 1b | Explicit regrid coverage in harness |
| Idle profile dominated by removable debug sync | 0 | 5.3 measurement before scoping 2 and 3 |
| Fine-NVTX build perturbs T0/T6b | 0, 3 | 12.2 two-build policy |
| Footprint exceeds 40 GB under device arena | 0 | 5.5; FP32 storage becomes structural |
| GPU-aware MPI silently inactive | 0.5 | Verify at runtime; multi-rank T3 unpassable otherwise |
| Phase 3 undoes Phase 1b | 3 | T4/T5 re-check as exit gate |
| Chamber is not the throughput bottleneck | 0 | Section 6, reopened |

---

## 18. Open questions

1. Does chamber pressure currently update host-side, and where does the reduction land?
2. Is the burning-surface reduction global or rank-local?
3. Ballistic scope — which of the three options in 9.2?
4. What does "current state" mean for T1: AMReX defaults, `MODE=bench`, or production launch practice?
5. Available hardware — determines how much of Phase 4 is testable, and whether the 40 GB budget is the operative one
6. Does `chamber-gpu` survive as fallback branch?
7. Does the harness cover pressure divergence, regrid, two ranks, and restart, or does 5.2 require building it?

---

## Note on specifics

Configuration parameter names, AMReX NVTX integration, arena API surface, and Nsight flags shift between versions. Treat every such name as a pointer to verify against the vendored AMReX and installed toolkit, not as a literal.

**Confidence:** ~90% on strategy, phase ordering, and the v2 target set. ~70% on parameter and flag spellings. The v1 target set was authored at similar stated confidence and was wrong in four places for a single structural reason — Section 3.1 exists so that confidence is earned by adversarial test rather than asserted.
