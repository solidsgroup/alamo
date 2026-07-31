# Phase 0 (v2) — deliverable status

Campaign PLAN v2 §5.7 lists nine deliverables. v1's Phase 0 record lives in
`docs/agent_plans/20260727-phase0-baseline/results/RESULT.md` and is not
superseded — this file tracks the v2 set, which is a different and larger set.

| # | Deliverable | § | State |
|---|---|---|---|
| 1 | Reproducible baseline capture | 5.1 | **DONE** — §A below |
| 2 | Fixed oracle with stated coverage | 5.2 | **BLOCKED, root cause diagnosed** — §H |
| 3 | `AbortIfDeviceError` disabled comparison | 5.3 | **HUMAN CHECKPOINT** — §H |
| 4 | Mechanical sync inventory | 5.4 | **ENUMERATED; aggregate cost measured; per-site ranking partial** — §H |
| 5 | Footprint budget including 40 GB | 5.5 | **DONE, measured** — §H |
| 6 | T0/T6b/T7 thresholds | 5.6 | **T0 DONE; T7 classified; T6b and numeric T7 pending human/instrumentation** — §H |
| 7 | Target set v2 false-pass validated in-tree | 3.1/3.3 | **PARTIAL** — §F |
| 8 | Gap table | — | **DONE** — §H |
| 9 | Revised cost estimate | 6 | **DONE; P4 decision pending** — §H |

---

## §A Capture reproducibility (P2) — DONE

`benchmark/phase0_capture.sh push` now ships three provenance artifacts instead
of a bare dirty-file count:

| Artifact | Contents |
|---|---|
| `_pushed_rev.txt` | host, time, branch, `local_head`, `local_dirty_files`, **`tree_hash`**, manifest size |
| `_pushed_tree.diff` | full `git diff HEAD` of tracked modifications |
| `_pushed_manifest.tsv` | sha256 of every file actually shipped |

`tree_hash = sha256(HEAD + sha256(tree.diff) + sha256(manifest))`, first 16 hex.
Content-derived, so the same source yields the same key committed or not.

Ledger key is now `(local_head, tree_hash, case, n_ranks, device, build_config)`
per §14.2. The env leg already `cat`s `_pushed_rev.txt`, so `tree_hash` reaches
every capture with no further change.

**Policy is RECORD, not REFUSE.** §5.1 offers both; refusing would make the
harness unusable exactly when it matters, since the working tree carries the
changes under test. `STRICT=1` restores refuse-on-dirty for release captures.

First run: `tree_hash=cda69b0d59e5d50d`, `local_head=ab3ea4c45`,
`local_dirty_files=573`, `manifest_files=3034`.

## §B Oracle integrity (P3) — PARTIAL

Done: `status.sh` split fast/`FULL=1`, honest leg names, `BLOCKED` distinct from
`FAIL`. Detail in campaign `NOTES.md` under "v2 precondition verification".

Not done: the four §5.2 coverage requirements — pressure history over enough
steps to diverge, at least one regrid, two ranks, at least one
checkpoint/restart cycle. None exist in either mode. This is harness
construction, not a tweak, and wants its own folder.

**The oracle is currently RED** (N11). §5.2 is not satisfiable while it is.

## §C `AbortIfDeviceError` early measurement — NOT DONE

Sites confirmed (see §D), measurement not run. It needs an A/B against a build
with the facility disabled, and there is no runtime switch — only a compile-time
edit, which Phase 0 is not authorized to make. Two ways forward:

1. Scratch patch, measure, revert, commit nothing but the numbers.
2. Fold it into the N10 fix (apply Elastic's existing `#ifdef AMREX_DEBUG`
   guard to Flame), which makes the A/B a by-product.

Option 2 is better: it delivers the measurement and the fix in one tier-2 task
instead of doing throwaway work.

## §D Mechanical sync inventory — enumerated, not yet ranked

§5.4 requires enumeration from source, per-step frequency, and a ranking by
**measured** cost. The first two are below. Ranking is blocked on the nsys
re-capture, and is deliberately not guessed.

### Explicit device-wide syncs

| Site | Frequency | Note |
|---|---|---|
| `Util/Util.H:156` | per `AbortIfDeviceError` call | `streamSynchronizeAll()`, then `flag.value()` which itself syncs + D2H |
| `Integrator/Base/Mechanics.H:239` | 1 per elastic solve | `elastic_op` lifetime UAF guard; comment at `:229-238` explains |
| `Solver/Nonlocal/Newton.H:860` | per NR iteration | `dsol_mf` overwrite UAF guard |
| `Solver/Nonlocal/Newton.H:1001` | per NR iteration | after `MultiFab::Add` into the solution |
| `Solver/Nonlocal/Newton.H:777` | per line-search **backtrack** | inside `restore_baseline` lambda |
| `Solver/Nonlocal/Newton.H:1294` | per fill call | after `FillBoundaryAndSync` loop |
| `Solver/Nonlocal/Newton.H:1300` | per `HasNonFinite` call | first statement of the function |
| `Util/BMP.H:62`, `IC/StarAftGrain.H:382` | IC/IO only | not in the step loop |

### Amplifiers — what actually multiplies the above

`AbortIfDeviceError` callers: `Flame.cpp:861` and `Flame.cpp:915`, both
**unguarded**, per level per `Advance` → the "two full-stream syncs per Flame
level per step" §5.3 predicts. `Elastic.cpp:120,391,517` and `Newton.H:439,569`
are `#ifdef AMREX_DEBUG`-guarded and compile out of release.

`AbortIfNonFinite` → `HasNonFinite` → sync at `Newton.H:1300`. **Ten call
sites**: `:728, 750, 771, 783, 790, 849, 854, 961, 972, 1002`.

Per Newton iteration, non-line-search path (`:956-1002`): `:961`, `:972`,
`:1001`, `:1002` = **4 device-wide syncs minimum**.

Line-search path (`:723-909`): `:728`, `:750`, `:771`, `:849`, `:854`, `:860`
= 6, **plus** `:777`, `:783`, `:790` per backtrack, **plus** `FieldNorm0` at
`:741`, `:822`, `:855` (each a device reduction + sync + MPI allreduce).

Newton iterations run to `m_nriters` per elastic solve; elastic solves fire on a
40-50 step interval. So the Newton syncs are bursty, not per-step — which is
exactly why a per-step average would misrank them and why §5.4 demands measured
cost rather than a count.

### Host-landing reductions in the step loop

| Site | Frequency |
|---|---|
| `Flame.cpp:1103` | per box, per level, per step (NOTES N1) |
| `Flame.cpp:663` + `:671-675` | per `thermo.interval`; 1 reduce + 5 `ReduceReal` back to back (NOTES N2) |
| `Base/Mechanics.H:564,587,610,629` | 4 per `Integrate` |
| `Newton.H:1333/1336` | per `HasNonFinite` (device reduce + `ParallelAllReduce::Or`) |
| `Newton.H:1429-1433` | per convergence-metric evaluation |

### Not a contributor — checked and cleared

`Operator.cpp:210` and `:222-231` are 11 host-landing reductions, but they sit
inside `AlamoPrintSpatialDiagMF`, called only under `spatial_diag_enabled` /
`diag_enabled` / `resid_diag_enabled` (`Operator.cpp:660,665,1095,1099`). Debug
diagnostics, off by default. Flagging them would have been the authored-inventory
error §5.4 warns about, in the opposite direction.

### Coverage gap found while doing this

`Newton.H` contains **two** `BL_PROFILE` markers, both named
`Solver::Nonlocal::Newton::DW()` (`:1127`, `:1171`). The entire Newton solve —
including every sync above — is unattributed in the NVTX timeline. F1 cannot
show where Newton time goes. Adding regions is a tier-2 edit; logged, not made.

## §E Footprint budget — 40 GB clears on arithmetic; measurement pending

Local anchor: 350-380 B per resident cell (2D, elastic, FP64), ~267 MiB fixed
CUDA context.

Deck geometry, from the decks themselves:

| deck | dim | base `n_cell` | `max_level` | base cells |
|---|---|---|---|---:|
| `input` | 2D | 64 64 (64) | 3 | 4,096 |
| `input_copy` | 2D | 128 128 | 1 | 16,384 |
| `input_3d_centre_bore_128_a2` | 3D | 128 128 64 | 1 | 1,048,576 |

Worst case is every level fully covered at `ref_ratio` 2, i.e. ×4 cells per level
in 2D and ×8 in 3D:

| deck | fully-refined cells | × 380 B | + context |
|---|---:|---:|---:|
| `input` | 348,160 | 132 MB | ~400 MB |
| `input_copy` | 81,920 | 31 MB | ~300 MB |
| `input_3d_centre_bore_128_a2` | 9,437,184 | 3.59 GB | ~3.9 GB |

**The 40 GB budget clears by an order of magnitude on every current deck**, and
would still clear at 3× the per-cell cost (3D worst case → ~10.8 GB).

Two honest caveats:

1. 350-380 B/cell was measured in **2D**. 3D carries a third displacement
   component, a larger `Matrix4` model, and bigger ghost volumes, so the true 3D
   per-cell figure is higher — unmeasured, and the ×3 sensitivity above is a
   guard, not a measurement.
2. This is the *model*, not a high-water. §5.5 prefers a measured arena
   high-water where one exists. It does not yet; the re-capture supplies it.

**FP32 storage is not structural** — the answer v1 gave, now supported at 40 GB
rather than only at 8/80 GB. The binding constraint is future decks, not these.

Note the earlier flip leg ran `the_arena_init_size` at 0.5 × 80 GB = 40 GiB, but
that is an *initial* size on an 80 GB card, not a cap, so it did not test the
budget. AMReX exposes no hard arena ceiling; the high-water measurement is the
real test.

## §F False-pass validation — PARTIAL

§3.1 admits a target only once its metric is demonstrated to catch a known
in-tree instance. T4's three named instances:

| Instance | Status |
|---|---|
| MLMG operator + solver constructed per elastic solve | **CONFIRMED** — `Base/Mechanics.H:211` constructs `Operator::Elastic<MODEL::sym> elastic_op(...)` inside `TimeStepBegin`, destroyed at scope exit |
| Per-step device error flag: `DeviceScalar` alloc + H2D + D2H | **CONFIRMED** — `Util/Util.H:77` wraps `amrex::Gpu::DeviceScalar<int>`; constructed unguarded at `Flame.cpp:747` |
| `FieldNorm0` allocating composite MultiFabs per call | **CONFIRMED** — `Newton.H:1351-1358` allocates a full `MultiFab` **per level** (`make_unique`, full `nComp` and `nGrowVect`) and `MultiFab::Copy`s every level into it, then refluxes (`:1362-1365`). Called at `:741, :822, :855, :1382`, including inside the line-search backtrack loop |

All three named instances exist. **The metric still does not.** Catching them
requires arena-level request counting, not `cudaMalloc` counting — which is the
whole point of the T4 redefinition, since AMReX arenas cache and reuse backing
allocations and would report zero `cudaMalloc` for all three. Until that
instrumentation exists, T4 is defined but unvalidated and **P1 does not clear**.

`FieldNorm0` is the worst of the three by inspection: it copies the entire
solution field, at every level, per norm evaluation, and the line-search path
evaluates it per backtrack. That is a per-call allocation proportional to the
whole field, not a scalar. It also carries the `Reflux` call, so it is not a
pure norm — removing the copies is a correctness-sensitive change, not a
cleanup.

---

## §G Re-capture — submitted 2026-07-28

Build `11774784` COMPLETED (11:33), all four sm_80 binaries rebuilt from the
current tree.

| Job | Deck | Dim |
|---|---|---|
| 11775088 | `input_copy` | 2D |
| 11775089 | `input` | 2D |
| 11775090 | `input_3d_centre_bore_128_a2` | 3D |

Provenance: `local_head=ee65247b7`, `tree_hash=5e0efcf7a5766c1b`,
`src_hash=73126499fd418d70`, `local_dirty_files=575`, 3,036 files manifested.

This run carries four harness changes that the previous one did not:

1. nsys no longer traces MPI — the previous three captures produced **no**
   `trace.nsys-rep` at all, so F1/F2/F3/F8/F9 were empty.
2. ncu runs `--kill yes`, so it stops the app after its launches instead of
   riding the full horizon and being OOM-killed.
3. F10 has a `max_step=1` startup calibration and 3 reps with a standard
   deviation, so it is admissible for the first time.
4. The N4 probe reads the MPI the binary actually links, not the module.

All three COMPLETED (8:53 / 8:04 / 6:44). CSVs and logs collected under
`results/figures/` (47 MB; the 173 MB `.ncu-rep` binaries stay on NOVA).

### Outcome per leg

| Leg | Result |
|---|---|
| env | ok, all three |
| timing | **ok — admissible for the first time** |
| flip | ok, `failures.txt` empty on all three |
| ncu | **FIXED** — 3/3 ranges on every deck, `mgVcycle` no longer OOM-killed |
| nsys | **STILL DEAD on all three** |

### F10 — startup-corrected, 3 reps, with uncertainty

| deck | arm | startup | wall mean | steady per step |
|---|---|---:|---:|---:|
| `input_copy` | managed | 3.921 s | 35.703 ± 0.984 | 0.35313 ± 0.01093 |
| `input_copy` | device | 2.462 s | 35.766 ± 0.539 | 0.37005 ± 0.00599 |
| `input` | managed | 3.196 s | 16.953 ± 0.849 | 0.12506 ± 0.00771 |
| `input` | device | 2.798 s | 16.382 ± 2.991 | 0.12349 ± 0.02719 |
| 3D | managed | 4.186 s | 15.814 ± 4.218 | 0.10571 ± 0.03835 |
| 3D | device | 2.600 s | 13.176 ± 0.750 | 0.09614 ± 0.00682 |

**No arena verdict.** Every per-step difference is inside 2 sd of run-to-run
scatter, and the scatter is large and asymmetric — one arm's sd is 3-5× the
other's on two of three decks, meaning individual reps were disturbed. n=3 is
not enough. Raise `REPS` and prefer a quiet node before anyone claims a
direction.

**Startup is the one clean signal:** the device arena initializes 0.4-1.6 s
faster than managed on every deck. On `input_copy` that offsets a slower steady
state so exactly that the raw wall means are 35.703 vs 35.766 — the
uncalibrated comparison would have reported "no difference" and hidden both
effects. This is why §12.1 requires the calibration.

Every earlier F10 number in this campaign, including "device −5.5% on
`input_copy`", was a single unrepeated run and is withdrawn.

### F5/F6/F7 — T7 is answered, and v1's assumption was wrong

3D, `Operator::Elastic::Fapply()`, 135 launches:

| metric | value |
|---|---|
| Registers per thread | **254** |
| Block Limit Registers | **2** blocks/SM |
| Theoretical occupancy | 12.50% |
| Achieved occupancy | 11.95% (96% of theoretical) |
| Compute (SM) throughput | 31.0% |
| **DRAM throughput** | **11.9%** |
| Duration | ~552 µs |

**These kernels are nowhere near the HBM roof.** DRAM throughput is 11.9% and
compute 31%; neither is the limiter. The limiter is occupancy, and occupancy is
capped by *registers* — 254 regs/thread gives a hard 2-block/SM ceiling and a
12.5% theoretical occupancy that the kernel already achieves 96% of.

§3.3's T7 rewrite predicted exactly this: "a kernel can sit far from the HBM
roof and be correctly optimized." Confirmed by measurement, not argument.
**T7's threshold for this class is a register/occupancy target, not a bandwidth
target.** A bandwidth-roof target would have declared Fapply an 88%-headroom
opportunity and sent Phase 3 chasing memory traffic that is not the constraint.

Setting a number still needs the missing F1: with no timeline, the fraction of
wall these kernels actually own is unmeasured here.

### F4 is confirmed broken as §14.1 says

`OperatorElasticFapply.csv` and `OperatorFsmooth.csv` contain the **same
kernel** — identical name, identical grid `(2315,1,1)`, identical 135 launches.
`Operator::Elastic::Fapply()` is nested inside `Operator::Fsmooth()`, and
`--nvtx-include "range/"` takes nested launches. The three hardcoded ranges
resolve to one distinct elastic kernel plus MLMG's extra `IsFabArray` kernel.

So the ncu numbers are sound but the labels are not: this is one kernel measured
three times, not three kernels. §14.1's F4 requirement — "**discovered** top 10,
must not be hardcoded to three ranges" — is now demonstrated rather than
asserted. Until F4 discovers kernels by time, the Pareto does not exist and
"top kernels" is an assumption.

### nsys — H1 fix failed, root cause restated

Dropping `mpi` from `-t` did not help. The injection library loads regardless of
the trace list. The real abort:

```
terminate called after throwing an instance of 'boost::wrapexcept<std::runtime_error>'
  what():  Expected shared object name, found a path delimiter
```

thrown from `libToolsInjection64.so`, during `ompi_mpi_init`, then SIGKILL. My
first reading blamed the MPI trace because MPI frames dominated the backtrace;
that was where it aborted, not why.

Two guard defects this exposed, both now fixed (`d9e2005be`):

- `CAPTURE_FAILED` tested `-f` on `trace.nsys-rep`. nsys creates that file up
  front, so a killed capture leaves a 0-byte file that passes. On the 3D deck
  the guard stayed silent and `ls` showed a full set of reports, all 0 bytes.
- The stats loop trusted the exit code. `nsys stats` exits 0 after writing a
  0-byte CSV when the sqlite export fails.

Next step is isolation, not another guess: `nsys profile -t cuda,nvtx` on a
trivial binary, with and without `srun`. Every run here is `-n1`, so dropping
`srun` for the nsys leg is the leading candidate.

Collect with `bash benchmark/phase0_capture.sh collect <remote_dir>`.

## Blocking summary

- **P1:** request-layer instrumentation detects churn, but live allocation
  count and direct `FieldNorm0` attribution are still missing.
- **P3:** the two-rank restart oracle is red because required thermal history is
  not checkpointed. The source repair is at the required human checkpoint.
- **P4:** the revised decision gate below has not been answered.
- **Local GPU availability:** the final rerun is blocked by CUDA error 803.
  NOVA evidence is complete, but the local FULL gate cannot currently run.

---

## §H Authoritative 2026-07-31 addendum

This section supersedes the earlier “submitted/pending” capture status. The
full evidence and invalid-retry history are in
`phase0-baseline/results/figures/PHASE0_CAPTURE_SUMMARY.md`; the concise
measurement narrative is in `phase0-baseline/results/RESULT.md` H17-H21.

### H.1 Configuration and provenance

No metric below combines incompatible configurations:

| Purpose | Runtime/build | Provenance |
|---|---|---|
| T0 | plain binary, paired short/long, balanced five-run order | `2057d3206 / 05fe0312f5798018` |
| F1-F3/F9 | fine-NVTX profile binary, managed arena, 90/110-step trace | `ce3e47acf / 073ff4f513fbcf63` |
| T4/T5a | profile binary, device arena, endpoint request tables | `ce3e47acf / 073ff4f513fbcf63` |
| T5a production stability | profile binary, device arena, 6,000 steps | `1433d55a0 / 293f2d1b0e208804` |
| NCU selector recovery | fine-NVTX profile binary, managed arena | supplemental captures, selector table retained per capture |

All rows run binaries built from `src_hash=c88836ce414b44cc`.

### H.2 T0 baseline

| Case | Device median (s/step) | MAD (s/step) | Managed median | Paired arena verdict |
|---|---:|---:|---:|---|
| `input_copy` | **0.42729** | 0.00156 | 0.42750 | inconclusive |
| `input` | **0.08290** | 0.00202 | 0.08512 | inconclusive |
| 3D centre-bore | **0.09271** | 0.00076 | 0.09562 | inconclusive |

These device medians/MADs are the Phase 0 T0 baseline and uncertainty band.
None licenses an arena-speedup claim.

### H.3 Target gap table

| Target | Phase 0 observation | State | Closure |
|---|---|---|---|
| T0 | Balanced five-run medians above; all paired 2-MAD bands overlap zero | **BASELINED** | Compare every phase against device median + uncertainty |
| T1 | Device-arena endpoints have zero application managed requests; production source has no managed call, while `src/Test/BC/Constant.H` does and diagnostic launchers intentionally enable managed | **HUMAN DEFINITION NEEDED** | Adopt production launch practice (`the_arena_is_managed=0`) as “current state,” recommended |
| T2 | UM reports empty/unavailable, not zero | **NA/BLOCKED ON T1** | Re-capture only after T1 configuration is fixed |
| T3a-c | Single-rank transfer bytes/counts measured; blocking memcpy count is zero; linked MPI has CUDA support disabled | **GAP** | Attribute small transfers; GPU-aware MPI build required before multi-rank pass |
| T4 | Primary 6,000-step run: 19,829.974 device requests/step and 15,615.913 pinned requests/step; known flag/operator churn detected | **FAIL / P1 PARTIAL** | Directly attribute `FieldNorm0`, then justify or remove every steady request class |
| T5a | Primary 135→138 MiB over 6,000 steps; 3D peak 26,569 MiB | **PASS** | Preserve through later phases; 40 GiB clears |
| T5b | AMReX table emitted after teardown, no live endpoint count | **UNAVAILABLE / P1 BLOCKER** | Human-approved in-evolution memory endpoint probe |
| T5c | `Ballistic::Advance` appends to unread `dpdt` vector every step | **FAIL** | Bound/remove history or expand Ballistic scope |
| T6a | 7 explicit stream sync sites, 10 device-result landings, 11 blocking collectives; profile traces show 1,559-33,656 CUDA sync calls/step | **FAIL / PARTIALLY RANKED** | Abort A/B first; use coarse ranges to attribute remaining gaps |
| T6b | Fine diagnostic idle is 65.9%, 60.9%, and 10.2% for the three decks | **NO GATE VALUE** | Fine NVTX is inadmissible; approve coarse instrumentation and threshold review |
| T7 | 2D is launch/underfill-limited; 3D `Fapply` and `prepareForSolve` are register-limited; `SetModel`/`Fsmooth` are memory-dominant | **CLASSIFIED, NUMERIC GATE PENDING** | Human signs off per-class targets |

### H.4 False-pass meta-gate

Closed:

- full payload provenance and source hashes;
- failed required SLURM legs propagate;
- missing scheduler/profile evidence differs from zero;
- semantic NCU owner/arity validation catches the known wrong kernels;
- legacy and unbalanced timing cannot enter T0;
- Tier 2 requires both sanitizer cleanliness and application exit zero;
- the top-level status script propagates FAIL and BLOCKED exit codes.

Open:

- T5b still cannot observe its claimed property;
- `FieldNorm0` is not directly attributable at the request layer;
- T6b has no admissible coarse capture;
- the two-rank/regrid/restart oracle is red.

P1 and P3 therefore remain blocking.

### H.5 Revised cost estimate for P4

These are engineering ranges, not elapsed cluster time:

| Scope | Estimated effort | Main risk |
|---|---:|---|
| Finish Phase 0 (oracle, live endpoint, coarse NVTX, Abort A/B, reruns) | 2-4 days | checkpoint compatibility and human review |
| Phase 1a arena hygiene | 1-3 days | launch/configuration coverage |
| Phase 1b lifetime redesign | 6-10 days | operator lifetime across regrid; persistent `FieldNorm0` scratch |
| Phase 2 without Ballistic port | 4-7 days | cannot fully clear T3/T5c |
| Phase 2 with Ballistic/global-reduction path | 7-12 days | device-result reduction, Allreduce shape, diagnostic contract |
| Phase 3 tuning | 4-8 days | register spilling, workload-dependent underfill |

Full target compliance including Ballistic is therefore roughly **20-33
engineer-days**, plus queue time and required human review. Omitting the
Ballistic port is cheaper but knowingly leaves the governing residency target
and T5c unmet; it is not an equivalent completion path.

### H.6 Required human decisions

Immediate correctness checkpoint:

1. Approve persisting `temps` and `thermal.has_exceeded_Tcutoff`, plus an
   explicit abort when an old checkpoint lacks a required current field.

Before Phase 1:

2. Decide chamber versus `multicomponent/FMA` at the 20-33 day estimate and
   state whether chamber currently gates SRM paper throughput.
3. Resolve Ballistic scope: port it, or explicitly accept that T3/T5c cannot
   pass. Duplicating its formula in Flame is not recommended.
4. Define T1 current state. Recommendation: production launch practice with
   `amrex.the_arena_is_managed=0`; managed diagnostic launchers are not the
   production state.
5. Approve the scratch-only Abort A/B, the in-evolution T5b endpoint, and a
   coarse-NVTX build/capture.
6. Sign off T6b/T7 numeric thresholds after the admissible coarse capture.

Until those checkpoints clear, **Phase 1 is not authorized**.
