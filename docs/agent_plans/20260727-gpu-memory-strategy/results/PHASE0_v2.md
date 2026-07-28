# Phase 0 (v2) — deliverable status

Campaign PLAN v2 §5.7 lists nine deliverables. v1's Phase 0 record lives in
`docs/agent_plans/20260727-phase0-baseline/results/RESULT.md` and is not
superseded — this file tracks the v2 set, which is a different and larger set.

| # | Deliverable | § | State |
|---|---|---|---|
| 1 | Reproducible baseline capture | 5.1 | **DONE** — §A below |
| 2 | Fixed oracle with stated coverage | 5.2 | **PARTIAL** — §B |
| 3 | `AbortIfDeviceError` disabled comparison | 5.3 | NOT DONE — §C |
| 4 | Mechanical sync inventory | 5.4 | **ENUMERATED, unranked** — §D |
| 5 | Footprint budget including 40 GB | 5.5 | **DONE (model); high-water pending** — §E |
| 6 | T0/T6b/T7 thresholds | 5.6 | BLOCKED — re-capture IN FLIGHT (11775088/89/90) |
| 7 | Target set v2 false-pass validated in-tree | 3.1/3.3 | **PARTIAL** — §F |
| 8 | Gap table | — | BLOCKED on 1, 6 |
| 9 | Revised cost estimate | 6 | BLOCKED — N11 moves it again |

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

Collect with `bash benchmark/phase0_capture.sh collect <remote_dir>`.

## Blocking summary

- **N11** — live 3D OOB in `Elastic::Diagonal`, oracle RED. Blocks §B, and
  every cost estimate downstream.
- **NOVA re-capture** — build `11774784` running; captures not yet submitted.
  Blocks §5.6 thresholds, the gap table, and F1-F9.
- **P1** does not clear without arena-level instrumentation (§F).
- **P4** is a human decision and is now against a cost that N11 has moved.
