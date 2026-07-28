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
| 5 | Footprint budget including 40 GB | 5.5 | PARTIAL — §E |
| 6 | T0/T6b/T7 thresholds | 5.6 | BLOCKED on re-capture |
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

## §E Footprint budget — PARTIAL

Local anchor stands: 350-380 B per resident cell (2D, elastic, FP64), ~267 MiB
fixed CUDA context. The 8 GB and 80 GB cases are covered.

**The 40 GB case is still untested** and is the one §5.5 names, because the
hardware in hand is an A100-SXM4-**80** GB. It can be forced with
`the_arena_init_size` at half of device memory, which the capture already
parameterizes via `ARENA_FRAC`. Not yet run.

FP32-is-structural remains answered NO on the 8/80 GB evidence; that answer is
provisional until the 40 GB case runs.

## §F False-pass validation — PARTIAL

§3.1 admits a target only once its metric is demonstrated to catch a known
in-tree instance. T4's three named instances:

| Instance | Status |
|---|---|
| MLMG operator + solver constructed per elastic solve | **CONFIRMED** — `Base/Mechanics.H:211` constructs `Operator::Elastic<MODEL::sym> elastic_op(...)` inside `TimeStepBegin`, destroyed at scope exit |
| Per-step device error flag: `DeviceScalar` alloc + H2D + D2H | **CONFIRMED** — `Util/Util.H:77` wraps `amrex::Gpu::DeviceScalar<int>`; constructed unguarded at `Flame.cpp:747` |
| `FieldNorm0` allocating composite MultiFabs per call | **NOT VERIFIED** — `Newton.H:1349`; called at `:741, :822, :855, :1382` |

The instrumentation that must catch these does not exist yet: arena-level
request counting, not `cudaMalloc` counting. Until it exists T4 is defined but
not validated, and P1 does not clear.

---

## Blocking summary

- **N11** — live 3D OOB in `Elastic::Diagonal`, oracle RED. Blocks §B, and
  every cost estimate downstream.
- **NOVA re-capture** — build `11774784` running; captures not yet submitted.
  Blocks §5.6 thresholds, the gap table, and F1-F9.
- **P1** does not clear without arena-level instrumentation (§F).
- **P4** is a human decision and is now against a cost that N11 has moved.
