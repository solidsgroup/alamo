# READ_FIRST_NEXT_STEP — Agent Entry Point

> **This file is for any agent or collaborator starting cold on the chamber-gpu
> work.** Read this before touching any code. It tells you where we are, where to
> look, and exactly what process to follow before taking any action.

---

## 1. Start here: the almanac index

All documentation for the GPU port lives in `benchmark/`. The canonical entry
points are:

| File | What it is |
|------|------------|
| `benchmark/README.md` | Build system, run commands, benchmarking harness |
| `benchmark/GPU_ROADMAP_V3.md` | **THE PLAN** — phases 1–5, each task (1.A, 1.B, …) with Goal/Steps/Done-when/Artifact |
| `benchmark/GPU_BRANCH_GUIDE.md` | Branch policy (never merge), build matrix, naming conventions |
| `benchmark/SUCCESS_BOOK.md` | Wins log — confirmed results, do not re-litigate |
| `benchmark/PERF_TRACKING.md` | Standing metric set, CSV schema, regression harness |
| `benchmark/ALPHA1_BASELINE.md` | Frozen alpha-1.0 reference (commit + inputs) |
| `benchmark/validate/` | **Physics validation suite** (v3 Phase 1) — one-command CPU/GPU compare with organized, diffable output bundles |
| `benchmark/archive/README.md` | History archive index — superseded v2 plan, Phase 0–5 + Phase A/C records. The measured evidence (`PHASE_A_FINDINGS.md`, `PHASE3_R3_crossover.md`, `G0_BASELINE_OF_RECORD.md`) lives here. |
| `benchmark/GPU_AUDIT_20260702.md` | Full-`src/` code audit (2026-07-02) — live `DeviceErrorFlag` defect, CI gap, Fapply perf anatomy, dormant landmines. Feeds `GPU_ROADMAP_V3.md` tasks 3.A/4.G/4.H. |
| `benchmark/GPU_STRUCTURAL_PLAN_20260703.md` | Structural-change plan (2026-07-03) — elastic kernel restructure options (incl. stored-stencil proposal 3.G), phase-field sync/launch findings (PF.1–PF.5), multi-GPU shortcoming re-diagnosis, local 3D register A/B for the C1 edits. Feeds v3 Phases 3/5. |

The auto-memory index at
`/home/jackplum/.claude/projects/-home-jackplum-Projects-alamo/memory/MEMORY.md`
has one-line pointers into each memory file and is the quickest way to get
current project state without reading every doc.

---

## 2. Read and understand first

Before proposing any next step:

1. **Read `benchmark/GPU_ROADMAP_V3.md` in full.** It is the live plan (phases
   1–5; tasks 1.A…5.B, each with explicit STATUS). Know which phase is active and
   why; §2 maps what v2 settled vs what v3 inherits. v3's governing rule: **no
   optimization ships without a physics-error budget pass** (Phase 1).
2. **Read `benchmark/archive/PHASE_A_FINDINGS.md`** to see what has actually been
   measured on A100. Many numbers that "seem like guesses" are now confirmed.
3. **Check `benchmark/SUCCESS_BOOK.md`** for wins that are settled fact — do
   not redo analysis the book already records.
4. **Scan the memory index** for anything that might be stale (memory can lag
   reality — cross-check against current code/files when in doubt).

---

## 3. Identify the next task

After reading, find the next actionable item:

- Open `benchmark/GPU_ROADMAP_V3.md` and scan each task's **STATUS** line.
- The active phase is shown in the roadmap header. Find the lowest-numbered task
  in the active phase that is **not** marked DONE.
- Confirm its **Depends** field — ensure all prerequisites are actually complete
  (do not trust "DONE" status alone; spot-check the named artifact file exists).
- If multiple tasks are unblocked, pick the one with the lowest ID **or** the one
  the user most recently discussed.

---

## 4. Flesh out the task into an actionable plan

Before touching any code, write out the following (in the conversation — not in
a file unless the user asks):

**What:** One sentence describing the concrete output.

**How:** Step-by-step actions. For each step name:
- The exact file(s) to change or commands to run.
- Which build configuration to use (strict/nofast vs fast-math; 2D vs 3D; A100
  vs local A1000).
- Which input config to pass (`input_3d_centre_bore_256_a2_tuned`, etc.).

**How to test:** The exact correctness gate to run after each change:
- CPU golden compare (`benchmark/golden_compare_flame.sh`) for source edits.
- GPU test suite (`tests/GPU/run_gpu_tests.py`) for anything touching device paths.
- **Phase 1 physics-budget gate (`benchmark/validate/`) for any change touching
  the elastic hot path or an input lever — see §8a below. Mandatory before
  Phase 2/3 work lands; the v3 governing rule.**
- Standing metric set (see `benchmark/PERF_TRACKING.md`) for any performance claim.

**Definition of done:** Quote the **Done-when** clause from the roadmap task
verbatim, then confirm you will hit it.

---

## 5. Present the plan for user approval

Do **not** start executing until the user has explicitly approved the plan.
Present it as a short proposal (what/how/test/done-when), name any risks or
assumptions, and ask if this is the right task to tackle now.

---

## 6. Execute, test, and write results

Once approved:

1. Make the smallest change that satisfies the task. No side-quests; no
   speculative cleanup.
2. Run the stated correctness gate. Fix failures before declaring done.
3. Run the standing metric set if the task touches a hot path (Fapply, Newton,
   MLMG) — capture the before/after row and append to `benchmark/perf_regression.csv`.
4. Write results to the artifact file named in the roadmap task's **Artifact**
   field. If that file already exists, append a dated section rather than
   overwriting.
5. Update the task's **STATUS** line in `benchmark/GPU_ROADMAP_V3.md` to DONE
   (or PARTIAL with a note if blocked).
6. Update the relevant memory file under
   `/home/jackplum/.claude/projects/-home-jackplum-Projects-alamo/memory/` and
   its entry in `MEMORY.md` so the next agent starts with accurate state.

---

## 7. Worktree isolation and concurrent agents

**Always work in your own git worktree.** The main working directory
(`/home/jackplum/Projects/alamo`) may have other agents or the user actively
editing files. Do not work directly in it — check out an isolated worktree
before making any changes:

```bash
git worktree add /tmp/alamo-<task-id> chamber-gpu
cd /tmp/alamo-<task-id>
# ... do all work here ...
git worktree remove /tmp/alamo-<task-id>   # clean up when done
```

Assume other agents are running concurrently in other worktrees on the same
branch. Coordinate by writing to separate files and by reading the current
state of the branch before starting (run `git log --oneline -10` to check for
commits that landed since you began reading).

**Commits are fine without asking.** A clean, scoped commit with a descriptive
message is always welcome. **Pushes to the remote require explicit user
approval** — never run `git push` without asking first.

---

## 8. Where to run tests

**Prefer local testing.** The local workstation (A1000, sm_86) can run:
- CPU correctness gates (`golden_compare_flame.sh`, `tests/GPU/run_gpu_tests.py`
  in CPU mode)
- GPU correctness gates in 2D (the A1000's 8 GB is sufficient for 2D configs)
- Build-smoke and unit tests

**NOVA is a last resort.** Use NOVA only when:
1. The local machine physically cannot run the test (e.g. a 3D 256³ run that
   OOMs the 8 GB A1000), **and**
2. The test is on the critical path — not "nice to verify" but blocking the
   task's Done-when criterion.

If a NOVA test is genuinely necessary, **do not submit it yourself.** Instead:
- Write a self-contained test script (a `.slurm` batch file or a shell wrapper
  that calls `sbatch`) with all paths, flags, and expected output clearly
  documented.
- Place it under `benchmark/` with a name matching the task (e.g.
  `nova_c1_ab_test.slurm`).
- Tell the user exactly what to run, what output to look for, and where to
  record the result.

The user will submit and relay the output back to you.

---

## 8a. Phase 1 physics-budget gate (mandatory for hot-path edits)

Before any Phase 2/3 change (input lever or kernel rewrite) lands, run the
validation suite and check it against the recorded golden reference:

```bash
# 1. run the suite locally (strict/no-fast-math by default; add gpu_fast for a
#    perf-labelled smoke row only -- never use it for a correctness claim)
bash benchmark/validate/run_validation_local.sh --profiles cpu,gpu_strict

# 2. compare your candidate bundle against the committed reference, in GATE
#    mode (fails on ANY CORRECTNESS regression, or any ENGINEERING regression
#    beyond its physical-% budget -- not just a plain PASS/FAIL report)
python3 benchmark/validate/compare_validation.py \
  benchmark/validate/references/cpu_strict \
  benchmark/validate/runs/<your-new-gpu_strict-bundle> \
  --gate
```

Exit code 0 = clear to proceed; non-zero = read `compare_report.md` in the
candidate bundle directory for the named offending observable(s) before
continuing. This is the literal instrument the v3 governing rule (§0 of the
roadmap) refers to: *"No optimization ships without a physics-error budget
that proves the coupled burn/pressure/stress behavior is preserved."*

- Observable definitions + thresholds: `benchmark/validate/physics_budget.md`.
- Case set: `benchmark/validate/cases.manifest.yaml`.
- NOVA leg (only needed when a change is GPU-architecture-specific and the
  local A1000 can't exercise it, e.g. multi-box at 256³): submit
  `benchmark/validate/run_validation_nova.slurm` per
  `benchmark/NOVA_SLURM_RUNBOOK.md`'s conventions — same "do not submit it
  yourself" rule as §8 above.
- **Known gap:** NOVA's build scripts (`build_alamo_nova[_3d].sh`) only
  produce fast-math `--profile` binaries today; there is no NOVA strict/
  no-fast-math build target yet. `run_validation_nova.py --profile gpu_strict`
  fails loudly with the exact `configure`/`make` command needed rather than
  silently substituting the fast-math binary — extending those build scripts
  to add a `--cuda-fp strict` target is an open follow-up, not yet done.

---

## 8b. Known pending hardening item (small, unblocked)

`Util::DeviceErrorFlag::value()` (`src/Util/Util.H:68-78`) syncs only the
current CUDA stream, but `MFIter` cycles boxes across a stream pool — the same
class of bug as the historical `interpolation()` `tmpfab` elixir UAF. Device
error flags (`Newton::prepareForSolve`, `Fapply`, `Diagonal`) can be read
before an earlier-streamed box's `SetDeviceError` write lands, producing a
silent false negative on multi-box runs. Not yet fixed. Small, independent,
parallel-ok fix: add `amrex::Gpu::streamSynchronizeAll()` inside
`Util::AbortIfDeviceError` before the flag read. Full detail:
`benchmark/GPU_AUDIT_20260702.md` §2; tracked as roadmap task **4.G**.

---

## Key invariants — never violate these

- **`chamber-gpu` is never merged to master.** Phase done = Phase 5 definition-
  of-done, not a merge. See `benchmark/GPU_BRANCH_GUIDE.md`.
- **No kernel optimization without a counter that justifies it.** Every hot-path
  change needs an ncu or nsys number. See §1 of the roadmap.
- **Two builds, two purposes.** Performance headlines use the fast-math binary;
  correctness claims use strict/no-fast-math only. Never mix them in a comparison.
- **Combined physics (flame + elastic) is the unit of measurement.** Flame-only
  or elastic-only numbers are diagnostics, not the headline.
