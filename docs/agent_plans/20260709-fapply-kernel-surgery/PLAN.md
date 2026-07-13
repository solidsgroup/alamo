# TASK: fapply-kernel-surgery
# Folder: docs/agent_plans/20260709-fapply-kernel-surgery/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 3 (solver: Operator/Elastic hot kernels)                     |
| Model        | sonnet (user-directed for this session; orchestrator + fresh verifier compensate) |
| Verification | full-oracle (bit-exact intent: golden compare + device lint + sanitizer) |
| Est. scope   | 2 files src/ (Operator/Elastic.cpp, Set/Matrix4_Major.H), ~150 lines |
| Parallel-safe| no (touches elastic solver used by all open tasks)           |

## Operating rules

1. Read ONLY the files listed in Context budget. docs/archive/ is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: STOP,
   report discrepancy, wait.
3. One commit per step unless stated. Message: `<area>: <what> (<task-folder>)`.
   IMPORTANT: working tree is dirty with UNRELATED edits
   (src/Integrator/Flame.{cpp,H}, .githooks/pre-commit, benchmark manifest,
   many input_* files). Stage ONLY files this task touches.
4. No scope expansion. New ideas go to NOTES.md in this folder, not into code.
5. If required knowledge is missing, STOP and report.

## Context budget

Read first: this PLAN.md
Read: src/Operator/Elastic.cpp:317-900 (Fapply + Diagonal),
      src/Set/Matrix4_Major.H (whole file),
      benchmark/fapply_register_ab.sh (header comments for usage),
      benchmark/build_alamo_local_gpu.sh (header for env vars)
Reference only if step names it: benchmark/status.sh, benchmark/lint_device_patterns.sh
Forbidden: docs/archive/*, other task folders, src/Integrator/*

## Objective

PLAN.md (live plan) task 3.2b: cheap kernel surgery in the elastic hot kernels.
Fapply is ~74.6% of GPU kernel time. Four bit-exact-intent edits reduce
redundant global loads and instruction count. This session is LOCAL-ONLY
(A1000, sm_86): correctness gates are fully valid locally; wall-time numbers
are recorded with the shared/50W-cap caveat; A100 wall judgment deferred.

## Oracle

Command(s), all from repo root, all must pass:
- `benchmark/lint_device_patterns.sh` (exit 0)
- `GOLDEN_MODE=cpu benchmark/ci_golden_compare.sh` (exit 0)
- `TIERS=2 benchmark/local_a100_gate.sh` (exit 0; Tier 2 = compute-sanitizer
  memcheck — required by the hard rule for kernel edits)
- Standalone bit-exactness check (Step 1 builds it) comparing old vs new
  Matrix4 x Matrix3 path on randomized inputs: exact `==` equality.
Covers: CPU bit-exactness, device access patterns, device memory errors.
Does NOT cover: A100 occupancy/wall (NOVA off-limits this session);
local wall numbers are indicative only.

## Steps

### Step 0 - Baseline capture
VERIFY: `git status --porcelain -- src/Operator src/Set` shows no changes.
DO:
- Build local GPU 2D binary with TinyProfiler:
  `PROFILE=1 DIM=2 SMOKE=0 benchmark/build_alamo_local_gpu.sh`
- Build 3D object as well: `PROFILE=1 DIM=3 SMOKE=0 benchmark/build_alamo_local_gpu.sh`
- Save static register baseline:
  `bash benchmark/fapply_register_ab.sh <3d-binary-or-Elastic.o> | tee results/regs_baseline.txt`
- Timing baseline: pick an elastic-heavy short deck that exists in repo root
  (prefer `input_nova_centre_bore` if present; else smallest input_confirm_*),
  cap steps (max_step or short stop_time) so one run is 2-5 min on the A1000.
  Run 3x, save TinyProfiler tables (grep `Fapply`, `Diagonal`, `MLMG`) to
  `results/wall_baseline_runN.txt`. Record GPU clocks/power state
  (`nvidia-smi -q -d CLOCK,POWER | head -40`) alongside.
CHECK: baseline files exist in results/.

### Step 1 - Matrix4_Major x Matrix3 hand-unroll (edit c)
DO: In src/Set/Matrix4_Major.H:552-565, replace the branchy-accessor loop in
`operator*(const Matrix4<AMREX_SPACEDIM,Sym::Major>&, const Set::Matrix3&)`
with direct `data[]` indexing. HARD CONSTRAINT: the accumulation order of the
floating-point sum for each ret(i) MUST be identical to the existing loop
order (J outer, then k, then L) — FP addition is non-associative and the
golden gate will catch reordering. Derive the (i,J,k,L)->data[] map
mechanically from the class's own operator() accessor (2D: 10-entry map in
Matrix4<2,Major>; 3D: 45-entry map in Matrix4<3,Major>). Handle BOTH 2D and
3D (the operator is compiled per AMREX_SPACEDIM; use overloads or
`#if AMREX_SPACEDIM` consistent with file style).
Also write a small standalone host test (scratch dir or test harness style
used in src/Test if trivial) that fills a Matrix4 with Increment()/Random()
and a Matrix3 with distinct values, computes the product via a reference
implementation using the ORIGINAL accessor loop, and asserts exact equality
with the new operator. Keep the reference loop inside the test only.
CHECK: test passes 2D and 3D; `make` still compiles both dims.

### Step 2 - Fapply DDW hoist + column-restricted contractions (edits a,b)
DO: In src/Operator/Elastic.cpp Fapply lambda (lines ~502-671):
(a) Load `MATRIX4 const ddw = DDW(i, j, k);` once after `sten` setup; replace
    the three uses at lines ~532, ~613, ~629 (and probe use ~607) with `ddw`.
(b) Replace `(Cgrad1 * gradu).col(0) + (Cgrad2 * gradu).col(1) [+ ...]`
    (lines ~621-623) with column-restricted evaluation: compute ONLY the
    needed column of each Matrix4 x Matrix product, preserving the exact
    per-entry summation order of the existing `operator*(Matrix4, Set::Matrix)`
    row expressions (Matrix4_Major.H:540-548 for 3D; 2D analog earlier in
    file). Implement as an inline helper in Matrix4_Major.H (e.g.
    `MulCol(const Matrix4&, const Set::Matrix&, int col)` or per-column
    functions), NOT by changing the existing operator*.
CHECK: compiles both dims; quick 10-step run of the timing deck produces
finite values.

### Step 3 - Diagonal DDW hoist (edit d)
DO: In Elastic.cpp Diagonal lambda: hoist `MATRIX4 const ddw = DDW(i, j, k);`
above the `for (int p ...)` loop (~line 851); replace uses at ~860 and ~868.
CHECK: compiles both dims.

### Step 4 - Gates (oracle)
DO: run all oracle commands. All must pass. If golden compare fails: STOP,
report diff, do not tune tolerances, do not weaken tests.
CHECK: all exit 0; save logs to results/.

### Step 5 - Post-edit measurement
DO: rebuild PROFILE=1 both dims; rerun static register A/B
(`results/regs_after.txt`) and the same timing deck 3x
(`results/wall_after_runN.txt`) with clock/power snapshot.
CHECK: files exist; write comparison table (baseline vs after: regs, Fapply
incl. wall, Diagonal incl. wall, total) into results/RESULT.md draft.

## Checkpoints (tier 3)

Orchestrator (main session) reviews line-by-line diff before commit;
fresh-context verifier pass mandatory before commit (run by orchestrator,
not this agent). This agent STOPS after Step 5 and reports; it does NOT
commit.

## Adversarial review

Orchestrator spawns fresh verifier: check device-lambda captures, FP
accumulation-order preservation, Matrix4 data[] index-map correctness,
2D-vs-3D operator coverage, whether any gate was weakened.

## Closeout (orchestrator)

- [ ] Oracle passes; status.sh all green
- [ ] results/RESULT.md: what changed, evidence, deviations
- [ ] touch results/DONE
- [ ] SESSION_LOG.tsv line
