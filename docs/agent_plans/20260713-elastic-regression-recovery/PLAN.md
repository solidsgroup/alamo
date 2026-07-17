# TASK: elastic-regression-recovery
# Folder: docs/agent_plans/20260713-elastic-regression-recovery/

---

## Header

| Field | Value |
|---|---|
| Risk tier | 3 |
| Model | GPT-5 |
| Verification | partial-oracle, escalating to full regression |
| Est. scope | Elastic/Newton/stencil contracts, focused/full CPU tests, 2-D/3-D CPU/CUDA chamber validation, performance evidence, and task records |
| Parallel-safe | no; all likely fixes share Elastic/Newton paths |

## Operating rules

1. Preserve the current failed artifacts and do not modify the original WIP
   worktree while constructing a baseline comparison.
2. Before a source edit, reproduce the relevant failure from the handoff and
   state the invariant the change will restore.
3. Do not use a psi floor, artificial stiffness, tolerance relaxation,
   iteration-cap increase, forced nonlinear acceptance, or reference refresh
   as a substitute for a correctness fix.
4. Keep residual, Fapply, Diagonal, smoother, and boundary closure consistent
   with one discrete operator.
5. Record every baseline and failed experiment in this task's NOTES.md. Do
   not commit unless the user explicitly asks for one.
6. Tier-3 checkpoints apply before each source change and before any proposed
   commit.

## Context budget

Read first: benchmark/status.sh output, this PLAN.md,
  ../20260709-elastic-void-robustness/results/REGRESSION_RECOVERY_HANDOFF.md,
  ../20260709-elastic-void-robustness/results/REVIEW.md.
Read: ../20260709-elastic-void-robustness/NOTES.md,
  src/Operator/Elastic.{H,cpp}, src/Solver/Nonlocal/Newton.H,
  src/Numeric/Stencil.H, src/Operator/Operator.cpp,
  src/Test/Numeric/Stencil.H, tests/ElasticSoftVoid/{input,test}.
Reference only if the relevant step names it:
  report/output_2026-07-13_15.29.01_kermit.{json,html},
  tests/SCPSpheresElastic/, tests/PlateHole/, tests/FracturePFCZM/,
  tests/RubberPressurizedHole/, historical June report JSON files.
Forbidden: docs/archive/* and unrelated task folders.

## Objective

Identify and repair the correctness regressions exposed by the full GCC suite
after the paired elastic-stencil work. Establish whether each failure belongs
to the current unstaged delta, the WIP base, or an independent environment
condition before changing solver code. Restore numerical correctness without
reintroducing hidden regularization. Validate the finite-soft-void chamber in
2-D and 3-D on CPU and CUDA, including a realized sub-MPa effective material
screen and reproducible local performance measurements, and leave a complete
validation record.

## Handoff status (2026-07-16)

**BLOCKED — source repairs and focused tests are preserved, but the final
finite-soft-void equilibrium oracle fails in both dimensions.  Do not commit
or create `results/DONE`.**

| Step | Status | Evidence |
|---|---|---|
| 1. Same-base attribution | Complete | Parent `d964cfab8` passes the trusted mechanics cases; clean paired rewrite `a06c6f15d` introduces false-zero/NaN/mechanics drift. |
| 2. First violations | Complete | First-order paired mixed derivative, Constant Neumann integer spacing, missing periodic/model mixed-corner halos, and diagonal/normalize ghost mismatch were isolated independently. |
| 3. Contract repairs | Partial | Parent operator restored; Newton, Constant BC, periodic model, diagonal, and normalize fixes compile and pass focused probes.  Residual/tangent consistency for the final soft-void AMR case remains unresolved. |
| 4. Oracle gaps | Partial | Unit and strengthened regression oracles exist; fresh 2-D serial/MPI and retained 3-D runs reach only the plotted equilibrium gate and fail there. |
| 5. Complete solver surface | Blocked | No final-source complete 2-D/3-D suite or fresh adversarial review. |
| 6. CPU/CUDA performance | Blocked | No accepted correctness case, CUDA toolkit is not currently available, and no timing campaign was run. |

Canonical closeout records are `results/RESULT.md`, `results/REVIEW.md`, and
`results/HANDOFF.md`.  `NOTES.md` remains the chronological experiment
almanac.

## Oracle

Command(s):

~~~
scripts/runtests.py --comp=g++ --no-backspace --no-clean
scripts/runtests.py --comp=g++ --no-backspace --no-clean --sections 2d-parallel-long tests/SCPSpheresElastic
scripts/runtests.py --comp=g++ --no-backspace --no-clean --sections 2D-serial tests/PlateHole
scripts/runtests.py --comp=g++ --no-backspace --no-clean --sections 2d-serial-notch tests/FracturePFCZM
./configure --dim=2 --comp=g++ --no-debug --cuda=86 --cuda-fp=strict --gpu-integrator=flame --build-amrex --offline
./configure --dim=3 --comp=g++ --no-debug --cuda=86 --cuda-fp=fast --gpu-integrator=flame --build-amrex --offline
~~~

Covers: execution failures, checked numerical regressions, serial/MPI
consistency, the finite-soft-void targeted path, CUDA correctness, and local
rank-one CPU/CUDA chamber performance.

Does NOT cover: a raw exact-zero psi formulation, independent experimental
constitutive validation, production-scale GPU scaling, or whole-node CPU/GPU
throughput. CUDA execution begins only after the CPU correctness issue is
resolved and a working GPU environment is recorded.

## Steps

### Step 1 - make same-base baseline evidence

VERIFY: the dated full-run report and preserved failure outputs exist.

DO: create separate clean worktrees for the parent/pre-WIP baseline and for
a clean a06c6f15d build. Reproduce SCPSpheresElastic, PlateHole, and
FracturePFCZM in each, recording compiler/configuration and first divergent
step.

CHECK: NOTES.md contains a table that distinguishes current-delta failures,
WIP-base failures, and unresolved cases.

### Step 2 - isolate the first correctness violation

VERIFY: one target has a deterministic first bad linear or nonlinear solve.

DO: add temporary diagnostics or a minimal operator probe outside the
acceptance test only as needed to compare residual directional derivatives,
Fapply, Diagonal impulse response, finite values, ghost availability, and
boundary/C-F behavior.

CHECK: the experiment yields a falsifiable cause, not just a parameter that
changes iteration count. Remove diagnostic-only code before the next step.

### Step 3 - repair one proven contract at a time

VERIFY: Step 2 identifies the violated discrete or data-lifecycle invariant.

DO: change the smallest Elastic/Newton/stencil code that restores that
invariant. Treat the direct SetPsi three-cell halo contract and masked
theta*DDW AMR restriction as separate hypotheses unless evidence combines
them.

CHECK: the minimal reproducer passes with the original tolerances and the
change has a direct low-level or regression test.

### Step 4 - close oracle gaps

VERIFY: repaired target passes and its behavior is understood.

DO: add durable tests for the repaired contract. Strengthen ElasticSoftVoid
to assert its material values, automate MPI coverage, prove C/F
heterogeneity, and add a mechanics-sensitive assertion. Ensure
SCPSpheresElastic rejects non-finite plots and false nonlinear success.

CHECK: the new tests fail on the demonstrated bad state and pass on the
repaired state without weakening an existing reference.

### Step 5 - validate the complete solver surface

VERIFY: all targeted reproducers and new tests pass.

DO: rebuild 2-D and 3-D GCC targets, run the complete suite, and have a fresh
agent adversarially review the final diff for data-hierarchy freshness,
residual/tangent/diagonal agreement, boundary consistency, and test-oracle
weakening.

CHECK: every existing run/check section passes; make test also passes when
the local doxygen/libclang dependency is repaired.

### Step 6 - validate chamber CPU/CUDA correctness and performance

VERIFY: final 2-D/3-D CPU suites pass and the CUDA toolkit/device provenance is
captured.

DO: build strict and fast CUDA Flame targets. Run plotted CPU, strict-CUDA, and
fast-CUDA chamber correctness cases with `elastic.use_psi=0` and
`elastic.psi_floor=0`. Screen the interface width until the plotted effective
bulk and shear moduli are below 1 MPa without weakening solver tolerances.
Then run one excluded warm-up and three rank-one timed repetitions per backend
and dimension with exactly five elastic solves and no plot I/O.

CHECK: correctness plots and solver diagnostics are finite and converged;
strict CPU/CUDA field metrics satisfy the repository comparison gate; every
timed log has equal work and accepted residual/update histories; raw times,
medians, speedups, iteration counts, and peak GPU memory are preserved with
binary/configuration hashes. Report the reduced 3-D case and rank-one/local
scope explicitly.

## Checkpoints

- [x] After plan restatement: the assignment explicitly authorized the
      baseline-first investigation before the recovery source edit.
- [ ] Before each source edit: minimal reproducer, measured invariant, and
      expected effect recorded in NOTES.md. The discarded FillPatch ghost
      experiment was recorded after, not before, its edit and remains an
      explicit process deviation.
- [x] Before each commit: N/A; the user did not request a commit and none will
      be created.
- [ ] Before closeout: full-suite report and fresh adversarial review attached.

## Adversarial review

After any candidate fix passes the complete suite, use a fresh Tier-3 reviewer
with no implementation context. The reviewer must inspect the final diff and
test changes for hidden regularization, stale AMR/MG coefficients, unsafe
Array4 accesses, residual/Jacobian/Diagonal mismatch, boundary inconsistency,
GPU lifetime/capture issues, and weakened tests. Record findings in this
task's results/REVIEW.md.

## Closeout

- [x] Targeted root-cause evidence and each deviation recorded in NOTES.md
- [ ] All existing GCC regression runs and checks pass
- [ ] New tests prove repaired contracts and no oracle was weakened
- [ ] 2-D/3-D CPU/CUDA chamber correctness and performance evidence complete
- [x] results/RESULT.md and results/REVIEW.md record evidence and limitations
- [ ] make test passes when the documentation environment is available
- [ ] touch results/DONE only after every gate passes
- [x] Session log line appended with the actual final state
