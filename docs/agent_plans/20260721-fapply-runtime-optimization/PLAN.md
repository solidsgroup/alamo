# TASK: fapply-runtime-optimization
# Folder: docs/agent_plans/20260721-fapply-runtime-optimization/

---

## Header

| Field         | Value |
|---------------|-------|
| Status        | DRAFT v2 — plan-antagonist verdict SAFE; awaiting human approval |
| Risk tier     | 3 |
| Model         | opus |
| Verification  | partial-oracle |
| Est. scope    | 2-5 source files and 3-6 benchmark/test files, split into independently judged commits |
| Parallel-safe | Read-only scouting/review only. GPU measurements and all `Elastic.cpp` mutations are serial. |

## Operating rules

1. Read ONLY the files listed in Context budget. `docs/archive/` is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: STOP,
   report the discrepancy, and wait.
3. One commit per retained optimization. Message:
   `elastic: <what> (20260721-fapply-runtime-optimization)`.
4. No scope expansion. New ideas go to `NOTES.md`; changes involving integration
   cadence, solver tolerances, or Newton require a new approval checkpoint.
5. Preserve all pre-existing work. Never clean, stage, or rewrite unrelated dirty
   files. Use explicit path lists for every diff and commit.
6. Tier 3: stop after every candidate's correctness and performance report. A
   candidate is retained only after human approval; otherwise revert only that
   candidate's isolated commit.
7. The primary device for all testing and profiling is the local NVIDIA A1000.
   Run GPU experiments serially and without unrelated GPU load. Do not use NOVA
   unless a local finalist exists and the human has supplied SSH/access.
8. Do not compare profiled and unprofiled wall times. All timing A/B runs use the
   same executable, input, environment, warmup, and repetition count.

## Context budget

Read first: `benchmark/status.sh` output, this `PLAN.md`, `docs/llm/PLAN.md`

Read:

- `src/Operator/Elastic.cpp:1-520`
- `src/Operator/Elastic.cpp:700-750`
- `src/Operator/Elastic.cpp:950-1075`
- `src/Operator/Elastic.H`
- `src/Set/Matrix4_Major.H:520-680`
- `src/Integrator/Flame.cpp:310-340`
- `src/Solver/Nonlocal/Linear.H:160-340`
- `src/Solver/Nonlocal/Newton.H:160-210`
- `src/Solver/Nonlocal/Newton.H:300-330`
- `src/Solver/Nonlocal/Newton.H:380-410`
- `src/Solver/Nonlocal/Newton.H:430-500`
- `src/Solver/Nonlocal/Newton.H:580-690`
- `benchmark/status.sh`
- `benchmark/lint_device_patterns.sh`
- `benchmark/ci_golden_compare.sh`
- `benchmark/local_a100_gate.sh`
- `benchmark/fapply_register_ab.sh`
- `benchmark/build_alamo_local_gpu.sh`
- `benchmark/compare_tinyprofiler.py`
- `benchmark/validate/README.md`
- `benchmark/validate/cases.manifest.yaml`
- `benchmark/validate/physics_budget.yaml`
- `benchmark/validate/run_validation_local.sh`
- `benchmark/validate/run_validation_local.py`
- `benchmark/validate/validation_common.py`
- `benchmark/validate/extract_metrics.py`
- `benchmark/validate/compare_validation.py`
- `tests/ElasticSoftVoid/input`
- `tests/ElasticSoftVoid/test`
- `src/Test/Operator/Elastic.H` — Step 3 focused four-state dispatch test
- `src/Test/Solver/Nonlocal/Newton.H` — focused-test registration pattern only
- `src/BC/Operator/Elastic/Elastic.H` — focused-test boundary fixture
- `src/BC/Operator/Elastic/Constant.H` — focused-test boundary fixture
- `src/Unit/Test.H`
- `src/Unit/Unit.H`
- `src/test.cc` — focused-test registration only
- `scripts/runtests.py`
- `input`
- `input_3d_centre_bore_128_a2`
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_FabArrayBase.H:640-720` — Step 7 FabArray launch metadata
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_MultiFabUtil.H:380-430` — Step 7 FabArray-wide `ParallelFor` interface
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_MFParallelForG.H` — Step 7 GPU implementation
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_MFParallelFor.H` — Step 7 FabArray-wide overloads
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_iMultiFab.H:600-620` — Step 7 owner-mask declaration
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_iMultiFab.cpp:699-780` — Step 7 owner-mask construction
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_MultiFab.cpp:1280-1370` — Step 7 nodal owner use
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_MultiFab.cpp:1475-1565` — Step 7 overlap/owner masks
- `ext/AMReX-Codes/amrex/Src/LinearSolvers/MLMG/AMReX_MLNodeLinOp.H:70-100` — Step 7 nodal owner-mask helper
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_MFIter.H:140-170` — Step 7 nodal tile-box interface
- `ext/AMReX-Codes/amrex/Src/Base/AMReX_MFIter.cpp:420-530` — Step 7 nodal tile-box construction
- `src/Operator/Operator.cpp:350-420` — Step 7 downstream nodal synchronization/consumption
- `src/Solver/Nonlocal/Newton.H:560-690` — Step 7 input/coefficient fill and synchronization

Reference only if the named step reaches its gate:

- `benchmark/NOVA_SLURM_RUNBOOK.md` — Step 9 only
- Installed AMReX FabArray-wide launch declarations found with `rg` — Step 7
  feasibility only; add exact paths to this plan before any implementation
- The exact calling/invalidation sites identified in Step 6 — add exact paths
  to this plan before implementing a psi cache

Forbidden: `docs/archive/*`, unrelated task folders, propellant parameter-sweep
artifacts, and unrelated source modules.

## Objective

Reduce end-to-end elastic MLMG GPU runtime by (a) shortening each `Fapply`
execution and (b) reducing `Fapply` calls without exceeding the established
physics error budget. The starting point is the already-merged kernel surgery:
on A100 it reduced `Fapply` wall time by 14.5% and `MLMG::solve` by 10.6%.
This campaign must determine which additional improvements survive controlled
A1000 A/B measurement and correctness gates, retain only demonstrated wins, and
leave speculative or architecture-specific ideas documented rather than merged.

## Scope and hypotheses

The campaign tests these hypotheses in order:

1. Solver configuration can reduce applications cheaply. Measure `4/4` versus
   `2/2` pre/post smoothing first; judge total solve time and physics, not iteration
   count. Tolerance, warm-start, solve interval, coarsening, and bottom-solver
   experiments are separately reported because some change solver or time-integration
   semantics.
2. A host-selected conservative kernel can avoid classification, full `gradu`,
   `psi_avg`, and cell `ddw` construction on conservative interior nodes while
   retaining the current boundary-row path.
3. The nonconservative/psi kernel can construct and contract coefficient-gradient
   tensors sequentially, reducing live ranges without changing contraction order.
4. Cached nodal `psi_avg` and `gradpsi`, correctly refreshed across AMR/MG levels,
   may trade modest memory for repeated interpolation/gradient work.
5. A FabArray-wide launch may reduce launch overhead, but it is attempted only if
   tracing shows meaningful launch overhead and an overlap-safe nodal mapping is
   proven first.

Not in the initial implementation scope: tensor cores, a new linear solver,
changes to physical equations, weakened validation budgets, or combined refactors
that prevent attribution.

## Oracle and performance contract

Every executable used by the oracle must be rebuilt from the candidate source,
and its hash/timestamp recorded; a passing stale binary is a failed gate. Step 1
freezes the exact build matrix for CPU 2D, strict-GPU 2D/3D, and performance-GPU
2D/3D (`sm_86`) and archives immutable baseline binaries/objects so candidate
builds cannot overwrite the A/B arm.

Correctness commands (all must exit 0 for every retained source change after the
required candidate builds):

```bash
: "${CANDIDATE_GPU2D_STRICT:?set exact absolute candidate binary path}"
: "${CANDIDATE_GPU2D_FAST:?set exact absolute candidate binary path}"
: "${CANDIDATE_GPU3D_STRICT:?set exact absolute candidate binary path}"
: "${CANDIDATE_GPU3D_FAST:?set exact absolute candidate binary path}"
: "${TASK_RUN_ROOT:?set a new task-local absolute artifact directory}"
: "${BASELINE_STRICT_BUNDLE_2D:?set frozen strict baseline bundle path}"
: "${BASELINE_STRICT_BUNDLE_3D:?set frozen strict baseline bundle path}"
: "${BASELINE_FAST_BUNDLE_2D:?set frozen fast baseline bundle path}"
: "${BASELINE_FAST_BUNDLE_3D:?set frozen fast baseline bundle path}"
: "${CANDIDATE_STRICT_BUNDLE_2D:?set new strict candidate bundle path}"
: "${CANDIDATE_STRICT_BUNDLE_3D:?set new strict candidate bundle path}"
: "${CANDIDATE_FAST_BUNDLE_2D:?set new fast candidate bundle path}"
: "${CANDIDATE_FAST_BUNDLE_3D:?set new fast candidate bundle path}"
test -x "${CANDIDATE_GPU2D_STRICT}"
test -x "${CANDIDATE_GPU2D_FAST}"
test -x "${CANDIDATE_GPU3D_STRICT}"
test -x "${CANDIDATE_GPU3D_FAST}"
test -d "${BASELINE_STRICT_BUNDLE_2D}"
test -d "${BASELINE_STRICT_BUNDLE_3D}"
test -d "${BASELINE_FAST_BUNDLE_2D}"
test -d "${BASELINE_FAST_BUNDLE_3D}"
test ! -e "${CANDIDATE_STRICT_BUNDLE_2D}"
test ! -e "${CANDIDATE_STRICT_BUNDLE_3D}"
test ! -e "${CANDIDATE_FAST_BUNDLE_2D}"
test ! -e "${CANDIDATE_FAST_BUNDLE_3D}"
sha256sum "${CANDIDATE_GPU2D_STRICT}" "${CANDIDATE_GPU2D_FAST}" \
    "${CANDIDATE_GPU3D_STRICT}" "${CANDIDATE_GPU3D_FAST}"

benchmark/lint_device_patterns.sh
GOLDEN_MODE=cpu benchmark/ci_golden_compare.sh
scripts/runtests.py --dim=2 --comp=g++ \
    --sections 2d-serial 2d-parallel tests/ElasticSoftVoid
BIN="${CANDIDATE_GPU3D_FAST}" TIERS='1 2' \
    benchmark/local_a100_gate.sh

SOFTVOID_COMMON=(
    max_step=2 amr.max_level=2 amr.n_cell=64 64 8
    explicitmesh.lo1=32 32 0 explicitmesh.hi1=95 95 0
    explicitmesh.lo2=72 72 0 explicitmesh.hi2=183 183 0
    pf.eta.ic.expression.constant.w=0.002
    model_void.kappa=0.2_MPa model_void.mu=0.2_MPa
)
run_softvoid_gpu () {
    local label="$1" binary="$2"; shift 2
    local out="${TASK_RUN_ROOT}/${label}"
    test ! -e "${out}"
    mkdir -p "${out}"
    "${binary}" tests/ElasticSoftVoid/input \
        "${SOFTVOID_COMMON[@]}" "$@" plot_file="${out}" \
        >"${out}/stdout" 2>"${out}/stderr"
    python3 tests/ElasticSoftVoid/test "${out}"
}
run_softvoid_gpu strict_single "${CANDIDATE_GPU2D_STRICT}"
run_softvoid_gpu strict_multibox "${CANDIDATE_GPU2D_STRICT}" \
    amr.max_grid_size=32
run_softvoid_gpu fast_single "${CANDIDATE_GPU2D_FAST}"
run_softvoid_gpu fast_multibox "${CANDIDATE_GPU2D_FAST}" \
    amr.max_grid_size=32

bash benchmark/validate/run_validation_local.sh \
    --profiles gpu_strict --case canonical_2d_elastic \
    --binary "${CANDIDATE_GPU2D_STRICT}" \
    --bundle-dir "${CANDIDATE_STRICT_BUNDLE_2D}"
bash benchmark/validate/run_validation_local.sh \
    --profiles gpu_strict --case centre_bore_3d_128_a2_converged \
    --binary "${CANDIDATE_GPU3D_STRICT}" \
    --bundle-dir "${CANDIDATE_STRICT_BUNDLE_3D}"
bash benchmark/validate/run_validation_local.sh \
    --profiles gpu_fast --case canonical_2d_elastic \
    --binary "${CANDIDATE_GPU2D_FAST}" \
    --bundle-dir "${CANDIDATE_FAST_BUNDLE_2D}"
bash benchmark/validate/run_validation_local.sh \
    --profiles gpu_fast --case centre_bore_3d_128_a2_converged \
    --binary "${CANDIDATE_GPU3D_FAST}" \
    --bundle-dir "${CANDIDATE_FAST_BUNDLE_3D}"
python3 benchmark/validate/compare_validation.py \
    "${BASELINE_STRICT_BUNDLE_2D}" "${CANDIDATE_STRICT_BUNDLE_2D}" \
    --case canonical_2d_elastic --gate --require-compatible-manifest
python3 benchmark/validate/compare_validation.py \
    "${BASELINE_STRICT_BUNDLE_3D}" "${CANDIDATE_STRICT_BUNDLE_3D}" \
    --case centre_bore_3d_128_a2_converged --gate \
    --require-compatible-manifest
python3 benchmark/validate/compare_validation.py \
    "${BASELINE_FAST_BUNDLE_2D}" "${CANDIDATE_FAST_BUNDLE_2D}" \
    --case canonical_2d_elastic --gate --require-compatible-manifest
python3 benchmark/validate/compare_validation.py \
    "${BASELINE_FAST_BUNDLE_3D}" "${CANDIDATE_FAST_BUNDLE_3D}" \
    --case centre_bore_3d_128_a2_converged --gate \
    --require-compatible-manifest
```

The local gate must be confirmed to use the A1000 before execution. Tier 2 is the
compute-sanitizer/multi-box smoke, but it deliberately uses 1/1 smoothing, one
Newton iteration, at most two MLMG iterations, and its Tier 1 stops before the
elastic solve. It proves memory safety, not normal solver semantics. Step 1 adds
the displayed `--binary` and `--bundle-dir` options to the local validation runner,
rejects non-absolute/missing binary paths and pre-existing bundle paths, and records
binary path/hash, build command/flags, HEAD, scoped-source-diff hash, GPU UUID,
profile/architecture, full executed commands, input/override hashes, and hashes of
the case manifest, physics budget, extractor, comparator, and property-oracle
scripts in each manifest. It also adds the displayed manifest-compatibility
preflight to `compare_validation.py`. The preflight rejects any baseline/candidate
mismatch in device, profile/arch, inputs/overrides, command shape, or oracle/parser
identity; only timestamp/output path and the explicitly compared binary/source
identity may differ. Step 1 also replaces all bundle/run
variables above with frozen absolute artifact paths before any edit. The 3D
comparison is candidate versus the same-GPU pre-change bundle: the repository
already documents a CPU-versus-GPU field-norm discrepancy for this case, so using
the immutable CPU reference would confound a new optimization with an existing
baseline difference.

Covers: device-lambda/lifetime lint, strict- and fast-build normal-setting behavior,
local A1000 execution of both conservative and psi paths, compute-sanitizer
coverage, multi-box behavior, and physics-budget validation. `scripts/runtests.py` itself
selects CPU-style binary names; it is not accepted as evidence that the
conservative path executed on the A1000. Any psi-cache or launch-mapping candidate
must additionally run `BIN="${CANDIDATE_GPU3D_FAST}" TIERS='1 2 3'
benchmark/local_a100_gate.sh`; Tier 3 supplies racecheck and initcheck.

Does NOT cover: every nodal overlap ordering, all AMR/MG hierarchy invalidation
orders, statistical performance certainty, or A100 architecture behavior. Those
gaps require focused tests, adversarial review, and (for a shipping finalist) a
separate NOVA confirmation.

Performance contract:

- Primary solve metric: median inclusive `MLMG::solve` region wall from matched
  profiler-enabled baseline/candidate runs. Co-primary application metric: median
  external end-to-end wall from matched unprofiled runs. Never compare the
  profiler-enabled value to the unprofiled value.
- Secondary metrics: TinyProfiler inclusive `Fapply` region wall and `NCalls`,
  time per call, nonlinear/MLMG iterations, residuals, registers, local/stack
  bytes, spills/replays, achieved occupancy, launch count, and summed GPU kernel
  duration from Nsight. Only the Nsight-derived kernel duration may be labeled
  exclusive device time.
- Regimes: a 2D conservative (`elastic.use_psi=0`) case and a 3D nonconservative
  (`elastic.use_psi=1`) case that both fit on the A1000. Step 1 fixes exact inputs.
- Timing protocol: one warmup plus at least five measurements per arm and mode,
  interleaved A/B when practical. Record clock/power/P-state and concurrent GPU
  processes. Profiler-enabled A/B runs use identical synchronization settings;
  unprofiled A/B runs use none of them.
- Keep criterion: correctness passes, neither primary metric regresses, and the
  solve gain exceeds the larger of 3% or twice the baseline median absolute
  deviation. A smaller candidate may be retained only if the human explicitly
  accepts it as an enabling change.
- Call-reduction experiments must report total solve time and physics output;
  fewer calls or iterations alone are never a pass.

## Multi-agent workflow

At most three subagents are used over the entire campaign:

1. **Scout / bounded worker (cheap):** read-only reconnaissance during planning;
   after approval, owns one explicitly assigned file set at a time for mechanical
   harness work, isolated edits, and commanded test runs. It never makes keep/revert
   decisions and never overlaps another writer.
2. **Plan antagonist (reviewer):** receives the draft plan and assumes it is wrong;
   attacks attribution, numerical equivalence, hierarchy lifetime, overlap safety,
   and measurement validity. It makes no edits and ends before implementation.
3. **Fresh implementation antagonist (reviewer):** spawned only after all local
   finalists pass. It receives commit hashes and the approved plan but no worker
   narrative; findings go to `results/REVIEW.md`.

The orchestrator retains complex reasoning: experiment design, source-path
decisions, numerical-order adjudication, GPU scheduling, keep/revert decisions,
integration, and all human checkpoints. CPU/read-only work may overlap, but no two
agents profile the GPU concurrently and no two agents edit `Elastic.cpp`.

## Plan-antagonist disposition

The independent plan reviewer judged v1 unsafe to present. This v2 incorporates
its critical/major findings: exact binary paths and hashes replace wildcard/latest
selection; the conservative A1000 single-/multi-box oracle is executable; sanitizer
tiers are described by their real coverage; the absent launch-bounds helper was
removed; conservative dispatch and mixed psi states are explicit; psi cache
generation/refresh order and memory headroom are gated; launch fusion requires a
pre-sync writer-count/value oracle plus racecheck; and profiler call/time semantics
are frozen. A second hostile pass then found that fast production binaries lacked
normal-setting semantic validation and that metrics comparison did not enforce
manifest identity; v2 now gates both strict and fast builds in both regimes and
requires manifest compatibility before comparing metrics. The A100 step remains a
transparent pre-ship gate because the live repository plan requires it, but it
cannot begin without the user's post-local authorization.
The final hostile audit found no remaining critical or major plan defects and
returned **SAFE** for human presentation.

## Steps

### Step 0 — Human plan approval and workspace freeze

VERIFY:

```bash
benchmark/status.sh
git status --short -- src/Operator src/Set src/Solver benchmark tests \
    docs/agent_plans/20260721-fapply-runtime-optimization
git branch --show-current
git rev-parse HEAD
git config --get core.hooksPath
```

DO: Record HEAD, branch, scoped dirty paths, compiler/CUDA/AMReX identity, and
A1000 identity in `NOTES.md`. Do not disturb the existing dirty worktree. The
launch-bounds helper named by the live plan is not present on this branch; do not
merge or recreate that separate campaign during an A/B series. The current
`status.sh` may remain nonzero because of pre-existing unrelated dirty artifacts;
this checkpoint passes only when its substantive gate output is recorded, device
lint passes, and the explicitly scoped source diff is clean.

CHECK: Human approves this plan in writing. This is the hard stop for the current
planning session; no implementation, build, profiling, or benchmark run precedes
approval.

### Step 1 — Freeze cases, runner commands, and baseline

VERIFY:

```bash
git diff --exit-code -- src/Operator/Elastic.cpp src/Operator/Elastic.H \
    src/Set/Matrix4_Major.H
GPU_QUERY=name,uuid,driver_version,memory.total,compute_cap,pstate
GPU_QUERY="${GPU_QUERY},clocks.current.sm,clocks.current.memory,power.limit"
nvidia-smi --query-gpu="${GPU_QUERY}" \
    --format=csv,noheader
nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory \
    --format=csv,noheader
```

DO: The orchestrator chooses exact A1000-sized conservative 2D and psi-enabled 3D
inputs from existing harnesses. The bounded worker may add only missing benchmark
runner/parsing support under `benchmark/`, with explicit ownership. Its first
bounded change is to make `build_alamo_local_gpu.sh` produce/copy an exact requested
artifact and make `run_validation_local.py` accept that exact absolute binary path
and a new exact bundle directory; wildcard/latest-binary discovery and implicit
timestamp discovery are forbidden for campaign evidence. The manifest
records the full identity listed in the Oracle section, and the comparator gets a
manifest-compatibility preflight that is self-tested with one valid pair and
deliberately mismatched device/profile/input/oracle-hash pairs. Record exact build
and validation commands, fixed parameters, repetitions, profiler versions, and
output paths. Archive immutable strict and fast baseline binaries/objects and their
hashes so later candidate builds cannot silently replace the baseline arm. Establish baseline calls,
timing distribution, physics output, register/stack/spill data, and a launch trace.
The conservative correctness/performance input is based on
`tests/ElasticSoftVoid/input`; `canonical_2d_elastic` remains a psi-path validation
case and must not be mislabeled conservative. Freeze a direct strict-GPU invocation
of both the `2d-serial` and chopped/multi-box `2d-parallel` ElasticSoftVoid layouts,
followed by `tests/ElasticSoftVoid/test` on each output. If a small task-local
wrapper is needed, it is the worker's only file ownership for this substep and is
validated against the section arguments before use.

Build a profiler-enabled A1000 executable with `PROFILE=1` and freeze its exact
path/hash too. Raw logs must contain the AMReX TinyProfiler row
`Operator::Elastic::Fapply()` and its `NCalls`. Step 1 adds a task-local parser or
extends `compare_tinyprofiler.py` to report the count and inclusive region wall
explicitly; it may not call that value exclusive time. External timed
whole-run wall is the unprofiled application metric; matched TinyProfiler
`MLMG::solve` wall is the solve metric, and Nsight data is diagnostic attribution
only.

CHECK: Repeating the baseline produces stable results within the recorded noise
envelope; all correctness oracle commands pass before source experiments.

Checkpoint report: case identities, command transcript, environment, baseline
table, noise/MAD, known coverage gaps. Human approves the frozen baseline.

### Step 2 — Reduce calls through configuration experiments

VERIFY: Step 1 baseline and physics outputs exist; no source modification is
present.

DO: On the conservative and psi regimes, first sweep pre/post smoothing `4/4`
versus `2/2` (optionally add `3/3` only if it clarifies a non-monotonic result).
Report `Fapply` calls, iterations, total solve time, and physics-budget results.
Then, as separate labeled experiments, evaluate coarsening depth and bottom solver
only where current inputs expose supported knobs. Pass every experimental setting
as a command-line override; do not edit canonical input decks. Do not combine
changes until each has an individual result.

Tolerance relaxation, previous-displacement warm start, and elasticity solve
interval are report-only feasibility experiments in this task. Because they can
change convergence or time-integration semantics, changing defaults or source code
requires a new human checkpoint and an amended context budget.

CHECK: Re-run the unmodified `4/4` arm after the sweep to detect drift. Retain a
configuration only if it meets the performance contract and physics budget.

Checkpoint report: full sweep table and recommended configuration, including any
case where fewer applications made total solve slower. Human decides whether any
configuration change proceeds beyond experiment status.

### Step 3 — Conservative-kernel specialization

VERIFY: Step 1 shows the conservative kernel is still material on A1000; scoped
source paths are clean at the accepted baseline commit.

DO: In `Elastic::Fapply`, select conservative versus general device lambdas on the
host using exactly `m_conservative_face_flux`—not `elastic.use_psi` or
`m_psi_set`—as the dispatch predicate. The conservative variant keeps identical
tile coverage and one launch per existing tile. It determines whether a row is a
boundary row before constructing boundary-only state; conservative interior rows
construct only the face-gradient/face-`ddw` quantities they actually consume.
Boundary rows must retain the current `m_psi_set`/`psi_avg`, displacement,
gradient, center-`DDW`, and BC behavior. Do **not** split interior and boundary
into separate launches, because nodal grown-box overlap and write ordering have
not yet been proven safe. Preserve current boundary conditions, face-gradient
operations, and output ownership. Before implementation, define a focused operator
test for the four representable states `{conservative on/off} × {psi set/unset}`;
if a state is truly unsupported, prove and enforce that invariant before launch
instead of silently dispatching by `use_psi`. Add the exact focused-test path to
the Context budget at the pre-edit checkpoint.

CHECK:

```bash
git diff --check -- src/Operator/Elastic.cpp
# Rebuild every required candidate executable, then run the full Oracle section.
```

Then perform isolated A/B timing and register/stack/occupancy profiling on the 2D
conservative case, plus a no-regression run on the psi case.

Checkpoint report: line-by-line diff, oracle output, per-call and solve statistics,
compiler resource table, and boundary/multi-box evidence. Human says keep or revert.

### Step 4 — Sequential coefficient-gradient construction

VERIFY: Step 1 shows register/local-memory pressure or per-call cost remains
material in the 3D psi case. Start from the accepted Step 3 baseline (or original
baseline if Step 3 was reverted).

DO: Construct and contract `Cgrad1`, `Cgrad2`, and `Cgrad3` in separate scopes so
only one tensor is live at a time. Preserve the mathematical and floating-point
addition order of the current expression; explicitly materialize only the smaller
vector accumulator needed to maintain `(v1 + v2) + v3` before the `psi_avg`
scaling. Do not mix this with psi caching or other Matrix4 changes.

CHECK: Run the full correctness oracle, isolated A/B timing on the 3D psi case,
compiler resource reporting, and a conservative no-regression run. Report any
bitwise differences separately from physics-budget acceptance rather than labeling
them automatically harmless.

Checkpoint report: diff, expression-order explanation, golden deltas, resource
table, and end-to-end timing. Human says keep or revert.

### Step 5 — Combined retained-wins evaluation

VERIFY: Every retained change has an isolated commit and isolated A/B report.

DO: Benchmark the accepted configuration and kernel commits individually and
together. Attribute improvements to call count versus time per call. Check for
negative interactions such as extra iterations, cache pressure, or a conservative
win hiding a psi regression.

CHECK: Full correctness oracle on the combined tree and the complete A1000 timing
protocol on both regimes.

Checkpoint report: additive/non-additive decomposition and recommended merge set.
Human approves the base for larger experiments.

### Step 6 — Conditional nodal-psi cache design and experiment

VERIFY: Profiling shows `CellToNodeAverage`/`CellGradientOnNode` work is still
material after Step 4, and the A1000 memory budget can hold per-AMR/MG nodal caches
with ghosts. The gate is numerical, not subjective: projected baseline peak plus
cache bytes must be below both 85% of physical memory and physical memory minus
1 GiB; the trial run must retain at least 1 GiB observed headroom and show no OOM
or managed-memory thrashing. If not, document and skip.

DO (design gate first): map allocation, fill, restriction, ghost-fill,
invalidation, define/regrid, and hierarchy-lifetime paths. `Elastic::SetPsi` has no
null overload: explicitly decide whether no in-process psi→no-psi transition is
supported or whether a separately reviewed `ClearPsi()` API is required. Do not
invent `SetPsi(nullptr)`. Use a per-AMR/MG generation scheme invalidated by every
`SetPsi`, define, and regrid; populate the cache only after `averageDownCoeffs` has
completed MG restriction and psi boundary fill/synchronization. Add the exact
source/test paths to Context budget before implementation. Account explicitly for
the separate Newton RHS use of psi quantities and repeated `SetPsi` during
line-search preparation; do not assume the `Fapply` cache can be shared before
`SetPsi`. Specify cache centering, components, ghost width, exact byte footprint,
and a focused cached-versus-uncached test that mutates the same psi field between
solves, exercises line-search re-preparation, AMR levels, MG levels, and both
single-/multi-box layouts.

Only after human approval of that design may the worker implement an Elastic-owned
cache. The implementation is a separate commit and is judged by the full oracle,
the psi-change/hierarchy test, A1000 memory high-water mark, 3D psi A/B timing, and
conservative no-regression timing.

Checkpoint report: lifetime diagram, memory calculation, invalidation test, and
performance table. Human says keep, revise, or abandon.

### Step 7 — Conditional multi-box launch feasibility gate

VERIFY: A dedicated Nsight Systems trace shows launch overhead is material after
all retained simpler work. Quantify it; do not infer it from launch count alone and
do not compare the traced wall time with an unprofiled wall measurement.

DO (read-only design spike): locate the exact installed AMReX FabArray-wide launch
API and add its paths to Context budget. Prove how `grownnodaltilebox()` coverage,
nodal overlaps, boundary masks, and duplicate writes map to the new launch. Build a
per-index writer-count **and value** oracle before proposing an implementation.
The oracle must compare the unsynchronized `a_f` result on decomposed single- and
multi-box layouts, because final plot fields after synchronization can hide duplicate
non-atomic writes. Specify a disjoint-owner mapping before launch fusion. If output
rows cannot be made disjoint or pre-sync equivalence cannot be executable-tested,
stop and record the blocker.

CHECK: The design demonstrates exactly one writer per output index, pre-sync value
equivalence, Tier 3 racecheck/initcheck success, and an upper bound on recoverable
launch time.

Checkpoint report: trace evidence, API mapping, coverage proof/test, and estimated
payback. Because this is a structural refactor, implementation requires a separate
tier-3 child plan and explicit human approval.

### Step 8 — Fresh antagonistic implementation review

VERIFY: All candidate commits and local result bundles are complete; no author is
still modifying the reviewed paths.

DO: Spawn the fresh implementation reviewer with commit hashes and this prompt:

> Assume these FApply optimizations contain a defect or a misleading benchmark.
> Find it. Check device-lambda captures and elixir lifetimes; conservative boundary
> rows; nodal grown-box ownership; exact coefficient-gradient accumulation order;
> psi cache freshness across SetPsi, restriction/boundary fill, line-search
> re-preparation, regrid, and any supported clear/no-psi transition; weakened
> tests or budgets; profiler/timing contamination; and whether fewer FApply calls
> actually reduce total solve time. Report findings only.

CHECK: Save findings to `results/REVIEW.md`. The orchestrator reproduces or disproves
each finding with evidence. Any code change triggered by review repeats its isolated
oracle and A/B gate.

Checkpoint report: findings, adjudications, fixes, and residual risks. Human
approves or rejects the local finalist.

### Step 9 — Pre-ship NOVA/A100 confirmation (conditional on access)

VERIFY: A local A1000 finalist exists, all local gates and adversarial review pass,
and the human has provided SSH/NOVA access. Otherwise stop with a locally validated
finalist that is explicitly **not cleared to ship**; missing remote access is not a
failure of the local campaign, but the live repository plan requires A100 A/B
evidence before a kernel optimization lands.

DO: Follow `benchmark/NOVA_SLURM_RUNBOOK.md` for baseline/finalist A/B wall timing,
Nsight resource/occupancy evidence, and physics compare. Do not tune on NOVA before
the local finalist is frozen. Never request or start this remote step until the
human reviews the local result and explicitly authorizes NOVA use.

CHECK: A100 results confirm correctness and no material end-to-end regression.
Architecture-specific differences are documented; they do not overwrite A1000
results.

### Step 10 — Closeout

VERIFY: Only approved commits are present; all rejected experiments are absent;
scoped status is understood.

DO: Write `results/RESULT.md` with retained and rejected changes, commands, raw
artifact paths, call-count/time decomposition, validation evidence, adversarial
adjudication, and remaining ideas. Update the live plan by hand without exceeding
its 100-line limit. Append (never rewrite) release/process records where required.

CHECK:

```bash
benchmark/status.sh
git diff --check
```

## Checkpoints

- [ ] Human approves this draft before any implementation, build, or profiling.
- [ ] Baseline cases, exact commands, and noise envelope are approved.
- [ ] Before each source edit, the worker restates ownership and approach.
- [ ] Before each commit, provide a line-by-line diff summary and oracle output.
- [ ] After each isolated candidate, human decides keep/revert.
- [ ] Psi cache design receives separate approval before implementation.
- [ ] Multi-box conversion, if justified, receives a separate child plan.
- [ ] Fresh adversarial findings are adjudicated before final acceptance.
- [ ] NOVA is used only after access is supplied for a local finalist.

## Closeout checklist

- [ ] Correctness oracle passes on every retained commit and the combined tree.
- [ ] A1000 performance contract passes with raw artifacts retained.
- [ ] `results/RESULT.md` records changes, evidence, and deviations.
- [ ] `results/REVIEW.md` records fresh adversarial findings/adjudications.
- [ ] Changelog entry is append-only; `VERSIONS.md` changes only if release-worthy.
- [ ] `results/DONE` exists only when all approved scope is complete.
- [ ] One session-log line is appended with date, task, model, tier, cost, outcome.
