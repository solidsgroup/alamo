# TASK: gpu-manual-reliability
# Folder: docs/agent_plans/20260721-gpu-manual-reliability/

---

## Header

| Field         | Value |
|---------------|-------|
| Risk tier     | 0 |
| Model         | opus for oracle/schema judgment; sonnet for deterministic scripts |
| Verification  | partial-oracle |
| Est. scope    | 12-18 manual/script/CSV files; no `src/` edits |
| Parallel-safe | no; schema, scanner, labels, and generated ledgers share contracts |

## Operating rules

1. Read only the files listed in Context budget. `docs/archive/` is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: stop and report the discrepancy.
3. Keep source branches read-only. All writes stay in `docs/gpu_manual/**` and this task's `results/`.
4. Do not import FEATURE, `[NUM]`, or performance changes into a correctness transform.
5. Do not call a scan, transform set, port, or optimization complete from a partial gate.
6. Preserve the untracked `docs/gpu_manual.zip`; it is user-owned and out of scope.
7. Treat `chamber-gpu` diffs outside their tested integrator closure as candidate evidence. Strict selected-closure compiler/runtime results outrank them.

## Context budget

Read first: `benchmark/status.sh` output, this `PLAN.md`
Read: `docs/gpu_manual/BRIEF.md`, `INDEX.md`, `recognizers/table.csv`, `recognizers/scan.py`, `COVERAGE.csv`, `build/phase6/GATE_RUNS.csv`, `build/phase6/COUNTS.txt`, `build/validate_patterns.py`, all `patterns/GPU-*.md`
Reference only when labeling evidence: `build/HUNK_MAP.csv`, exact supporting BASE/chamber-gpu diff hunks, and the corresponding target file
Forbidden: `docs/archive/*`, unrelated task folders, source changes, and chamber-gpu implementation outside exact evidence adjudication

## Objective

Replace the current instance-fitted, positive-only recognizer workflow with a site-level, held-out-validated system that distinguishes transform correctness, recognizer quality, end-to-end port completion, baseline GPU efficiency, and later profile-driven optimization. Expand the manual plan so manual/build-error patterns enter the port plan, scanner reruns converge, correctness validation covers distinct failure classes, and a correct port cannot pass while retaining obvious hot-path host fallback, transfer, synchronization, dispatch, or launch-multiplication pathologies.

## Oracle

Commands after their implementing steps:

```bash
python3 docs/gpu_manual/build/validate_patterns.py
python3 docs/gpu_manual/recognizers/validate.py \
  --table docs/gpu_manual/recognizers/table.csv \
  --labels docs/gpu_manual/recognizers/labels.csv \
  --metrics docs/gpu_manual/recognizers/metrics.csv
python3 docs/gpu_manual/recognizers/scan.py \
  --root docs/gpu_manual/build/gpu-scan \
  --table docs/gpu_manual/recognizers/table.csv \
  --out docs/gpu_manual/build/target-sites.csv
python3 docs/gpu_manual/build/build_port_plan.py \
  --sites docs/gpu_manual/build/target-sites.csv \
  --inspections docs/gpu_manual/INSPECTIONS.csv \
  --out docs/gpu_manual/PORT_PLAN.csv
```

Covers: v2 schema, labeled-site accounting, known-defect recall, post-state convergence, mixed-state reporting, manual inspection inclusion, and honest per-pattern status reporting.

Does not cover: unknown defects outside the labeled corpus, correctness of human confirmation answers, transitive CUDA call chains beyond compiler diagnostics, runtime races, multi-box behavior, or performance. Those require the compiler/runtime matrix in Step 6.

## Steps

### Step 1 - Freeze the v1 failure corpus

VERIFY:

```bash
test "$(awk -F, 'NR>1 && $2=="manual" {n++} END {print n+0}' docs/gpu_manual/recognizers/table.csv)" -eq 8
grep -q 'VALIDATED_PATTERNS=GPU-007;GPU-016' docs/gpu_manual/build/phase6/COUNTS.txt
```

DO: Record every known v1 miss/self-match as a labeled regression, including all five BASE `Operator/Elastic.cpp` GPU-007 sites, the intentional GPU-016 wrapper implementation, dotted and arrow receiver forms for GPU-023, renamed-variable variants for GPU-013/019/030, and mixed converted/unconverted files. GPU-007 candidates must not depend on the identifier `DX`, an exact scalar spelling, or a single-subscript geometry expression. GPU-023 candidates include `.` and `->`, but confirmation must distinguish `BaseFab`/`FArrayBox` explicit execution from GPU-aware `MultiFab` operations. Labels must include exact tree/ref, file, line, state, and evidence.

CHECK: every regex pattern has a positive and converted label; high-risk patterns also have intentional/negative labels; labels are deterministic and duplicate-free.

### Step 2 - Migrate to the v2 pattern and recognizer schema

VERIFY:

```bash
python3 docs/gpu_manual/build/validate_patterns.py
```

DO: Split `Status` into `Transform status`, `Detection type`, and `Recognizer status`; add Candidate, Converted, Exclusions, and Confirmation. Migrate `table.csv` to candidate/converted/exclusion fields and update the validator. Manual/build-error patterns use recognizer status `not-applicable`. Do not mark a migrated pattern verified unless its corresponding old gate proves that specific transform; no regex recognizer inherits verification.

CHECK: all 25 active patterns validate; only GPU-007/GPU-016 retain transform verification; no recognizer is held-out-pass before Step 4.

### Step 3 - Implement site-level scanning and explicit inspections

VERIFY:

```bash
python3 -m py_compile docs/gpu_manual/recognizers/scan.py
```

DO: Emit line-addressable candidate/converted results and file-pattern states `needs-work`, `converted`, `mixed`, and `inspect`. Apply path/type exclusions without hiding other sites. Generate `INSPECTIONS.csv` rows for all manual/build-error patterns using concrete anchors and questions. Generate `PORT_PLAN.csv` from both sources; retain `COVERAGE.csv` only for compatibility.

CHECK: the GPU-016 wrapper is intentional, a mixed fixture reports mixed, pointer receivers appear as candidates without being presumed defects, and all eight manual plus GPU-001 inspection classes appear in the plan.

### Step 4 - Add the held-out recognizer oracle

VERIFY:

```bash
python3 -m py_compile docs/gpu_manual/recognizers/validate.py
```

DO: Compare scan sites with frozen labels and emit TP/FP/FN, precision, recall, and status per pattern. Correctness recognizers require zero known false negatives and recall 1.0. Keep precision visible but never trade away a labeled correctness site to improve it.

CHECK: deliberately restoring each v1 regex makes the oracle fail; GPU-007 finds all labeled pointer shapes; GPU-016 converges without self-reporting its wrapper as a defect.

### Step 5 - Replace the convenience-file closed-book gate

VERIFY:

```bash
test "$(cut -d, -f2 docs/gpu_manual/build/phase6/GATE_RUNS.csv | tail -n +2 | sort -u | wc -l)" -eq 3
```

DO: Create `VALIDATION_MATRIX.csv` covering device call chain, static value dispatch, BC dispatch, storage lifetime, Elixir lifetime, reductions, dimension guards, explicit Fab execution, and FFT API routing. A fresh session validates each transform before that transform becomes verified. Report `TRANSFORMS_VALIDATED=x/25` and `RECOGNIZERS_VALIDATED=y/25`; remove unqualified `GATE=pass` language.

CHECK: the matrix cannot turn green from GPU-007/GPU-016-only runs, and draft status is visible beside each unexercised pattern at point of use.

### Step 6 - Make compiler/runtime completion explicit

VERIFY:

```bash
grep -q 'device-lint.*PASS' docs/agent_plans/20260721-gpu-manual-reliability/results/STATUS_BASELINE.txt
grep -q 'golden-compare.*PASS' docs/agent_plans/20260721-gpu-manual-reliability/results/STATUS_BASELINE.txt
```

DO: Add exact gates for strict nvcc closure builds in 2-D/3-D, preserved CPU golden results, representative multi-box GPU execution, and compute-sanitizer on a named GPU. Compiler diagnostics become port-plan rows. Record unavailable hardware/toolchains as blocked. The current baseline sanitizer failure is environmental (`bin/alamo_gpu-3d-cuda86-g++` absent), not a pass or a source regression.

CHECK: Tier 0 reports scanner, compiler, runtime, sanitizer, baseline efficiency, and profiled-optimization status independently; no missing axis can be summarized as a GPU-native baseline.

### Step 7 - Require a credible GPU-native baseline

VERIFY:

```bash
grep -R '^Class: performance' docs/gpu_manual/patterns/GPU-*.md
```

DO: Define a mandatory `BASELINE_EFFICIENCY` gate separately from `GPU_SAFETY` and `PROFILED_OPTIMIZATION`. It requires device-resident steady-state fields, GPU execution for hot domain loops, static/value dispatch in device call paths, device reductions with bounded host results, no per-tile global synchronization, no hot host quarantine, and no straightforward component-loop launch multiplication. Require one named-GPU timeline/TinyProfiler capture showing that the intended kernels execute on device without avoidable transfer or synchronization domination. No speedup or occupancy threshold is required.

CHECK: a merely launchable but host-staged/serialized implementation fails baseline efficiency; a safe GPU-native implementation may pass without speculative fusion or tensor surgery.

### Step 8 - Reserve “optimized” for measured improvements

VERIFY:

```bash
grep -R '^Class: performance' docs/gpu_manual/patterns/GPU-*.md
```

DO: Keep kernel-specific GPU-015/018/019/020 variants in `PERFORMANCE.md` as hypotheses unless evidence contains named GPU, exact binary/command, dimensions, box layout, kernel time, registers, occupancy, meaningful bandwidth, and preserved golden output. Use “GPU-native baseline” for a safety-plus-baseline-efficiency pass and reserve “optimized” for measured improvements.

CHECK: no unmeasured kernel-specific pattern uses verified/optimized language; later tuning remains an explicit backlog rather than a prerequisite for the base conversion.

## Closeout

- [ ] All v2 schema and recognizer oracles pass
- [ ] Known correctness false negatives are zero
- [ ] Manual/build-error rows are present in PORT_PLAN.csv
- [ ] VALIDATION_MATRIX.csv reports per-pattern, per-axis state
- [ ] BASELINE_EFFICIENCY rejects host fallback, transfer/sync pathologies, runtime device dispatch, and avoidable launch multiplication
- [ ] Tier 0 contains no unqualified project-wide pass
- [ ] results/RESULT.md records metrics, remaining gaps, and hardware-blocked gates
- [ ] results/REVIEW.md contains a fresh adversarial review
- [ ] `results/DONE` exists only when required compiler/runtime gates have evidence or are explicitly scoped to a separate target-port task
- [ ] Session log line appended to `docs/llm/SESSION_LOG.tsv`
