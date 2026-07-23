# TASK: gpu-manual-hostile-review
# Folder: docs/agent_plans/20260721-gpu-manual-hostile-review/

---

## Header

| Field         | Value |
|---------------|-------|
| Risk tier     | 0 |
| Model         | top-tier manager for pattern/contract judgment; grunts for source reconnaissance and fixtures |
| Verification  | partial-oracle |
| Est. scope    | GPU manual docs, recognizer rules/tests, schemas, and task evidence; no source conversion |
| Parallel-safe | yes: recognizer/source audit and performance/source audit are read-only and disjoint; implementation is integrated centrally |

## Operating rules

1. Preserve all existing user/manual changes and `docs/gpu_manual.zip`; do not edit source.
2. Treat the hostile review as hypotheses until reproduced against the current generalized manual and current source.
3. Prefer high-recall shape candidates plus explicit confirmation over identifier lists. Where regex cannot establish reachability, receiver type, or control flow, say so and require compiler/AST/inspection evidence.
4. Coverage is generated per port and source revision. A shipped snapshot may be retained only as task evidence, never as the port plan.
5. Do not claim Hydro or Fracture is ported, safe, or performant. Source snippets are recognizer/contract fixtures only.
6. No optimization strategy or numerical change is authorized. Architecture evidence may fail or block a port without prescribing a physics implementation.
7. No commits in the shared dirty worktree; report verification and remaining pilot dependencies.

## Context budget

Read first: current `benchmark/status.sh`, this plan, hostile review
Read: current `docs/gpu_manual/**` excluding archived build diffs; `src/Integrator/Hydro.{cpp,H}`, `src/Integrator/Fracture.H`, `src/Model/Gas/**`, `src/Model/Interface/Crack/**`, `src/Solver/Local/Riemann/**`
Reference: official AMReX, CUDA, and Compute Sanitizer documentation; current task results
Forbidden: `docs/archive/*`, source edits, integrator-specific conversion plans, numerical/optimization implementation

## Objective

Repair the current manual's live transferability defects exposed by Hydro and Fracture. Replace identifier-fitted recognition with high-recall shapes and honest manual/compiler boundaries, generalize value/static dispatch, add a host-only numerical-kernel pattern, make coverage per-port and revision-bound, and add mandatory GPU-native architecture evidence before baseline/performance claims.

## Oracle

```bash
python3 docs/gpu_manual/build/validate_patterns.py
python3 -m unittest discover -s docs/gpu_manual/recognizers -p 'test_*.py'
python3 docs/gpu_manual/recognizers/scan.py --root . --table docs/gpu_manual/recognizers/table.csv --port-id hostile-current --source-revision "$(git rev-parse HEAD)" --out /tmp/hostile-current.csv
python3 docs/agent_plans/20260721-gpu-manual-hostile-review/verify_current_coverage.py /tmp/hostile-current.csv
cmp /tmp/hostile-current.csv docs/agent_plans/20260721-gpu-manual-hostile-review/results/current-source-coverage.csv
git diff --check
```

Covers: pattern/schema consistency, absence of banned identifier-fitted candidates, Hydro/Fracture-shaped recognizer fixtures, per-port coverage metadata, candidate precedence/disposition behavior, architecture/template completeness, status semantics, and Tier 0 budget.

Does not cover: actual device-call closure, correct tuple/value design for Hydro or Fracture, numerical validation, unseen-target transform success, hardware resource results, or speedup. Those remain open port work.

## Steps

### Step 1 - Reproduce and classify the review

VERIFY: current status and source anchors are recorded.

DO: separate stale claims from live defects; record exact Hydro/Fracture/Gas/Riemann evidence.

CHECK: no manual change is justified only by the old review's quoted text.

### Step 2 - Generalize detection and value dispatch

VERIFY: current table still contains the named identifiers and current scanner misses the reproduced blockers.

DO: use broad declaration/call/write/expression/control-flow shapes, explicit confirmation, and manual/compiler detection where lexical matching is unsound. Generalize GPU-002 and add a self-contained value/tuple-dispatch example.

CHECK: production-table fixtures flag host/container state, indirect device calls, aggregate writes, repeated expressions, branch-only work, host sentinels, and uninitialized aggregates without relying on Flame/Elastic identifiers.

### Step 3 - Add host-only numerical-kernel and native-shape contracts

VERIFY: Hydro Riemann and Fracture crack paths demonstrate the missing class; existing baseline lacks layout/kernel-graph/resource policy.

DO: add GPU-031 plus a mandatory GPU-native shape contract and templates for field layout, kernel graph, residency/transfers, call-chain complexity, compiler resources, divergence, shared memory, atomics, and arenas.

CHECK: Hydro/Fracture-shaped fixtures have required inspection questions; the contract has evidence-based pass/fail without universal performance thresholds.

### Step 4 - Make coverage and verification non-misleading

VERIFY: the root static coverage artifact and file-repaired gate status are visible.

DO: move coverage to per-port evidence with required port/revision metadata; append the historical GPU-002 correction; change status to draft/file-verified/transfer-verified/cross-family; require frozen unseen targets and prohibit repair/retest promotion.

CHECK: no root `COVERAGE.csv` is an operational input; GPU-007/GPU-016 are file-verified only; current-source coverage is task evidence.

### Step 5 - Validate and review

VERIFY: scoped oracle passes.

DO: run current-source scan, quantify Hydro/Fracture findings, obtain a focused reviewer audit, fix findings, and record result/session evidence.

CHECK: live blocker items 1-5 are implemented in the manual while stale review claims are not reintroduced as fact.

## Closeout

- [x] Shape recognizers replace named-instance candidate rules
- [x] GPU-002 covers polymorphic and container-backed state with worked example
- [x] GPU-031 covers host-only numerical call chains
- [x] Per-port coverage schema/CLI is revision-bound and tested
- [x] GPU-native shape contract and templates cover all architecture topics
- [x] File-specific versus unseen/cross-family status is honest
- [x] Scoped oracle and focused reviewer pass
- [x] Results, DONE marker, and session log recorded
