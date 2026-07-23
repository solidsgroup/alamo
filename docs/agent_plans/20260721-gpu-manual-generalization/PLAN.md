# TASK: gpu-manual-generalization
# Folder: docs/agent_plans/20260721-gpu-manual-generalization/

---

## Header

| Field         | Value |
|---------------|-------|
| Risk tier     | 0 |
| Model         | opus for contracts and status judgment; cheap agents for bounded docs/scripts |
| Verification  | partial-oracle |
| Est. scope    | GPU manual docs, templates, CSV schemas, and recognizer scripts; no source changes |
| Parallel-safe | yes: Tier 1 halves, recognizer tooling, and framework documents have disjoint ownership |

## Operating rules

1. Read only the files listed in Context budget. `docs/archive/` is forbidden.
2. Preserve user-owned changes in `BRIEF.md`, `BUILD_LOG.md`, `docs/llm/PLAN.md`, the task plan preceding this one, and `docs/gpu_manual.zip`.
3. All product writes stay under `docs/gpu_manual/**`; this task folder records verification and review. No `src/`, branch, worktree, or numerical changes.
4. The supplied generalization plan supersedes held-out-corpus and formal precision/recall work. Retain separate transform, recognizer, port, runtime, and efficiency states without team-scale metrics.
5. Treat scanner output as advisory. Compiler diagnostics and dispositioned inspection always outrank regex output.
6. Do not invent a pilot port. Templates may show clearly labeled illustrative rows, while pilot-dependent acceptance remains pending until a separately authorized port instantiates it.
7. No commits in this shared dirty worktree; report the complete diff and verification instead.

## Context budget

Read first: `benchmark/status.sh` output, this `PLAN.md`, the user's GPU Manual Generalization Plan
Read: `docs/gpu_manual/BRIEF.md`, `INDEX.md`, all `patterns/GPU-*.md`, `recognizers/table.csv`, `recognizers/scan.py`, `build/validate_patterns.py`, `FEATURES.md`, `ONE_OFFS.md`, `BUILD_LOG.md`, `evidence/*.md`, `build/phase6/*`
Reference only for corpus evidence: `build/HUNK_MAP.csv`, exact source paths already named by Tier 1 evidence, official AMReX/CUDA documentation
Forbidden: `docs/archive/*`, source edits, new integrator closures/oracle decks, optimization implementation

## Objective

Generalize the GPU manual from a Flame-derived corpus into a reusable single-maintainer system. Readers must be able to distinguish invariant device semantics, port-supplied contracts, and corpus examples; a new port must be able to start from templates for scope, compiler-first closure, inspection, validation, and efficiency without `chamber-gpu` access. The manual must also define a stateful advisory-recognizer loop and a required per-port harvest/status process.

## Oracle

Commands:

```bash
python3 docs/gpu_manual/build/validate_patterns.py
python3 -m unittest discover -s docs/gpu_manual/recognizers -p 'test_*.py'
python3 docs/gpu_manual/recognizers/scan.py --root docs/gpu_manual/recognizers --table docs/gpu_manual/recognizers/table.csv --out /tmp/gpu-manual-coverage.csv
```

Covers: generalized Tier 1 schema, pattern inventory/status/evidence fields, versioned recognizer and coverage schemas, scanner state preservation, template/schema presence, Tier 0 token budget, and policy references.

Does not cover: a new port's analytic oracle, physics-specific tolerances, real compiler closure, real-GPU runtime/sanitizer evidence, or named-GPU profiling. Those remain `pending-pilot`/`blocked`, never pass.

## Steps

### Step 1 - Define contracts and authority

VERIFY: repository is `manual-build`; baseline status is recorded; no required source write exists.

DO: define the invariant/port/corpus legend, primary-source evidence map, validation contract, onboarding schemas/templates, architecture policies, efficiency contract, status/harvest model, and blind-spot register.

CHECK: every required artifact is linked from Tier 0 or the brief and can be instantiated without a source branch.

### Step 2 - Migrate Tier 1

VERIFY: exactly 25 active patterns exist; only GPU-007 and GPU-016 have prior closed-book verification.

DO: migrate each pattern to the generalized schema, remove port-specific commands from correctness `Verify`, mark before/after material as corpus examples, add primary evidence for invariant claims, and point architectural questions at their single policy home.

CHECK: validator passes; no correctness pattern asserts a Flame-specific constraint or verification command as universal procedure.

### Step 3 - Make the scanner a maintained advisory instrument

VERIFY: v1 table and coverage schemas are captured in `BUILD_LOG.md`.

DO: version the table/coverage schemas, emit line-addressable candidate/converted states, preserve reviewed `not-applicable`/`false-positive` dispositions on rerun, provide the feedback procedure, and add focused scanner tests.

CHECK: candidate, converted, mixed-site, false-positive, and not-applicable cases converge; no regex result is described as proof of completion.

### Step 4 - Integrate the brief and Tier 0

VERIFY: contracts, patterns, and scanner agree on status names and precedence.

DO: rewrite the builder/maintainer brief around generic port phases; update INDEX within budget; formalize cross-family status, closed-book onboarding gate, per-port harvest, and pilot dependencies; retain the three-tier layout.

CHECK: a reader cannot mistake Flame evidence for universal procedure; no unqualified project-wide pass remains.

### Step 5 - Validate and review

VERIFY: scoped oracle commands pass.

DO: run documentation consistency searches, inspect the diff, obtain a fresh adversarial review, record results and deviations, and append the session log.

CHECK: findings are fixed or explicitly listed as pilot-blocked; create `results/DONE` only for the documentation implementation, not for pilot-dependent runtime acceptance.

## Adversarial review

A fresh reviewer checks the completed diff for instance-as-invariant leakage, conflicting policy homes, scanner overclaiming, schema drift, hidden Flame commands in correctness verification, missing retirement/harvest obligations, and misleading pass language.

## Closeout

- [x] Scoped oracle passes
- [x] All 25 patterns use the generalized schema
- [x] Scope, closure, inspection, validation, efficiency, and harvest templates exist
- [x] Stateful recognizer lifecycle is versioned and tested
- [x] Tier 0 is within budget and labels invariant/port/corpus content
- [x] Pilot-dependent gates are explicit, not inferred
- [x] `results/RESULT.md` and `results/REVIEW.md` record evidence and open dependencies
- [x] `results/DONE` marks the docs/tooling task only
- [x] Session log line appended
