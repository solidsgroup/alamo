# TASK: chamber-fma-adoption
# Folder: docs/agent_plans/20260731-chamber-fma-adoption/

---

## Header

| Field         | Value |
|---------------|-------|
| Risk tier     | 3 |
| Model         | Codex (GPT-5) |
| Verification  | partial-oracle |
| Est. scope    | 2 existing source files via two cherry-picks; task provenance and target-bound validation results |
| Parallel-safe | no; one agent owns the target source state, binaries, and performance adjudication |

## Operating rules

1. Read only the files listed in Context budget. `docs/archive/` is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure, stop and record the discrepancy.
3. Preserve the two reviewed source commits as separate commits; use scoped, imperative messages for adoption-only documentation.
4. Do not merge the campaign branch or import its benchmark/Phase-0 history.
5. Retain MGS32/BF8, 4/4 smoothing, and synchronized MLMG; no rejected configuration may be adopted.
6. Record target commit, source diff, binary/input hashes, exact overrides, raw logs, output hashes, and timing samples.
7. No source edits beyond cherry-picking `bcc1407a2` followed by `d08d6e1ba`.

## Context budget

Read first: `benchmark/status.sh` output, this `PLAN.md`

Read:
- `AGENTS.md`
- `CLAUDE.md`
- `docs/llm/PLAN.md`
- `docs/llm/TASK_TEMPLATE.md`
- `src/BC/BC.H`
- `src/Integrator/BaseField.H`
- `benchmark/ci_golden_compare.sh`
- `benchmark/local_a100_gate.sh`
- `benchmark/build_alamo_local_gpu.sh`
- `benchmark/build_alamo_nova.sh`
- `benchmark/nova_flame_gpu_3d_a2.slurm`
- `benchmark/validate/**`
- `input_copy`
- `/home/jackplum/Projects/alamo-chamber-fma-test/docs/agent_plans/20260730-chamber-fma-transfer/PLAN.md`
- `/home/jackplum/Projects/alamo-chamber-fma-test/docs/agent_plans/20260730-chamber-fma-transfer/results/RESULT.md`
- `/home/jackplum/Projects/alamo-chamber-fma-test/docs/agent_plans/20260730-chamber-fma-transfer/results/REVIEW.md`
- `/home/jackplum/Projects/alamo-chamber-fma-test/docs/agent_plans/20260730-chamber-fma-transfer/results/KEY_EVIDENCE_SHA256SUMS`

Reference only when the named step requires it:
- `/home/jackplum/Projects/alamo-chamber-fma-test/benchmark/chamber_fma_nova_matrix.slurm`
- `benchmark/NOVA_SLURM_RUNBOOK.md`
- `docs/llm/BUG_PATTERNS.md`

Forbidden: `docs/archive/*`, unrelated task folders, propellant sweep data

## Objective

Adopt the reviewed OPT-10 interior physical-BC skip and its mandatory index-type
repair onto the current `chamber-gpu` target without importing the campaign
branch. Bind the adoption to fresh target-branch correctness, sanitizer, and
retained 800-step production evidence before declaring it ready.

## Oracle

Commands:
- `GOLDEN_MODE=cpu benchmark/ci_golden_compare.sh`
- `GOLDEN_MODE=gpu benchmark/ci_golden_compare.sh`
- `TIERS="1 2 3" benchmark/local_a100_gate.sh`
- matched MGS32/BF8, 4/4, synchronized 800-step `input_copy` comparison
- `benchmark/status.sh`

Covers: build integrity, strict golden outputs, device memory/init/race defects,
matched-layout production equivalence, convergence, physics observables, and
target-bound performance.

Does not cover: every production geometry, alternative grid layouts, reduced
smoothing, or native no-sync. The frozen campaign evidence remains the broader
rejection record for those arms.

## Steps

### Step 1 - Freeze target and adoption provenance

VERIFY:
```bash
git status --short --branch
git merge-base --is-ancestor cc22ff9ee HEAD
git diff --quiet cc22ff9ee HEAD -- src/BC/BC.H src/Integrator/BaseField.H
```

DO: Record the clean current target and confirm the two source files are
unchanged from the frozen campaign base. Commit this adoption plan.

CHECK: Clean worktree at a plan-only commit; target base recorded in results.

### Step 2 - Transplant the reviewed source pair

VERIFY: Step 1 passes and both reviewed commits resolve locally.

DO: Cherry-pick `bcc1407a2` and then `d08d6e1ba`, preserving both commits.

CHECK: The target source diff is exactly 13 insertions in `src/BC/BC.H` and
`src/Integrator/BaseField.H`, with patch SHA-256
`6369167a6b823b60760d6a258cdbc71bc6d2dd72a29f96418e2614183f634585`.

### Step 3 - Run local target correctness and sanitizer gates

VERIFY: Clean source state and exact Step 2 diff.

DO: Run strict CPU/GPU golden comparisons and all three local A1000 gate tiers.

CHECK: All commands exit zero; memcheck, initcheck, and racecheck report no
findings; scientific comparisons pass.

### Step 4 - Run target-bound retained production comparison

VERIFY: Step 3 passes; baseline and candidate identities are frozen.

DO: Build matched baseline/candidate binaries from the current target before
and after the source pair. Run `input_copy` for 800 steps at MGS32/BF8, 4/4,
and synchronized MLMG with output enabled, plus interleaved timing repetitions.

CHECK: Scientific output and `thermo.dat` are bitwise identical; solver and
physics gates pass; the candidate shows a reproducible performance benefit.

### Step 5 - Review and close out

VERIFY: All required evidence exists and checksums verify.

DO: Perform a fresh adversarial review, record adjudication in
`results/REVIEW.md`, summarize in `results/RESULT.md`, create `results/DONE`,
and append the session log.

CHECK: `benchmark/status.sh` is green, the worktree is clean, and the adoption
commit history contains only the plan, reviewed source pair, and closeout.

## Checkpoints

- [x] After plan restatement: user accepted the selective-transplant and target-validation recommendation.
- [x] Before source transplant: target is clean and exact source diff is reviewed.
- [ ] Before final retention: strict gates, sanitizer, 800-step correctness, and performance evidence pass.
- [ ] Before closeout commit: fresh adversarial review is adjudicated.

## Adversarial review

After implementation and gates pass, use a fresh reviewer with no campaign
context. Assume the transplant contains a defect; focus on physical-vs-interior
BC semantics, index types, ghost widths, periodic exchange, target-branch
interactions, evidence identity, and whether any test was weakened. Record
findings in `results/REVIEW.md`.

## Closeout

- [ ] Oracle passes; required target gates are green
- [ ] `results/RESULT.md` records changes, evidence, and deviations
- [ ] No changelog/version change unless separately requested
- [ ] `results/DONE`
- [ ] Session log line appended to `docs/llm/SESSION_LOG.tsv`
