# TASK: thermoelastic-chain-rule
# Folder: docs/agent_plans/20260730-thermoelastic-chain-rule/

---

## Header

| Field         | Value |
|---------------|-------|
| Risk tier     | 3 (constitutive stress and tangent correction) |
| Model         | GPT-5 |
| Verification  | partial-oracle (analytic finite differences plus repository gates) |
| Est. scope    | 4 source/test/build files plus this task record |
| Parallel-safe | yes for independent review; builds remain sequential |

## Operating rules

1. Work only in the clean `chamber-gpu` worktree created for this repair.
2. Preserve unrelated branches and user-owned changes. Stage only files listed
   in this plan.
3. Do not weaken tolerances or golden checks to make a gate pass.
4. Use one commit after the complete oracle and adversarial review pass.
5. Commit message: `solid: fix predeformation chain rule
   (20260730-thermoelastic-chain-rule)`.

## Context budget

Read:
- this PLAN.md
- `docs/llm/PLAN.md`
- `src/Model/Solid/Finite/NeoHookeanPredeformed.H`
- `src/Model/Solid/Finite/PseudoAffine/Cubic.H` as the chain-rule analogue
- `src/Test/Model/Solid/Finite/NeoHookeanPredeformed.H`
- `src/test.cc`
- gate/build scripts needed to run the oracle

Do not inspect `docs/archive/`.

## Objective

Apply the accepted thermoelastic repair to `chamber-gpu`:

- transform `DW` and `DDW` through both required factors of `F0^-1`;
- make `Random()` return the predeformed model with a nonsingular perturbation
  of identity;
- export all nine `F0` components in 3-D;
- add deterministic 2-D/3-D derivative, free-expansion, and field-layout
  regression tests.
- keep the Flame CUDA source closure link-complete after the branch's
  concurrent development re-merge added `InputScraper` and `OutputLog`
  references.

## Oracle

Commands from the isolated worktree root:

- build and run the full test executable in 2-D;
- build and run the full test executable in 3-D;
- `bash benchmark/status.sh`;
- build the expected local CUDA 3-D binary if absent;
- `TIERS="1 2" bash benchmark/local_a100_gate.sh`;
- `git diff --check`.

The focused derivative test requires relative errors below `1e-6` for `DW`
and `1e-4` for `DDW`, in both dimensions. It also requires stress-free
expansion and the complete dimension-specific field layout.

Covers: analytic chain-rule consistency, host 2-D/3-D compilation and tests,
CPU golden/smoke behavior, device-pattern lint, strict GPU golden behavior,
and GPU memory safety.

Does not cover: a new production-scale chamber solve or a casing-resolved
refinement study.

## Steps

### Step 0 - Baseline and isolation

VERIFY: worktree is clean on local `chamber-gpu`.

DO: run `benchmark/status.sh` before edits and record environmental failures.

CHECK: baseline source is the unfixed implementation; unrelated worktrees are
untouched.

### Step 1 - Constitutive and export repair

DO: apply the accepted `DW`/`DDW` chain rule, shared inverse helper,
`Random()` return/type correction, and complete 3-D field names and copies.

CHECK: line-by-line diff contains no unrelated source changes.

### Step 2 - Deterministic regression coverage

DO: add the focused test header and register its derivative, free-expansion,
and field-layout tests in `src/test.cc`.

CHECK: tests compile and pass in both 2-D and 3-D.

### Step 3 - Repository correctness gates

DO: run every oracle command, building the missing local CUDA binary first if
needed. If the branch's GPU source closure omits a directly referenced
translation unit, add only that translation unit to the existing closure.

CHECK: all code-relevant gates pass. A baseline environment failure must be
resolved before commit.

### Step 4 - Adversarial review and closeout

DO: obtain a fresh reviewer pass over the exact diff, emphasizing tensor
indices, major symmetry, GPU compatibility, dimension coverage, and test
adequacy. Review the diff line by line in the main session.

CHECK: no material issue remains. Write `results/RESULT.md`,
`results/REVIEW.md`, `results/DONE`, and append `docs/llm/SESSION_LOG.tsv`.

## Checkpoints

The user's explicit acceptance of the verified repair and request to apply,
test, and commit authorizes all implementation checkpoints in this scoped
plan. Commit remains conditional on successful tests and review.

## Adversarial review

Use a fresh reviewer with no edit authority. Required questions:

- Are the indices in `DW` and `DDW` the exact chain rule for
  `W(F F0^-1)`?
- Does storage as `Sym::Major` remain valid and correctly assigned?
- Is every edited path valid in both 2-D and 3-D host/device builds?
- Do deterministic tests fail on the original bug and cover the export fix?
- Was any unrelated behavior or gate changed?

## Closeout

- [x] 2-D and 3-D tests pass
- [x] repository status/golden gate passes
- [x] strict GPU build, valid golden cases, smoke, and memory gates pass;
      the documented stale `rod_and_tube` reference is adjudicated
- [x] adversarial review passes
- [x] `results/RESULT.md` and `results/REVIEW.md`
- [x] `results/DONE`
- [x] `docs/llm/SESSION_LOG.tsv` line
- [x] one scoped commit (this task's commit)
