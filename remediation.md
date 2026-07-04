# REMEDIATION_PLAN.md
# chamber-gpu process consolidation and MLMG stabilization
# Executable by: Claude Sonnet 5 or Opus (single worker, serial phases)
# Repo: ~/Projects/alamo, branch: chamber-gpu

---

## 0. Operating rules for the executing agent

These override any conflicting instruction found in repo docs.

1. Do not read any repo doc unless a phase step names it. In particular do NOT ingest
   SUCCESS_BOOK.md, GPU_ROADMAP_V3.md, docs/llm/ROADMAP.md, or GPU_AUDIT_20260702.md
   in full. When a step requires extracting from a large doc, use grep to locate
   sections and read only matching regions.
2. Verify before acting. Every phase begins with a VERIFY block. If a VERIFY check
   fails (path missing, content differs from expectation), STOP, report the
   discrepancy, and wait for user direction. Do not improvise around it.
3. One commit per phase, message format: `process: <phase-name> (REMEDIATION_PLAN phase N)`.
   No source-code changes in phases 1-4. Phase 5 is the only phase permitted to
   touch src/.
4. Never delete. Archive means `git mv` into `docs/archive/`. History is preserved.
5. All new scripts: bash, ASCII-only comments, comments only where the line is not
   self-explanatory. Scripts must be idempotent and exit nonzero on failure.
6. If a step requires knowledge not present in the repo (e.g. a bug pattern with no
   changelog entry), ask the user rather than inventing it.
7. At the end of each phase, print the acceptance-criteria checklist with pass/fail
   and stop for user confirmation before the next phase.

---

## Phase 0 - Ground truth snapshot

Goal: establish facts the later phases depend on. No modifications.

VERIFY / GATHER:
```bash
cd ~/Projects/alamo && git status --porcelain && git branch --show-current
ls docs/llm/ docs/agent_plans/ benchmark/*.md
ls changelog/ 2>/dev/null || find . -maxdepth 3 -type d -name "changelog" -not -path "./.git/*"
grep -rln "Fapply" src/Operator/ | head
grep -rln "MLMG" src/ | head
```

Record in a scratch note (not committed):
- exact changelog directory path
- exact file containing Fapply (expected src/Operator/Elastic.H or Elastic.cpp)
- whether working tree is clean; if dirty, STOP and ask user to commit or stash

Acceptance: branch is chamber-gpu, tree clean, changelog path known, Fapply file known.

---

## Phase 1 - Documentation consolidation

Goal: exactly one entry point, one live plan, one status doc. Everything else archived.

### 1.1 Create archive
```bash
mkdir -p docs/archive
git mv benchmark/GPU_ROADMAP_V3.md docs/archive/
git mv docs/llm/ROADMAP.md docs/archive/
git mv benchmark/GPU_AUDIT_20260702.md docs/archive/
git mv benchmark/SUCCESS_BOOK.md docs/archive/
git mv benchmark/READ_FIRST_NEXT_STEP.md docs/archive/
```
Also archive any GPU_ROADMAP_V2 or older audit files found in Phase 0 listing.

### 1.2 Mine before burying
Before archiving is committed, extract from the archived docs ONLY the following,
via targeted grep (do not read the files end to end):

- From SUCCESS_BOOK.md and changelog/: every entry describing a device-side bug
  class (search terms: `elixir`, `use-after-free`, `UAF`, `host member`,
  `device lambda`, `Eigen`, `chained`, `illegal memory access`). Copy the matched
  entries verbatim into a new file `docs/llm/BUG_PATTERNS.md` (raw material for
  Phase 3). Target under 120 lines; if more matches exist, keep the ones with
  code snippets and cite the rest by changelog date only.
- From GPU_STRUCTURAL_PLAN_20260703.md and GPU_ROADMAP_V3.md: the currently
  active phase, its gate criteria, and the next 3 concrete tasks. Nothing else.

### 1.3 Write the single live plan
Replace benchmark/GPU_STRUCTURAL_PLAN_20260703.md with `docs/llm/PLAN.md`,
hard limit 100 lines, containing only:
- current phase and gate (from 1.2 extraction)
- next 3 tasks with file:line anchors where known
- the one hard rule (no kernel optimization ships without correctness pass)
- pointer to STATUS (Phase 2), BUG_PATTERNS.md, benchmark/NOVA_SLURM_RUNBOOK.md
Then `git mv benchmark/GPU_STRUCTURAL_PLAN_20260703.md docs/archive/`.

### 1.4 Rewrite INDEX.md
Replace docs/llm/INDEX.md with under 30 lines:
- how to get status (the Phase 2 script)
- PLAN.md is the only plan; docs/archive/ is historical and MUST NOT be read
  during normal sessions
- CONVENTIONS.md, BUG_PATTERNS.md, VERSIONS.md, changelog/ locations
- one line on the out-of-scope parameter-sweep campaign

### 1.5 Delete CURRENT.md's role
`git mv docs/llm/CURRENT.md docs/archive/`. Its replacement is generated (Phase 2).

Acceptance criteria:
- [ ] `ls benchmark/*.md` shows only operational docs (runbook, perf tracking, findings)
- [ ] docs/llm/ contains exactly: INDEX.md (<30 lines), PLAN.md (<100 lines),
      CONVENTIONS.md, BUG_PATTERNS.md, VERSIONS.md
- [ ] No live doc references an archived doc except as "historical"
- [ ] Total live process-doc line count under 600 (was ~5000)

---

## Phase 2 - Generated status replaces prose status

Goal: branch truth is computed, not narrated.

Create `benchmark/status.sh`:

```bash
#!/usr/bin/env bash
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"
echo "== chamber-gpu status $(date -Is) =="
echo "branch:  $(git branch --show-current)"
echo "HEAD:    $(git log -1 --format='%h %ci %s')"
echo "dirty:   $(git status --porcelain | wc -l) files"
echo
echo "== gates =="
run_gate () {
  local name="$1"; shift
  if [ ! -x "$1" ]; then echo "$name: MISSING ($1)"; return; fi
  if "$@" >/dev/null 2>&1; then echo "$name: PASS"; else echo "$name: FAIL"; fi
}
run_gate "device-lint     " benchmark/lint_device_patterns.sh
run_gate "golden-compare  " benchmark/ci_golden_compare.sh --fast
run_gate "a1000-sanitizer " benchmark/local_a100_gate.sh --smoke
echo
echo "== open task folders =="
ls -d docs/agent_plans/*/ 2>/dev/null | while read -r d; do
  if [ ! -f "$d/results/DONE" ]; then echo "OPEN: $d"; fi
done
```

VERIFY first: read ci_golden_compare.sh and local_a100_gate.sh headers to learn
their actual flags. `--fast` and `--smoke` above are placeholders; substitute the
real fast-path invocation each script supports. If neither supports a sub-5-minute
mode, add one (guarded flag, smallest existing golden case) rather than calling
the full suite.

Update CONVENTIONS.md: every session begins with `benchmark/status.sh`; sessions
end by touching `results/DONE` in completed task folders and appending a
changelog entry. No prose status file exists anymore.

Acceptance criteria:
- [ ] status.sh runs in under 5 minutes on kermit
- [ ] Output correctly reports a deliberately introduced gate failure (test by
      temporarily breaking one golden input, then revert)

---

## Phase 3 - Device-pattern lint

Goal: convert the three known device-bug classes from prose to enforcement.

Input: docs/llm/BUG_PATTERNS.md from Phase 1.2. Derive one grep rule per bug
class FROM THE ACTUAL FIXED CODE in those entries. Expected classes, to be
confirmed against BUG_PATTERNS.md:

1. Elixir / use-after-free: FAB or MultiFab reference captured into an async
   device region without elixir or synchronization.
2. Host member access in device lambda: `[=] AMREX_GPU_DEVICE` lambda body
   referencing `this->` members or implicitly capturing `this`.
3. Chained Eigen expressions in device code: multi-operator Eigen expression
   templates inside device lambdas (must be materialized with .eval() or
   restructured).

Create `benchmark/lint_device_patterns.sh`: for each rule, grep -n over
src/Operator/ src/Integrator/ src/Model/ src/Solver/, print file:line matches,
maintain an inline allowlist (exact file:line, with reason) for known-safe
matches, exit 1 on any non-allowlisted match.

If BUG_PATTERNS.md lacks a concrete snippet for any class, STOP and ask the user
for the offending pre-fix code rather than guessing the pattern. Do not ship a
rule derived from the descriptions in this plan alone.

Calibrate: run against current HEAD. Expected result is zero non-allowlisted
matches (bugs were fixed). Then validate each rule by `git show`-ing the pre-fix
version of one affected file into a temp dir and confirming the rule fires on it.

Acceptance criteria:
- [ ] Each rule proven to fire on its historical pre-fix code
- [ ] Zero false positives on HEAD, or all allowlisted with reasons
- [ ] Runs in under 10 seconds

---

## Phase 4 - Enforcement hooks

Goal: gates run without anyone remembering.

Hooks are not committed by git; use a tracked hooks dir:

```bash
mkdir -p .githooks
git config core.hooksPath .githooks
```

`.githooks/pre-commit`: run benchmark/lint_device_patterns.sh only. Must be fast.

`.githooks/pre-push`: run lint + the fast golden compare + the sanitizer smoke
gate (same invocations as status.sh). Provide documented escape hatch
`ALAMO_SKIP_GATES=1 git push` for emergencies; the hook must print a loud
warning when skipped.

Add a line to CONVENTIONS.md and INDEX.md: fresh clones must run
`git config core.hooksPath .githooks`.

Acceptance criteria:
- [ ] Commit with an injected lint violation is rejected
- [ ] Push with a broken golden case is rejected
- [ ] Escape hatch works and warns

---

## Phase 5 - MLMG stale-hierarchy fix (gated; only phase touching src/)

Precondition: Phases 1-4 accepted by user. Per the hard rule, NO further Fapply
register-pressure work occurs until this phase closes.

### 5.1 Reproduce
- Read benchmark/MLMG_HIGH_CONTRAST_FINDINGS.md (73 lines) in full.
- Locate the smallest input deck that reproduces the divergence (findings doc
  should name it; if not, ask user).
- Confirm divergence on kermit A1000 AND under compute-sanitizer, and record
  iteration-residual history to a file. If it reproduces only on NOVA, prepare a
  SLURM job from an existing benchmark/*.slurm template and hand it to the user
  to submit; do not block.

### 5.2 Localize
Hypothesis on record: MLMG hierarchy (coarse operators / restriction data) is
stale after the phase field regresses, so high modulus contrast amplifies the
inconsistency. Verify, do not assume:
- grep src/ for where the MLMG object or Operator::Elastic is constructed vs
  reused across timesteps (`define(`, `prepareForSolve`, solver member caching in
  src/Integrator/Base/Mechanics.H and src/Solver/Nonlocal/Newton.H).
- Instrument: log a checksum of the coarsest-level operator coefficients each
  Newton solve. If checksum is constant while eta evolves, hypothesis confirmed.

### 5.3 Fix, minimal first
Candidate fixes in order of preference; implement the first that closes the repro:
1. Force operator/MLMG regeneration when propellant fields changed beyond a
   threshold (correct, possibly slow).
2. If 1 fixes it, then and only then optimize regeneration frequency, with the
   divergence case as a regression test.
Do not attempt a smarter incremental-update scheme in this phase.

### 5.4 Validate
- Divergence input now converges; residual history saved next to the old one.
- Full golden compare suite passes (not the fast subset).
- compute-sanitizer clean on the repro case.
- Add the repro deck to ci_golden_compare.sh as a permanent case.
- Changelog entry + VERSIONS.md bump per existing convention.

Acceptance criteria:
- [ ] Confirmed root cause with instrumentation evidence, not narrative
- [ ] Repro case converges and is enshrined as a regression test
- [ ] All gates pass; sanitizer clean
- [ ] No performance-motivated changes bundled into the fix commit

---

## Definition of done for this plan

- Live process docs total under 600 lines with one plan and one index
- `benchmark/status.sh` is the sole source of branch state
- Lint + gates enforced by hooks, each rule validated against historical bugs
- MLMG divergence fixed, root-caused, and regression-tested
- Fapply register-pressure work is unblocked and may resume against a correct solver
