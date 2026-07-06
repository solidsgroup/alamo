# TASK: meta-workflow-install
# Folder: docs/agent_plans/YYYYMMDD-meta-workflow-install/

## Header

| Field        | Value                                   |
|--------------|-----------------------------------------|
| Risk tier    | 0                                       |
| Model        | sonnet                                  |
| Verification | full-oracle (each step has a CHECK)     |
| Est. scope   | 3 new files, 1 edit, ~60 lines          |
| Parallel-safe| no (touches shared process files)       |

## Operating rules

1. Read only files named below. docs/archive/ forbidden.
2. VERIFY must pass before each DO. On failure: STOP and report.
3. One commit per step: `process: <what> (meta-workflow-install)`.
4. No scope expansion.

## Context budget

Read: this PLAN.md, .githooks/pre-commit, docs/llm/INDEX.md
Provided by user alongside this plan: CLAUDE.md, TASK_TEMPLATE.md

## Objective

Install the session entry point (CLAUDE.md), the reusable task template, the
session telemetry log, and a doc-growth guard in the pre-commit hook.

## Steps

### Step 1 - Install CLAUDE.md and TASK_TEMPLATE.md
VERIFY:
```bash
test ! -f CLAUDE.md && test ! -f docs/llm/TASK_TEMPLATE.md && \
test -f docs/llm/PLAN.md && test -f benchmark/status.sh
```
If PLAN.md or status.sh are missing or named differently, STOP and report the
actual names; CLAUDE.md references must match reality before install.
DO: copy user-provided CLAUDE.md to repo root and TASK_TEMPLATE.md to
docs/llm/. Fix any path in CLAUDE.md that VERIFY showed to differ.
CHECK:
```bash
test -f CLAUDE.md && test -f docs/llm/TASK_TEMPLATE.md && \
grep -q status.sh CLAUDE.md
```

### Step 2 - Create SESSION_LOG.tsv
VERIFY: `test ! -f docs/llm/SESSION_LOG.tsv`
DO: create docs/llm/SESSION_LOG.tsv with single header line:
```
date	task	model	tier	cost_usd	outcome
```
Tab-separated. outcome is merged|abandoned|partial.
CHECK: `head -1 docs/llm/SESSION_LOG.tsv | awk -F'\t' '{print NF}'` prints 6.

### Step 3 - Doc-budget guard in pre-commit
VERIFY: `test -x .githooks/pre-commit` and read it; confirm it currently runs
the device lint and note its structure.
DO: append to .githooks/pre-commit (do not modify existing checks):
```bash
DOC_BUDGET=600
live_docs=$(find docs/llm -name '*.md' ! -name 'BUG_PATTERNS.md')
total=$(cat $live_docs | wc -l)
if [ "$total" -gt "$DOC_BUDGET" ]; then
  echo "pre-commit: live process docs ${total} lines > budget ${DOC_BUDGET}" >&2
  exit 1
fi
new_md=$(git diff --cached --name-only --diff-filter=A | grep '\.md$' | \
  grep -vE '^(docs/archive/|docs/agent_plans/|changelog/)' || true)
if [ -n "$new_md" ]; then
  echo "pre-commit: new .md outside allowed dirs:" >&2
  echo "$new_md" >&2
  echo "add to archive/agent_plans/changelog, or edit allowlist knowingly" >&2
  exit 1
fi
```
CHECK:
```bash
touch STRAY_NOTES.md && git add STRAY_NOTES.md
git commit -m tmp 2>&1 | grep -q "outside allowed dirs" && \
git reset STRAY_NOTES.md && rm STRAY_NOTES.md && echo GUARD_OK
```
Must print GUARD_OK. Then confirm a normal commit of Steps 1-2 files passes.

### Step 4 - Point INDEX.md at the new pieces
VERIFY: read docs/llm/INDEX.md, confirm under 30 lines.
DO: add (or amend) lines referencing CLAUDE.md as the session entry point,
TASK_TEMPLATE.md as mandatory for new tasks, SESSION_LOG.tsv location. Keep
INDEX.md under 30 lines total; remove redundancy rather than exceeding it.
CHECK: `wc -l docs/llm/INDEX.md` under 30; grep finds all three references.

## Closeout

- [ ] status.sh all green
- [ ] results/RESULT.md written; touch results/DONE
- [ ] First real line appended to SESSION_LOG.tsv for this session itself
- [ ] changelog/ entry noting workflow-install (no VERSIONS.md bump; process only)
