# RESULT: meta-workflow-install

Date: 2026-07-06 | Model: opus-4.8 | Risk tier: 0 (process only, no src/ changes)

## Outcome: installed & verified — DONE held on one unrelated gate

All four plan steps are complete and each CHECK passed. `results/DONE` is
**not** touched because closeout item #1 ("status.sh all green") is not
satisfiable: `golden-compare` is a pre-existing FAIL unrelated to this task
(see below). Flip to DONE once that gate is green or explicitly waived.

## Steps

| Step | What | State |
|------|------|-------|
| 1 | Install CLAUDE.md + docs/llm/TASK_TEMPLATE.md | already committed `ec8c47b06` |
| 2 | Create docs/llm/SESSION_LOG.tsv (6-col header) | already committed `c8e44736c` |
| 3 | Doc-budget guard in .githooks/pre-commit | committed `18f758adb` |
| 4 | Point INDEX.md at CLAUDE.md/TASK_TEMPLATE/SESSION_LOG | committed `42896dc2e` |

## Deviations from the plan as written (both necessary; reported per rule 2)

1. **Doc-budget scope.** The plan's guard counted `find docs/llm -name '*.md'`
   (recursive) against `DOC_BUDGET=600`, but the repo already holds 1534 such
   lines — dominated by append-only `changelog/` (~578) and the historical
   `perf/` report (623). Installed verbatim it would `exit 1` on *every* commit,
   including its own install commit. Per user decision, the find now excludes
   `docs/llm/changelog/*` and `docs/llm/perf/*` (append-only / historical),
   leaving 313/600 with headroom. `changelog/` under a shrinking budget also
   contradicts CLAUDE.md's append-only rule.

2. **Allowlist path fix.** The guard's `new_md` allowlist permitted `changelog/`
   but the real path is `docs/llm/changelog/`; as written it would have *blocked*
   this task's own closeout changelog entry. Corrected to `docs/llm/changelog/`
   (same reality-correction the plan sanctions for paths in Step 1).

3. **Insertion point.** The plan said "append" to pre-commit, but the hook ended
   in `exit $?` (a literal append is dead code). The guard is inserted after the
   device-lint runs, preserving the lint's exit status on failure. No existing
   check was altered.

## Verification

- Step 3 CHECK: stray `STRAY_NOTES.md` blocked with "outside allowed dirs" → GUARD_OK.
- Normal commits (the hook change, the INDEX.md change) pass the guard: device-lint
  0 violations + doc-budget 313/600 + no stray new .md.
- Step 4 CHECK: INDEX.md = 29 lines (<30); grep finds CLAUDE.md, TASK_TEMPLATE.md,
  SESSION_LOG.tsv.

## Pre-existing failure (out of scope, NOT introduced here)

`benchmark/status.sh` reports `golden-compare: FAIL`. Root cause:
`rod_and_tube_step2/cpu failed with exit 6` in `benchmark/baseline_suite.py`
(the known rod_and_tube stability deck; stale CPU binary / Mode-D recipe). This
task touched only `.githooks/pre-commit` and `docs/llm/INDEX.md` — neither feeds
that CPU sim. `device-lint: PASS`, `a100-sanitizer: PASS`.

## Untouched strays noted (not part of this task)

Root-level untracked files left as-is: `claude1.md`, `TASK_TEMPLATE.md` (leftover
root copy; the installed one is `docs/llm/TASK_TEMPLATE.md`), `chamber_gpu_changes.diff`.
The new guard will flag `claude1.md` / root `TASK_TEMPLATE.md` if ever staged.

## Post-close note (2026-07-06, efficiency-campaign)

The withheld gate is now green: golden-compare FAIL was a stale-binary
selection bug in ci_golden_compare.sh (lexicographic sort picked Jun-29
alamo-2d-perf-clang++ over the just-built alamo-2d-g++), fixed in
docs/agent_plans/20260706-efficiency-campaign/. DONE granted.
