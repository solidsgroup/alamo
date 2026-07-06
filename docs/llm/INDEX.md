# docs/llm — chamber-gpu map (read this first)

Session entry point is the repo-root `CLAUDE.md` (mandatory session checklist);
read it, then this file and `CONVENTIONS.md`, at session start.

## Status

Run `benchmark/status.sh` — sole source of branch truth (build gates, golden
compare, sanitizer smoke, open task folders). No prose status file exists.
Fresh clones: run `git config core.hooksPath .githooks` once (gates in `.githooks/`).

## Plan

`PLAN.md` is the only live plan (current phase, gate, next 3 tasks). Under
100 lines by design — prune it when it grows, don't append. Every new task
copies `TASK_TEMPLATE.md` into `docs/agent_plans/YYYYMMDD-<name>/PLAN.md`.

## References

- `CONVENTIONS.md` — session start/end checklist.
- `BUG_PATTERNS.md` — device bug classes already fixed once, enforced by
  `benchmark/lint_device_patterns.sh`.
- `VERSIONS.md` — semver ledger. `changelog/` — one entry per closed work item.
- `SESSION_LOG.tsv` — one tab-separated line per session (date/task/model/tier/cost/outcome).

## Historical & out of scope

`docs/archive/` — superseded docs, historical only, **never read in normal sessions**.
Out of scope: propellant sweep campaign (sims 030-085, `~/Desktop/*campaign*`) — ignore here.
