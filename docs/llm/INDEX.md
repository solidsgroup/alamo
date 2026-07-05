# docs/llm — chamber-gpu map (read this first)

Read this file, then `CONVENTIONS.md`, at session start.

## Status

Run `benchmark/status.sh` — sole source of branch truth (build gates, golden
compare, sanitizer smoke, open task folders). No prose status file exists.
Fresh clones: run `git config core.hooksPath .githooks` once (gates in `.githooks/`).

## Plan

`PLAN.md` is the only live plan (current phase, gate, next 3 tasks). Under
100 lines by design — prune it when it grows, don't append.

## References

- `CONVENTIONS.md` — session start/end checklist.
- `BUG_PATTERNS.md` — device bug classes already fixed once, enforced by
  `benchmark/lint_device_patterns.sh`.
- `VERSIONS.md` — semver ledger. `changelog/` — one entry per closed work item.

## Historical

`docs/archive/` holds every superseded roadmap/audit/status doc — historical
record only, **must not be read during normal sessions**.

## Out of scope

Propellant parameter-sweep campaign (sims 030-085, `~/Desktop/*campaign*` files) — unrelated to `chamber-gpu`/GPU-port work, ignore it here.
