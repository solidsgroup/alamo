# Conventions — how to work this map

Read this once at session start, alongside `INDEX.md`. It's short and stable;
it should rarely need to change.

## Read discipline

Orientation order when starting any task:

1. `INDEX.md` — pointers only.
2. `CURRENT.md` — what's in flight, with an explicit "Next / Test / Write
   results to" pointer.
3. The one named task file (e.g. `docs/agent_plans/<phase>/tasks/NNN-*.md`).

Stop there. Do not re-read `ROADMAP.md`, `benchmark/GPU_BRANCH_GUIDE.md`, or
any `changelog/` entry unless `CURRENT.md` is empty (nothing in flight — go
read `ROADMAP.md` to pick the next unstarted phase) or the task file is
missing something it needs. If you find yourself opening a third doc just to
figure out *what* to do or *how* to test it, the task file or `CURRENT.md`
was underspecified — fix the source file, don't just push through.

## Self-containment requirement

Every `tasks/*.md` must be actionable on its own, without requiring another
file to be opened. Use this shape:

```
## Task
## Why (one paragraph max)
## Test command (exact, runnable)
## Acceptance criteria
## Write results to: results/<NNN>-RESULT.md
```

A task file is not "done" (as a spec) if it requires opening `PLAN.md` or a
sibling task to be actionable.

## Write discipline

- `results/*.md` — append-only. One new file per task. Never edit a past
  result.
- `CURRENT.md` — overwrite only. Always ends with a concrete "Next:" pointer,
  or "phase complete."
- `benchmark/perf_regression.csv` (or equivalent) — append one row per
  measured change. Mechanical, no prose.
- `changelog/*.md` and `VERSIONS.md` — written once per phase boundary or
  version bump, not once per task. Never edited after the fact.
- If it isn't written to one of these files, it didn't happen. Don't leave
  load-bearing data (perf numbers, root causes, decisions) only in chat
  output.

## Data-recording discipline

- Every perf number in a `RESULT.md`, `changelog/`, or `perf/` entry must
  cite the exact command and source file that produced it — no unsourced
  numbers. (`benchmark/PERF_TRACKING.md` already does this well; match it.)
- A number that changes a roadmap decision (the D1/D2/D3/D4 pattern already
  in use — see `ROADMAP.md`) gets a one-line decision record in `VERSIONS.md`
  or the relevant `changelog/` entry, not just a buried mention in a results
  file.
- Version bumps (`VERSIONS.md`) require a real git tag on the commit being
  cited — a "version" is a checkout-able point, not just a narrative label.

## What never gets retroactively edited

`results/*.md`, `changelog/*.md`, and any file already cross-referenced by
path from one of those (this is why moving `benchmark/PHASE*.md` or
`docs/gpu_*.md` out of their current locations was deliberately avoided when
this map was built on 2026-06-22 — ~25 historical task/result files cite
them by their current path). New navigational docs point *into* these files
in place; they don't get relocated later just for tidiness.
