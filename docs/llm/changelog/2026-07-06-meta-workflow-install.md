# 2026-07-06 — meta-workflow-install (process)

Installed the session workflow scaffold. Process only; no `src/` changes, no
VERSIONS.md bump.

- `CLAUDE.md` (repo root) — session entry point / mandatory checklist (Step 1, `ec8c47b06`).
- `docs/llm/TASK_TEMPLATE.md` — copied for every new task folder (Step 1).
- `docs/llm/SESSION_LOG.tsv` — 6-column append-only session telemetry (Step 2, `c8e44736c`).
- `.githooks/pre-commit` — added a doc-budget guard after the device-lint (Step 3, `18f758adb`):
  live top-level `docs/llm` docs capped at 600 lines (excludes append-only `changelog/`
  and historical `perf/`; 313/600 now), and blocks new `.md` outside
  `docs/archive/ | docs/agent_plans/ | docs/llm/changelog/`.
- `docs/llm/INDEX.md` — now points at the three new pieces, kept at 29 lines (Step 4, `42896dc2e`).

Deviations (both necessary, see `docs/agent_plans/20260705-meta-workflow-install/results/RESULT.md`):
the guard's find scope was narrowed (would have bricked all commits at 1534>600) and its
allowlist path corrected `changelog/`→`docs/llm/changelog/` (would have blocked this entry).

Note: `status.sh` `golden-compare` is a pre-existing FAIL (`rod_and_tube_step2/cpu` exit 6),
unrelated to this task; `DONE` held pending that gate.
