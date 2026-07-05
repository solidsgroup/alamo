# CLAUDE.md - chamber-gpu

Every session, in order:
1. Run `benchmark/status.sh`. Do not trust any prose claim it contradicts.
2. Read `docs/llm/PLAN.md` (the only live plan) and the task folder you were
   assigned, nothing else.
3. Obey the task's Context budget. `docs/archive/` is historical and forbidden.

Rules:
- All new tasks start by copying `docs/llm/TASK_TEMPLATE.md` into a new
  `docs/agent_plans/YYYYMMDD-<name>/PLAN.md`. No untemplated work in src/.
- Hard rule: no kernel/perf optimization ships without a correctness pass
  (device lint, golden compare, compute-sanitizer on A1000).
- Hooks are required: `git config core.hooksPath .githooks` on fresh clones.
- Never edit changelog/ entries in place; append only.
- The propellant parameter-sweep campaign (~/Desktop, sims 030-085) is out of
  scope for this branch. Do not touch it.
- On session end: results/RESULT.md, touch results/DONE if complete, append one
  line to docs/llm/SESSION_LOG.tsv (date, task, model, tier, cost, outcome).

Key paths: src/Integrator/Flame.{cpp,H}, src/Model/Propellant/,
src/Model/Chamber/Ballistic.H, src/Operator/Elastic.*, src/Solver/Nonlocal/,
benchmark/ (harnesses), docs/llm/ (process), docs/agent_plans/ (tasks).
