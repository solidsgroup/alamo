# 2026-07-06 — efficiency-campaign

Task folder: `docs/agent_plans/20260706-efficiency-campaign/`

- **golden-compare gate un-broken (no src/ change needed).** The FAIL was a
  stale-binary selection bug: `ci_golden_compare.sh` picked its test binary by
  lexicographic `sort | tail -1`, which preferred Jun-29
  `bin/alamo-2d-perf-clang++` (predating `1cc3c32db`, so
  `elastic.solver.{line_search,resync_coeffs,nr_convergence}` were unparsed →
  ParmParse abort on `rod_and_tube_step2`) over the `bin/alamo-2d-g++` it had
  just built. Fixed: both legs now select the deterministic just-built name;
  `baseline_suite.py find_binary` picks newest-mtime. All 4 golden cases +
  NaN smoke PASS.
- **status.sh gate logs.** Gate output now tees to
  `benchmark/_gate_logs/<gate>.log` (gitignored); FAIL lines print the log path
  so the next session doesn't rediscover the failure reason from scratch.
- **Hygiene sweep.** Root strays deleted (`TASK_TEMPLATE.md`, `claude1.md`
  duplicates; stale `chamber_gpu_changes.diff`); `git gc` cleared 954 MiB pack
  garbage; SUPERSEDED headers on `remediation.md`,
  `docs/gpu_elastic_device_port_plan.md`, root `tasks/00*.md`;
  `benchmark/archive/README.md` live-plan pointer corrected to
  `docs/llm/PLAN.md`; 5 pre-template task folders retroactively closed
  (RESULT.md + DONE); `20260705-meta-workflow-install` DONE granted (its only
  blocker was this gate).
- **PLAN.md backlog.** Task 3.2b extended with the `Diagonal` per-component
  DDW-load hoist (Elastic.cpp:860/868); new backlog items 3.I (fuse Newton
  `norm0` reduction storm, Newton.H:419/572/816 → single ReduceOps pass) and
  an arena-policy A/B for Phase 5.
- **A100 A/B (PLAN task 3.1) NOT run** — NOVA rejects non-interactive auth.
  Copy-paste handoff:
  `docs/agent_plans/20260706-efficiency-campaign/results/A100_AB_HANDOFF.md`.
