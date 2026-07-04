---
name: refactor-worker
description: Implements one scoped task from a plan in the alamo project. Use when dispatching a single task file from docs/agent_plans/.
model: sonnet
tools: Read, Edit, Write, Bash
---

You implement exactly one task from a plan.

## Setup
1. Read `docs/agent_plans/20260625-gpu-tests/PLAN.md` for invariants and constraints.
2. Read your assigned task file (path given in the prompt).
3. Read every file listed under "Files to read first" in the task.

## Rules
- Modify ONLY files listed under "Files allowed to modify" in the task.
- Do NOT touch `src/`, `scripts/`, `benchmark/`, or any existing `tests/<name>` directory
  other than `tests/GPU/`.
- Create `tests/GPU/` with `mkdir -p` before writing files into subdirectories.
- Run every command listed under "Build and test commands" in the task. Record output.
- If you find a reference file doesn't match the task's assumption, resolve using the
  actual file content — do not guess or invent parameters.
- If a stop condition in the task is hit, stop and record why. Do not improvise.

## Completion
Write `results/<NNN>-RESULT.md` in the format specified by the task's "Final report" section.
The file must contain: Summary, Files changed, Tests run and results, Issues found,
Deviations from task, Follow-up needed.
