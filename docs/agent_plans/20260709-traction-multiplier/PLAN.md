# TASK: traction-multiplier
# Folder: docs/agent_plans/20260709-traction-multiplier/

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 2                                                             |
| Model        | sonnet                                                        |
| Verification | partial-oracle                                                |
| Est. scope   | 2 files (Flame.H, Flame.cpp), ~4 lines                        |
| Parallel-safe| yes, disjoint from 20260708-multifin-frontier-rootcause       |

## Operating rules

1. Read ONLY the files listed in Context budget. docs/archive/ is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: STOP,
   report discrepancy, wait.
3. One commit per step unless stated. Message: `<area>: <what> (<task-folder>)`.
4. No scope expansion. New ideas go to a NOTES.md in this folder, not into code.
5. If required knowledge is missing, ask the user. Do not wing it.
6. Tier >= 2: stop at each checkpoint and print the checklist.

## Context budget

Read first: benchmark/status.sh output, this PLAN.md
Read: src/Integrator/Flame.H:130-150, src/Integrator/Flame.cpp:255-275,475-490
Reference only if step names it: none
Forbidden: docs/archive/*, unrelated task folders

## Objective

`elastic.traction` was found to be a dead parameter whenever
`elastic.traction_from_chamber = 1` (the branch at Flame.cpp:481 ignores it
entirely, using `chamber.pressure` instead). The user wants a way to scale
the *applied* traction regardless of its source (constant or chamber-derived)
via a new scalar `elastic.traction_multiplier`, default 1.0, set to 1.16 for
the next confirm batch.

## Oracle

Command(s): `benchmark/status.sh` (device-lint, golden-compare, a100-sanitizer)
Covers: bit-identical behavior at default multiplier=1.0 (golden compare),
device-lambda capture correctness (device-lint), no illegal device memory
access (compute-sanitizer).
Does NOT cover: physical correctness of traction*1.16 runs — that's a new
confirm batch the user will judge by eye (P-trace / stability), not the gate.

## Steps

### Step 1 - add traction_multiplier param
VERIFY: `grep -n "traction_from_chamber" src/Integrator/Flame.H src/Integrator/Flame.cpp`
DO:
  - Flame.H:141-144 area: add `Set::Scalar traction_multiplier = 1.0;` next to
    `traction_from_chamber`.
  - Flame.cpp:269 area: add
    `pp_query_default("elastic.traction_multiplier", value.elastic.traction_multiplier, 1.0);`
  - Flame.cpp:481: change
    `const Set::Scalar traction = elastic.traction_from_chamber ? chamber.pressure : elastic.traction;`
    to
    `const Set::Scalar traction = (elastic.traction_from_chamber ? chamber.pressure : elastic.traction) * elastic.traction_multiplier;`
CHECK: `grep -n "traction_multiplier" src/Integrator/Flame.H src/Integrator/Flame.cpp`

### Step 2 - gate pass
VERIFY: build succeeds
DO: rebuild, run `benchmark/status.sh`
CHECK: all three gates PASS (default multiplier=1.0 must stay bit-identical
to pre-change golden output since 1.0 is a no-op multiply)

## Checkpoints (tier >= 2 only)

- [ ] After plan restatement: agent restates approach; human confirms before edit
- [ ] Before commit: diff summary + oracle output

## Adversarial review

Skipped by user request (small, mechanical, reviewed inline this session).

## Closeout

- [ ] Oracle passes; status.sh all green
- [ ] results/RESULT.md written
- [ ] touch results/DONE
- [ ] Session log line appended to docs/llm/SESSION_LOG.tsv
