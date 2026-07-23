# TASK: multifin-frontier-rootcause

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 0 (diagnosis + probe runs; no src/ edits unless promoted)    |
| Model        | fable orchestrator, scout/executor legwork                   |
| Verification | partial-oracle (run_search_multifin.sh PASS/FAIL per point)  |
| Est. scope   | 0 src files; probe decks + logs only                         |
| Parallel-safe| yes (disjoint from GPU Phase 3 tasks)                        |

## Context budget

Read first: benchmark/status.sh output, this PLAN.md
Read: search_results.log, run_search_*.log (tails), input_half_multifin,
  run_search_multifin.sh, src/ files identified by scouts (psi_floor path,
  void modulus path, Newton/MLMG solver knobs)
Forbidden: docs/archive/*

## Objective

Frontier mapped 2026-07-08: half-multifin fails below void 18/14.5 MPa +
psi_floor 0.06. Two failure modes observed: P-peak cliff (floor-sensitive)
and late-burn web (void-sensitive). Root-cause WHY these fail — solver
(Newton overshoot / MLMG conditioning) vs physics (timestep, CFL) — then
push the frontier as low as physically possible using solver-side levers
(line_search, resync_coeffs, bottom solver, tolerances) rather than raising
the artificial stiffness floor.

## Oracle

Command: ./run_search_multifin.sh <K> <M> <F>  — PASS t=6.5 = point survives.
Covers: full-burn survival at np=8 half-domain.
Does NOT cover: physical fidelity of result (P trace sanity needs eyeball).

## Steps

1. Scout evidence: failure signatures per failed point (log tails).
2. Scout code: psi_floor + void modulus mechanism, available solver knobs
   on this branch (line_search? resync_coeffs? bottom solver?).
3. Analyze: assign each failure to a mechanism.
4. Probe runs (executor): targeted knob changes at failing points.
5. Report new frontier + mechanism writeup in results/RESULT.md.

## Closeout

- [ ] results/RESULT.md with mechanism + new frontier
- [ ] touch results/DONE
- [ ] SESSION_LOG.tsv line
