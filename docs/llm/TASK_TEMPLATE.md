# TASK: <short-name>
# Folder: docs/agent_plans/YYYYMMDD-<short-name>/
# Copy this template to PLAN.md in a new dated folder. Fill every <field>.
# Delete instructional comments (lines starting with #>) before handing to agent.

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | <0 / 1 / 2 / 3>                                              |
| Model        | <haiku / sonnet / opus>                                      |
| Verification | <full-oracle / partial-oracle / judgment>                    |
| Est. scope   | <files touched, rough line count>                            |
| Parallel-safe| <yes / no> <if yes: disjoint file set from what?>            |

#> Risk tiers:
#>   0 = docs, scripts, no src/         -> unattended, gates only
#>   1 = tests, benchmark harnesses     -> unattended, gates + spot-check diff
#>   2 = src/ non-solver                -> checkpoint at plan + before commit
#>   3 = solver/physics (Elastic, MLMG, -> checkpoint every step, adversarial
#>       Flame model terms, Newton)        review mandatory, line-by-line diff
#> Model routing: haiku=recon/search, sonnet=tier 0-2 execution, opus=tier 3
#>   or any task requiring root-cause reasoning.
#> Verification: full-oracle = an existing script proves correctness;
#>   partial-oracle = script covers some properties, rest needs judgment;
#>   judgment = no executable check exists. judgment + tier 3 tasks should be
#>   reconsidered before delegation.

## Operating rules

1. Read ONLY the files listed in Context budget. docs/archive/ is forbidden.
2. Every step's VERIFY must pass before acting on that step. On failure: STOP,
   report discrepancy, wait.
3. One commit per step unless stated. Message: `<area>: <what> (<task-folder>)`.
4. No scope expansion. New ideas go to a NOTES.md in this folder, not into code.
5. If required knowledge is missing (no source doc, no snippet, no repro), ask
   the user. Do not wing it without explicit permission.
6. Tier >= 2: stop at each checkpoint and print the checklist. Do not continue
   on your own assessment of success.

## Context budget

#> List exact paths. This is the agent's entire allowed reading set.
Read first: benchmark/status.sh output, this PLAN.md
Read: <file:line-range where possible>
Reference only if step names it: <files>
Forbidden: docs/archive/*, SUCCESS_BOOK-era docs, unrelated task folders

## Objective

<2-4 sentences. What exists now, what must exist after, why.>

## Oracle

#> The executable definition of done. If none exists, first step is to build it.
Command(s): <benchmark/..., exit 0 = pass>
Covers: <what properties this actually proves>
Does NOT cover: <known verification gaps; these require human review>

## Steps

### Step 1 - <name>
VERIFY:
```bash
<checks that preconditions hold>
```
DO: <precise action, file:line anchors>
CHECK: <command proving this step done>

### Step 2 - <name>
...

#> Size steps so each is one commit and one clean context. If a step needs
#> more than ~5 files in context, split it.

## Checkpoints (tier >= 2 only)

- [ ] After plan restatement: agent restates approach in its own words; human
      confirms before any edit
- [ ] Before each commit: diff summary + oracle output
- [ ] <additional tier-3 checkpoints per step>

## Adversarial review (tier 3 mandatory, tier 2 optional)

After implementation is complete and gates pass, human spawns a FRESH session
(no context from this folder) with:
  "Review commit <hash> on chamber-gpu. Assume it contains a defect. Find it.
   Check: device-lambda captures, elixir lifetimes, Eigen expression chaining,
   MLMG hierarchy freshness, tolerance changes, and whether tests were weakened
   to pass. Report findings only."
Reviewer findings go to results/REVIEW.md. Human adjudicates.

## Closeout

- [ ] Oracle passes; status.sh all green
- [ ] results/RESULT.md: what changed, evidence, deviations from plan
- [ ] changelog/ entry (append-only) + VERSIONS.md bump if release-worthy
- [ ] touch results/DONE
- [ ] Session log line appended: date, task, model, tier, cost, merged/abandoned
      -> docs/llm/SESSION_LOG.tsv
