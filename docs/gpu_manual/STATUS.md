# Manual and port status model

Status is evidence-scoped. No single `pass` summarizes the manual or a port.

## Pattern status

- `draft`: reviewed guidance without a successful closed-book transform.
- `file-verified`: successfully applied and reviewed on a target that was used
  to author or repair the pattern. This proves fit to that file only.
- `transfer-verified`: successfully applied to a frozen, previously unseen
  target that was not used to author or repair the pattern.
- `cross-family`: transfer-verified through completed-port gates in at least two
  distinct physics families, with both harvest records linked.

Historical phase-6 evidence supports `file-verified` only for GPU-007 and
GPU-016. The targets were repaired and retested, so neither result demonstrates
transfer. The other 24 patterns are draft; none is transfer-verified or
cross-family.

## Port status axes

Each port copies `templates/PORT_STATUS.csv` and reports `not-run`, `blocked`,
`fail`, or `pass` independently for `SCOPE`, `CLOSURE`, `INSPECTION`,
`VALIDATION`, `GPU_SAFETY`, `GPU_NATIVE_SHAPE`, `BASELINE_EFFICIENCY`, and
`HARVEST`. Optional optimization reports `not-attempted`, `measured-no-win`, or
`measured-win`. Scanner status is the count of remaining per-port `candidate`
rows, never a safety verdict.

## Frozen unseen-target gate

Before a run, record the target path, source revision/hash, applicable pattern
IDs, risk class, and expected artifacts. The session receives only INDEX,
ONBOARDING, applicable Tier 1 files, templates, toolchain contract, and unported
target material; it has no `chamber-gpu` access.

Stage A produces scope, compiler-first closure, per-file inspection, field map,
and proposed kernel graph. Stage B converts under the validation and shape
contracts. The lead checks artifacts, semantics, and invented behavior.

If the manual or a pattern is repaired using that run, the target joins the
authoring corpus. Retesting it may earn or retain `file-verified` only. A new
frozen target is required for `transfer-verified`; repair-then-retest on the
same file never qualifies. Gate selection must cover call/dispatch and
host-numerical chains, ownership/lifetime, reductions, dimensional/data
execution, API wrappers, scaffolding/onboarding, and measured optimization—not
only the smallest convenient files.

## Harvest obligation

Every completed port fills `templates/HARVEST.md`, updates `BLIND_SPOTS.md`, and
changes at least one manual artifact with its new evidence. Compiler/inspection
misses feed recognizer candidates; false positives receive durable dispositions;
new architecture questions update the policy home. A port with no harvested
manual change is incomplete.

