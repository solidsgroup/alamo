# GPU-024: Const-qualify captured model accessors
Transform status: draft
Class: correctness
Detection: Advisory manual rule in `recognizers/table.csv`; inspect const receiver, device reachability, and accessor mutation.
Invariant: A device-reachable accessor must be device-callable and touch only device-safe state; declared read-only semantics must be preserved.
Port contract: Supply accessor call sites, mutability proof, device annotations, supported dimensions, and tolerance rationale. Record assumptions and dispositions in the port ledger.
Transform:
  Before: A device-reachable accessor lacks trailing `const` or performs hidden mutation.
  After: Mark the read-only accessor host/device callable and `const`; stage mutable state separately.
Corpus example: chamber-gpu Propellant/PNG edits are evidence of one accessor cleanup, not a Flame-specific recipe.
Constraints: Do not cast away mutation or alter model semantics; split genuinely mutable behavior from the device path.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: Compile errors, unsafe captures, or races from fake constness indicate failure; inspect the full call chain.
Evidence: Primary: `evidence/primary-sources.md` (CUDA callable/member rules). Corpus: chamber-gpu commit `587853d79432c0394928edb457683961ef5f08da`; cited Propellant/PNG paths.
