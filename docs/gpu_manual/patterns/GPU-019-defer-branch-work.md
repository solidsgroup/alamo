# GPU-019: Defer branch-only expensive work
Transform status: draft
Class: optimization
Detection: Advisory regex in `recognizers/table.csv`; inspect branch predicate, consumers, and device control flow.
Invariant: Device control flow must assign every consumed value and preserve synchronization; deferral is a measured hypothesis.
Port contract: Supply branch/consumer map, preserved assignment proof, baseline profile, and correctness/tolerance evidence. Record assumptions and dispositions in the port ledger.
Transform:
  Before: Expensive tensor/stencil work is computed before a branch that may not consume it.
  After: Evaluate the cheap predicate first and construct the value only in consuming branches.
Corpus example: chamber-gpu Elastic branch surgery is evidence of one optimization instance, not an integrator procedure.
Constraints: Preserve numerical ordering and all branch assignments; do not defer shared values without equivalent initialization.
Verify: Use `PERFORMANCE.md` baseline/optimization gates and `VALIDATION.md` correctness categories; optimization follows baseline-efficiency pass.
Failure modes: Undefined values, changed ordering, or no measured gain fail the hypothesis; inspect before widening it.
Evidence: Primary: `evidence/primary-sources.md` (CUDA control-flow semantics). Corpus: chamber-gpu commit `332ecffdd93a83def543d0348d592f815950cc5f`; cited Elastic/result paths.
