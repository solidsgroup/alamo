# GPU-030: Initialize device-local aggregates before use
Transform status: draft
Class: correctness
Detection: Advisory regex in `recognizers/table.csv`; confirm conditional writes, active dimensions, and device reachability.
Invariant: Every consumed device-local value must be initialized on every control-flow path; neutral initialization must preserve semantics.
Port contract: Supply aggregate type, active components, control-flow proof, and validation tolerance rationale. Record assumptions and dispositions in the port ledger.
Transform:
  Before: A local vector/aggregate is conditionally populated after declaration.
  After: Initialize it to the mathematically neutral value before conditional writes.
Corpus example: chamber-gpu PNG sampling changes are evidence of one uninitialized-local fix, not a PNG procedure.
Constraints: Initialize before every path, preserve clamps/assignments, and do not zero values whose neutral element differs.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: Uninitialized components cause nondeterminism, NaNs, or sanitizer failures; wrong neutral values alter boundary behavior.
Evidence: Primary: `evidence/primary-sources.md` (CUDA initialization/device execution). Corpus: chamber-gpu commit `f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da`; PNG path.
