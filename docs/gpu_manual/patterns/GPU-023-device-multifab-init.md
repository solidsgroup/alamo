# GPU-023: Initialize device-backed MultiFabs on device
Transform status: draft
Class: correctness
Detection: Advisory regex in `recognizers/table.csv`; inspect receiver type, data arena, launch reachability, and execution mode.
Invariant: Operations on device-resident data must execute in a compatible space with preserved ordering, ranges, and synchronization.
Port contract: Supply receiver/type evidence, selected `RunOn` policy, supported dimensions, and validation tolerance rationale. Record assumptions and dispositions in the port ledger.
Transform:
  Before: Fab initialization or algebra relies on implicit execution mode.
  After: Select explicit device execution for initialization and dependent operations where the data path requires it.
Corpus example: chamber-gpu Operator/Elastic edits are evidence of one MultiFab execution-mode fix, not a procedure for every Fab.
Constraints: Confirm receiver semantics; do not change component ranges, ordering, or synchronization, and do not apply to host-only data.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: Host execution on device memory, stale values, or race-dependent residuals indicate failure; inspect type and arena first.
Evidence: Primary: `evidence/primary-sources.md` (AMReX Fab execution policies). Corpus: chamber-gpu commit `d522e1ac08729b306a215ad896d8a304983d55de`; cited source paths.
