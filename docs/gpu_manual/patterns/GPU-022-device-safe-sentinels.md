# GPU-022: Replace host globals with device-safe sentinel access
Transform status: draft
Class: correctness
Detection: Advisory regex in `recognizers/table.csv`; confirm host-global reachability from device code and sentinel intent.
Invariant: Device code cannot dereference host-only storage; diagnostic sentinel behavior must remain observable and side-effect free.
Port contract: Supply sentinel ownership, host/device storage strategy, call-chain evidence, and oracle tolerance rationale. Record assumptions and dispositions in the port ledger.
Transform:
  Before: A device-reachable function returns a host-only global sentinel.
  After: Use a host/device accessor backed by storage valid for the active execution space.
Corpus example: chamber-gpu Set/Matrix changes are evidence of one sentinel migration, not a universal value choice.
Constraints: Preserve diagnostic semantics; do not mask invalid states or dereference host storage in device code.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: Device linking or illegal-address errors, or changed sentinel initialization, indicate failure; inspect the execution space.
Evidence: Primary: `evidence/primary-sources.md` (CUDA host/device memory rules). Corpus: chamber-gpu commit `f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da`; Set/Matrix paths.
