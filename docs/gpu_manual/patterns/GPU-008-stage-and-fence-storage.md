# GPU-008: Stage storage before its lifetime ends and fence the stream
Transform status: draft
Class: correctness
Detection: Advisory manual inspection; see `recognizers/table.csv`; inspect staging, ownership, and async lifetime boundaries.
Invariant: Asynchronous device work must not outlive the device-safe storage it references.
Port contract: The port supplies storage owners, launch streams, copy ordering, and the exact lifetime boundary requiring a fence. Record ownership and fence decisions for later review.
Transform: Stage host data into retained device-safe storage and fence outstanding work at local-storage destruction boundaries.
Corpus example: Flame/chamber-gpu Elixir and stream fixes are evidence of observed races, not a required synchronization schedule.
Constraints: Fence at ownership boundaries; do not substitute repeated global synchronization for correct ownership.
Verify: Instantiate `VALIDATION.md`; require `golden-regression`, `multi-box`, and `sanitizer` rows, plus lifetime evidence and a written tolerance rationale.
Failure modes: Intermittent illegal access, use-after-free, host faults, or stream-specific missed diagnostics. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (CUDA async/lifetime and AMReX stream anchors). Corpus: commits f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da, 54a941433b7582578cb5d56db794b1d648fb03cc, dd4054056aa5ef24bde948ac1b04013b5b106eb6; `docs/llm/BUG_PATTERNS.md`; `docs/llm/changelog/2026-07-02-gpu-audit.md`.
