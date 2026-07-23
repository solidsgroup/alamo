# GPU-012: Report kernel errors through device flags
Transform status: draft
Class: correctness
Detection: Advisory regex plus compiler/runtime inspection; see `recognizers/table.csv`; apply `docs/gpu_manual/ARCHITECTURE_POLICIES.md`.
Invariant: Device code cannot invoke host-only diagnostics; host observation of asynchronous failure state occurs only after required synchronization.
Port contract: The port supplies flag ownership/reset, race-safe writes, stream coverage, and host reporting policy. Record the policy decision and failure evidence for review.
Transform: Set a device error flag in kernels, synchronize all relevant work, then invoke the host-side diagnostic path.
Corpus example: Flame/chamber-gpu `DeviceErrorFlag` fixes are evidence of stream races, not a complete error-policy specification.
Constraints: Never suppress, clamp, or silently continue on invalid data; check every relevant launch.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, `multi-box`, and `sanitizer` rows, plus injected-failure evidence and a written tolerance rationale.
Failure modes: Host/device compile errors, false negatives from partial synchronization, stale flags, or silent NaN propagation. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (CUDA synchronization/error anchors). Corpus: commit 54a941433b7582578cb5d56db794b1d648fb03cc; `docs/llm/changelog/2026-07-02-gpu-audit.md`; `docs/llm/BUG_PATTERNS.md`.
