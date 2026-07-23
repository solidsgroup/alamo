# GPU-004: Remove virtual boundary-condition dispatch from kernels
Transform status: draft
Class: correctness
Detection: Advisory regex; see `recognizers/table.csv`; manually confirm virtual BC reachability.
Invariant: Device kernels cannot rely on host virtual dispatch or host-owned BC state.
Port contract: The port supplies BC alternatives, POD state, CPU behavior, and component/face/ghost semantics. Record the mapping and boundary decisions for later review.
Transform: Stage concrete BC state before launch and call static device-safe evaluation; retain virtual dispatch on the CPU path.
Corpus example: Flame/chamber-gpu elastic BC staging is evidence, not a universal BC API.
Constraints: Keep construction/parsing host-side; preserve ordering and boundary semantics. Record unresolved evidence explicitly always.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: Device vtable or host-state faults, wrong condition selection, or altered boundary ordering. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (CUDA device-callable/state anchors). Corpus: commits 1da57132e730ebaf97014ff63bc74f0b72cbf870, 54a941433b7582578cb5d56db794b1d648fb03cc.
