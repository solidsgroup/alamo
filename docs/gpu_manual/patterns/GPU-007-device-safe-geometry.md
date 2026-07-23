# GPU-007: Copy cell sizes into a device-safe value
Transform status: file-verified
Class: correctness
Detection: Advisory regex; see `recognizers/table.csv`; confirm raw `CellSize()` pointer reaches a launch.
Invariant: Device code must read geometry through device-safe values, not host-owned pointer state.
Port contract: The port supplies level, dimension ordering, units, and the value's lifetime. Record these choices and their inspection evidence for later review.
Transform: Copy cell sizes into a device-safe fixed-size value before launch and pass that value or its device-safe data.
Corpus example: Flame/chamber-gpu geometry conversion is evidence of one exercised path, not a universal discretization choice.
Constraints: Preserve configured dimensions, level, units, and numerical discretization. Record unresolved evidence explicitly.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, `multi-box`, and `sanitizer` rows, plus dimensional checks and a written tolerance rationale.
Failure modes: Host-pointer access, wrong level/index ordering, or pointer/value type errors. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (AMReX Geometry/device-data anchors). Corpus: commit d522e1ac08729b306a215ad896d8a304983d55de.
