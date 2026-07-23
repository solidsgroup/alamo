# GPU-006: Hoist kernel captures into device-safe locals
Transform status: draft
Class: correctness
Detection: Advisory regex; see `recognizers/table.csv`; confirm member/host capture rather than a shadowed local.
Invariant: Device closures capture only device-safe values; hidden `this` and host pointers are not device state.
Port contract: The port supplies capture inventory, ownership/lifetime proof, and explicit host/device boundaries. Record each capture disposition for later review and reuse.
Transform: Copy scalar or device-safe aggregates before launch and capture locals by value.
Corpus example: Flame/chamber-gpu member-hoisting diffs are evidence of one closure, not a port recipe.
Constraints: Exclude intentional device-safe locals; never capture allocators, polymorphic objects, or host references.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: nvcc access errors, host-pointer faults masked by HMM, or changed results from wrong ownership/dimension. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (CUDA capture and AMReX launch-region anchors). Corpus: commit d522e1ac08729b306a215ad896d8a304983d55de; `docs/gpu_device_capture_conventions.md`; `docs/llm/BUG_PATTERNS.md`.
