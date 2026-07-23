# GPU-002: Replace device-reached host indirection with values and views
Transform status: draft
Class: correctness
Detection: Advisory high-recall declaration candidates in `recognizers/table.csv`; compiler/call-graph inspection must confirm device reachability, receiver type, allocation space, and runtime dispatch.
Invariant: Device code may use only device-valid object representations and addresses; host containers and host-created polymorphic receivers are not device views.
Port contract: Inventory every reachable pointer, owning/smart container, indexed model collection, dispatch alternative, variable-sized buffer, and lifetime. Prove which state becomes a copied value, retained device buffer/view, static alternative, or explicit host boundary.
Transform:
  Before: A device path dereferences a host-created base pointer or indexes container-backed model state and invokes runtime-selected behavior.
  After: Keep parsing/ownership on host; stage variable data into retained device storage; capture POD values/views; select concrete alternatives through a tuple/switch or separate static launch path.
Corpus example: The Gas-family chamber-gpu change is one successful instance; current-target Gas and Fracture containers are unported findings. The complete illustrative contract is `evidence/GPU-002-value-dispatch-example.md`.
Constraints: Do not blindly copy a `std::vector`, smart pointer, vtable-bearing object, or host address. Preserve selection, ordering, coefficients, species/material identity, and CPU behavior; unrelated host-only state is not a defect.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `analytic-exact` where available, `golden-regression`, `multi-box`, and `sanitizer` rows, plus selection/lifetime assertions and a written tolerance rationale.
Failure modes: HMM can mask host addresses; shallow copies retain invalid pointers; missing alternatives or stale buffers silently select wrong physics; a broad regex hit can be a legitimate host-only owner and requires disposition.
Evidence: Primary: `evidence/primary-sources.md` (CUDA storage/callability; AMReX device buffers/views). Corpus: chamber-gpu commit 3d382a5655c715c7483a15d53c3f8e5c908330af; `evidence/GPU-002-value-dispatch-example.md`; current-target evidence in `evidence/host-only-numerical-kernels.md`.
