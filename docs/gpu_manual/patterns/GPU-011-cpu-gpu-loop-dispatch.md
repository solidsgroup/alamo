# GPU-011: Dispatch loops explicitly between CPU and GPU
Transform status: draft
Class: correctness
Detection: Advisory manual inspection; see `recognizers/table.csv`; inspect backend and ordering semantics.
Invariant: GPU-reachable work must execute in a backend compatible with its data, while the declared CPU path remains valid.
Port contract: The port supplies backend conditions, launch geometry, reduction/ordering semantics, and CPU fallback behavior. Record backend decisions and evidence for later review.
Transform: Use AMReX launch-region/ParallelFor dispatch with a device-callable body and an explicit CPU-safe path.
Corpus example: Flame/chamber-gpu loop dispatch is evidence of one backend split, not a universal launch shape.
Constraints: Do not move I/O, ordering-dependent work, or reductions without semantic proof. Record unresolved evidence explicitly.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: Races, divergent reductions, omitted CPU branches, or CPU/GPU numerical mismatch. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (AMReX ParallelFor/backend anchors). Corpus: commit d522e1ac08729b306a215ad896d8a304983d55de.
