# GPU-020: Specialize only the consumed tensor contraction
Transform status: draft
Class: optimization
Detection: Advisory regex in `recognizers/table.csv`; inspect hot-kernel reachability, consumed symmetry, and profile evidence.
Invariant: Index mapping, accumulation order, and scalar semantics must be preserved; specialization is not a universal requirement.
Port contract: Supply the consumed contraction shape, proof of symmetry assumptions, baseline/profile metrics, and tolerance rationale. Record assumptions and dispositions in the port ledger.
Transform:
  Before: A generic matrix contraction computes a broad operator surface in a hot device path.
  After: Implement only the consumed contraction with explicit terms or a narrow helper.
Corpus example: chamber-gpu Matrix4 edits are evidence of one specialization, not a general constitutive prescription.
Constraints: This is a profiled hypothesis. Do not specialize unused paths or change symmetry assumptions; reject register/code-size regressions.
Verify: Use `PERFORMANCE.md` baseline/optimization gates plus `VALIDATION.md` correctness categories; require baseline-efficiency pass first.
Failure modes: Wrong indices or accumulation order corrupt results; specialization without measured gain is unsuccessful.
Evidence: Primary: `evidence/primary-sources.md` (CUDA arithmetic/device execution). Corpus: chamber-gpu commit `9470889b14f10a902dab6dd9573deefc09972e8d`; Matrix4/result paths.
