# GPU-013: Use ReduceOps for aggregate kernel results
Transform status: draft
Class: correctness
Detection: Advisory regex; see `recognizers/table.csv`; confirm aggregate receiver, identity, and device reachability.
Invariant: Device aggregates require race-safe reductions with explicit identities; host assignment follows completion.
Port contract: The port supplies quantities, identities, ownership, ordering/tolerance rationale, and any required fusion decision. Record reduction choices and evidence for later review.
Transform: Use `ReduceOps`/`ReduceData` for device aggregation, then assign bounded results on host; fuse only as a measured variant.
Corpus example: Flame/chamber-gpu traction and displacement reductions are evidence of one application, not physics procedure.
Constraints: Preserve identities and semantics; do not fuse incompatible operators or change numerical behavior implicitly.
Verify: Instantiate `VALIDATION.md`; require `golden-regression`, `multi-box`, and `sanitizer` rows, plus decomposition checks and a written tolerance rationale.
Failure modes: Direct host writes can race or fault; wrong identities/types yield nondeterministic totals despite HMM success. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (AMReX ReduceOps/CUDA reduction anchors). Corpus: commits 76ae550e0388c96fd315724628942c2c5394426b, 02772a36997539b183c70ba195e469478b894007; `docs/llm/BUG_PATTERNS.md`; `docs/llm/changelog/2026-07-02-gpu-audit.md`.
