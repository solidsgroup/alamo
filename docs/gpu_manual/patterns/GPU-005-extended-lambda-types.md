# GPU-005: Make captured parameter types visible to extended lambdas
Transform status: draft
Class: correctness
Detection: Advisory manual inspection; see `recognizers/table.csv`; compiler access diagnostics are corroborating evidence.
Invariant: Extended device lambdas require accessible, device-safe parameter and capture types.
Port contract: The port supplies launch scope, parameter visibility, ownership, and encapsulation decisions. Record those decisions in the inspection artifacts for later review.
Transform: Expose or copy a device-safe parameter aggregate before launch, capture it by value, and use local fields.
Corpus example: Flame/chamber-gpu lambda visibility fixes are evidence of one compiler limitation, not an API mandate.
Constraints: Do not broaden public API unnecessarily; never retain hidden host `this` or inaccessible state.
Verify: Instantiate `VALIDATION.md`; require `strict-build` and `golden-regression` rows, and record diagnostics plus a written tolerance rationale.
Failure modes: nvcc can reject private/protected captured types or runtime can dereference host state. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (CUDA extended-lambda anchors). Corpus: commits d522e1ac08729b306a215ad896d8a304983d55de, 4afea70999397683597dc2fd8dbf4ad004e86bb1.
