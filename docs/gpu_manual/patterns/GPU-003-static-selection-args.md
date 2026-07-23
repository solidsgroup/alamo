# GPU-003: Forward context through static selection
Transform status: draft
Class: correctness
Detection: Advisory regex; see `recognizers/table.csv`; inspect parser recursion and argument flow.
Invariant: Static selection must preserve required context and value lifetime through every recursion step.
Port contract: The port supplies constructor/context arguments, their ownership, and selected alternatives requiring them. Record the complete mapping for later review and reuse.
Transform: Forward `std::forward<Args>(args)...` through static-selection recursion to each concrete parser/model.
Corpus example: Flame/chamber-gpu parser forwarding is evidence of one context-loss defect, not a prescribed API.
Constraints: Preserve argument order and lifetime; do not introduce dynamic dispatch or host pointers. Record unresolved evidence explicitly.
Verify: Instantiate `VALIDATION.md`; require `strict-build` and `golden-regression` rows, plus context assertions and a written tolerance rationale.
Failure modes: Dropped arguments produce default/uninitialized context or unresolved overloads; stale recognizer hits must stop after forwarding exists. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (CUDA value/lifetime anchors). Corpus: commit 3d382a5655c715c7483a15d53c3f8e5c908330af.
