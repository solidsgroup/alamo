# GPU-010: Quarantine host-only loops during incremental porting
Transform status: draft
Class: scaffolding
Detection: Advisory manual inspection; see `recognizers/table.csv`; apply the policy in `docs/gpu_manual/ARCHITECTURE_POLICIES.md`.
Invariant: Host-only work must remain outside device closures and behind an explicit, observable boundary.
Port contract: The port supplies scope, user approval, fallback observability, and a retirement owner/date for every quarantine. Record the approval and retirement evidence for review.
Transform: Isolate an unconvertible loop or callback behind an explicit host-only branch while keeping device dispatch visible.
Corpus example: Flame/chamber-gpu quarantines are evidence of incremental-porting choices, not permission to add fallbacks.
Constraints: Temporary only; never expand the quarantine or hide required production work.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, and `multi-box` rows, plus the scaffolding disposition and a written tolerance rationale.
Failure modes: Hidden fallback, missing updates, or host calls reached from device code. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (CUDA host/device boundary anchors). Corpus: commits d522e1ac08729b306a215ad896d8a304983d55de; `docs/gpu_safe_ic_bc_matrix.md`.
