# GPU-001: Make the complete device call chain callable
Transform status: draft
Class: correctness
Detection: Advisory build-error; see `recognizers/table.csv`; confirm the complete nvcc-reported call chain manually.
Invariant: Device code may call only host/device-callable functions and device-safe transitive callees.
Port contract: The port supplies the selected closure, call-chain inventory, and host-only boundary decisions. Record them for later review and reuse.
Transform: Add `AMREX_GPU_HOST_DEVICE` through every device-safe inline/helper path; split I/O, allocation, virtual, and other host-only work at an explicit boundary.
Corpus example: Flame/chamber-gpu annotation diffs are evidence of one exercised chain, not a universal procedure.
Constraints: Do not annotate routines with host-only effects; preserve semantics and host paths.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `golden-regression`, `multi-box`, and `sanitizer` rows, plus the port's written tolerance rationale.
Failure modes: Host/device diagnostics can remain in deeper template instantiations; host-only calls may fault at runtime or alter results. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (CUDA/AMReX callability anchors). Corpus: commits f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da, 3d382a5655c715c7483a15d53c3f8e5c908330af, 4afea70999397683597dc2fd8dbf4ad004e86bb1.
