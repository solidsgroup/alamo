# GPU-009: Keep GPU temporaries alive through the launch
Transform status: draft
Class: correctness
Detection: Advisory manual inspection; see `recognizers/table.csv`; inspect local Fab ownership and asynchronous consumers.
Invariant: AMReX temporary storage must remain owned until all dependent asynchronous launches complete.
Port contract: The port supplies temporary owners, dependent streams, and the lifetime point at which retention ends. Record each dependency and release decision for later review.
Transform: Retain an AMReX Elixir for each locally created Fab consumed asynchronously, through completion of all dependent work.
Corpus example: Flame/chamber-gpu interpolation Elixir repair is evidence of a multi-box race, not a universal workload requirement.
Constraints: Do not add Elixir to CPU-only paths; ownership, not incidental synchronization, is the proof.
Verify: Instantiate `VALIDATION.md`; require `golden-regression`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: Stream-pool reuse can cause silent corruption or intermittent invalid access; single-box success is insufficient. Unresolved cases remain blocked; never infer correctness from a clean scan.
Evidence: Primary: `docs/gpu_manual/evidence/primary-sources.md` (AMReX Elixir/lifetime anchors). Corpus: commit c00f69086c6d1dc67bdf71e2bd174dbdcc953c85; `benchmark/archive/elixir_race_audit.md`; `docs/llm/BUG_PATTERNS.md`.
