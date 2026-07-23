# GPU-017: Materialize chained matrix expressions
Transform status: draft
Class: correctness
Detection: Advisory regex in `recognizers/table.csv`; inspect expression-template type, device reachability, and scalar semantics.
Invariant: Device evaluation must preserve operation order, scalar type, and value lifetime; materialization is a device-compatibility tactic, not physics.
Port contract: Supply affected expression sites, supported matrix type, CPU/golden oracle, and tolerance rationale. Record assumptions and dispositions in the port ledger.
Transform:
  Before: A chained expression such as `F.inverse().transpose()` is evaluated in device code.
  After: Materialize each intermediate (`Finv`, then `FinvT`) before reuse.
Corpus example: chamber-gpu NeoHookean edits are corpus evidence for one expression-template failure mode, not a universal constitutive recipe.
Constraints: Apply only when compiler/runtime evidence identifies the chain; preserve mathematical order and do not change formulas.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `analytic-exact`, `golden-regression`, `restart-parity`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: Nested device expressions may fail compilation or launch; altered order changes results. Diagnose the exact signature first.
Evidence: Primary: `evidence/primary-sources.md` (CUDA device-expression constraints). Corpus: chamber-gpu commit `a5c1b2ddfd63e09848d69980d83496654d490ceb`; NeoHookean and BUG_PATTERNS paths.
