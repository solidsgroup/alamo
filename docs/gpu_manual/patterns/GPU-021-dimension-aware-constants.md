# GPU-021: Use dimension-aware compile-time constants
Transform status: draft
Class: correctness
Detection: Advisory manual rule in `recognizers/table.csv`; inspect shared 2-D/3-D code and active component/index use.
Invariant: `AMREX_SPACEDIM` and dimension-aware declarations must make every compiled index valid; active-dimensional formulas remain unchanged.
Port contract: Declare supported dimensions, list guarded constants/terms, and provide 2-D/3-D builds plus oracle tolerance rationale. Record assumptions and dispositions in the port ledger.
Transform:
  Before: Shared device code unconditionally defines or indexes an inactive dimension.
  After: Guard declarations and terms with `AMREX_SPACEDIM`/`AMREX_D_DECL` and supply dimension-valid values.
Corpus example: chamber-gpu Stencil changes are evidence of one dimension audit, not a physics procedure.
Constraints: Do not hide required components behind runtime branches or alter active-dimensional formulas.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `analytic-exact`, `golden-regression`, `multi-box`, and `sanitizer` rows in every supported dimension with a written tolerance rationale.
Failure modes: Unprotected inactive terms cause compile/index errors; wrong constants change results. Diagnose before widening scope.
Evidence: Primary: `evidence/primary-sources.md` (AMReX dimensional macros). Corpus: chamber-gpu commit `5efce7cff26a13c2741b011d08735d2a679ef036`; Stencil path.
