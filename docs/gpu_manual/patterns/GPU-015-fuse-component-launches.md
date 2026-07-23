# GPU-015: Fuse equivalent component launches
Transform status: draft
Class: optimization
Detection: Advisory regex in `recognizers/table.csv`; confirm loop bounds, component independence, and launch reachability by inspection.
Invariant: AMReX launch geometry and component ownership must remain valid; fusion is not a universal requirement.
Port contract: Supply the component ranges, dependency/order proof, CPU/golden oracle, and baseline/profile comparison. Record assumptions and dispositions in the port ledger.
Transform:
  Before: A loop launches identical index work once per component.
  After: Use one component-aware launch when bounds, dependencies, and reduction behavior are equivalent.
Corpus example: chamber-gpu `src/Operator/Operator.cpp:116-150` is evidence of one fusion instance, not procedure.
Constraints: This is a profiled hypothesis. Preserve ordering, synchronization, and component indexing; do not fuse differing bounds or reductions.
Verify: Use `PERFORMANCE.md` baseline and optimization gates plus `VALIDATION.md` correctness categories and a tolerance rationale; optimization is gated behind passing baseline efficiency.
Failure modes: Wrong component indexing, changed race/order, or register pressure can alter results or regress time; diagnose before broadening scope.
Evidence: Primary: `evidence/primary-sources.md` (AMReX launch semantics). Corpus: chamber-gpu commit `dc02baf1779b6668a6ea5725ade2e79b22118dfb`; cited source and result report.
