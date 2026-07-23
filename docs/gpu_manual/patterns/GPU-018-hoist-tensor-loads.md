# GPU-018: Hoist reused field and tensor loads
Transform status: draft
Class: optimization
Detection: Advisory regex in `recognizers/table.csv`; inspect same-launch repetition, immutability, bounds, and device reachability.
Invariant: A load may be reused only while its value and memory validity remain unchanged; hoisting is not universally required.
Port contract: Supply repeated-load sites, proof of immutability/bounds, baseline profile, and correctness/tolerance evidence. Record assumptions and dispositions in the port ledger.
Transform:
  Before: A kernel repeatedly loads the same tensor/member for one cell.
  After: Load once into a device-local value and reuse it.
Corpus example: chamber-gpu Elastic edits show one `DDW` hoisting instance as evidence, not procedure.
Constraints: This is a profiled hypothesis. Do not hoist mutable or potentially out-of-bounds values; account for register pressure.
Verify: Use `PERFORMANCE.md` baseline and optimization gates plus `VALIDATION.md` correctness categories; optimization requires a passing baseline-efficiency result.
Failure modes: Stale values, illegal access, increased registers, or unchanged time indicate a failed hypothesis.
Evidence: Primary: `evidence/primary-sources.md` (CUDA memory/device execution). Corpus: chamber-gpu commit `9470889b14f10a902dab6dd9573deefc09972e8d`; cited Elastic/result paths.
