# GPU-016: Route spectral transforms through the GPU wrapper
Transform status: file-verified
Class: correctness
Detection: Advisory regex/manual review in `recognizers/table.csv`; confirm direct FFT construction and wrapper applicability.
Invariant: Device code may call only device-safe APIs, and transforms must preserve layout, level, synchronization, and normalization semantics.
Port contract: Supply the wrapper API version, selected single/full-level overload, FFT guards, spectral allocation, and oracle tolerance rationale. Record assumptions and dispositions in the port ledger.
Transform:
  Before: Direct AMReX FFT construction duplicates layout, spectral loops, or inverse scaling.
  After: Use the project FFT wrapper for allocation, Forward/Backward, and wrapper-owned `ParallelFor`/scaling; retain explicit level semantics.
Corpus example: chamber-gpu Cahn-Hilliard/PFC changes document one validated wrapper migration, not a universal physics procedure.
Constraints: Do not add AMR levels or alter advancing policy through overload selection; preserve unsupported-feature behavior.
Verify: Instantiate `VALIDATION.md`; require `analytic-exact`, `golden-regression`, `restart-parity`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: Wrong overload/layout or duplicate scaling causes compile or numerical failures; missing guards breaks non-FFT builds.
Evidence: Primary: `evidence/primary-sources.md` (AMReX FFT/device semantics). Corpus: chamber-gpu commit `214780d1c3a8fcbe64e1ed4244dcdeb94ca2d8ca`; wrapper/source paths.
