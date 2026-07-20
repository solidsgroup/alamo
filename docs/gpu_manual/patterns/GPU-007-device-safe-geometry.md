# GPU-007: Copy cell sizes into a device-safe value
Status: draft
Class: correctness
Recognizer: regex: `(?:const )?(?:Set::Scalar|amrex::Real)\s*\*\s*DX\s*=\s*geom\[[^]]+\]\.CellSize\(\)`
Applies: A launch captures the raw pointer returned by `Geometry::CellSize()`.
Transform:
  Before:
    `const Real* DX = geom[lev].CellSize();` then pass `DX` to a device lambda.
  After:
    `Set::Vector DX(geom[lev].CellSize());` and pass `DX.data()` (or indexed value) in device-safe calls.
Constraints: Mandatory when recognizer matches. Preserve AMREX_SPACEDIM ordering and units; do not change discretization or silently copy a different level.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected GPU build succeeds and gradient/divergence smoke values agree with CPU.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. illegal host-pointer access, wrong dimensional indexing, or compile failure from pointer/value mismatch. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise. Include launch-region behavior, ownership, and dimensional assumptions in review; the transform is complete only when those remain explicit. Confirm dimensional and reduction identities explicitly.
Evidence: commits d522e1ac08729b306a215ad896d8a304983d55de
