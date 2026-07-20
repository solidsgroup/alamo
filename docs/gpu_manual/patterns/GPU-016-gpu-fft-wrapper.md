# GPU-016: Route spectral transforms through the GPU wrapper
Status: draft
Class: correctness
Recognizer: regex: `amrex::FFT::R2C|AMReX_FFT\.H`
Applies: Integrators that construct AMReX FFT objects directly and duplicate spectral-layout/parallel-loop plumbing.
Transform:
  Before:
    `amrex::FFT::R2C fft(geom.Domain());` plus local layout and transform kernels.
  After:
    `Operator::Spectral::FFT fft(geom, refRatio(), lev);` and call its `Forward`, `Backward`, or `ParallelFor` API.
Constraints: Mandatory when the recognizer matches and the wrapper supplies the required hierarchy, layout, and boundary semantics; retain `ALAMO_FFT` compile guards and explicit unsupported behavior.
Verify: `make -j4`; run FFT-enabled Cahn-Hilliard/PFC/HeatConduction tests; expect CUDA and CPU builds to compile and spectral regression checks to pass. Repeat in the configured 2-D and 3-D modes where applicable, and retain the command/output in the build log.
Failure modes: Bypassing wrapper guards can fail when FFT support is disabled; wrong level/layout handling yields invalid transforms or numerical mismatch. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit 214780d1c3a8fcbe64e1ed4244dcdeb94ca2d8ca; `src/Operator/Spectral/FFT.H:0-0`, `src/Integrator/CahnHilliard.cpp`, `src/Integrator/PFC.cpp`.
