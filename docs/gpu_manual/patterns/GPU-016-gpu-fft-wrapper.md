# GPU-016: Route spectral transforms through the GPU wrapper
Status: verified
Class: correctness
Recognizer: regex: `amrex::FFT::R2C|AMReX_FFT\.H`
Applies: Integrators that construct AMReX FFT objects directly and duplicate spectral-layout/parallel-loop plumbing.
Transform:
  Before:
    `amrex::FFT::R2C fft(geom.Domain());` plus local layout and transform kernels.
  After:
    Include `Operator/Spectral/FFT.H`; construct `Operator::Spectral::FFT fft(geom, refRatio(), lev)`; allocate with `fft.MakeSpectralFab()`. An explicitly single/full-level method uses `fft.Forward(*field[lev],out)` and `fft.Backward(out,*field[lev])`; otherwise use hierarchy overloads `Forward(field,lev,out,0,0,time)` and `Backward(out,field,lev)`.
    Replace the local spectral-coordinate kernel with `fft.ParallelFor(box, lambda(m,n,p,omega2))`; the wrapper supplies `omega2`, layout, and inverse-transform scaling.
Constraints: Mandatory when the recognizer matches and the wrapper supplies the required hierarchy and boundary semantics. Retain `ALAMO_FFT` guards and explicit unsupported behavior. Remove direct layout access, local wave-number/DX state, and manual inverse scaling; `Backward` already normalizes. Do not use overloads to add AMR capability or change advancing levels; those are opt-in FEATURE changes.
Verify: `make -j4`; run FFT-enabled Cahn-Hilliard/PFC/HeatConduction tests; expect CUDA and CPU builds to compile and spectral regression checks to pass. Repeat in the configured 2-D and 3-D modes where applicable, and retain the command/output in the build log.
Failure modes: Bypassing wrapper guards can fail when FFT support is disabled. Calling nonexistent AMReX accessors on the wrapper fails compilation; retaining local scaling after `Backward` double-scales the solution; wrong level/layout handling yields numerical mismatch. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit 214780d1c3a8fcbe64e1ed4244dcdeb94ca2d8ca; `src/Operator/Spectral/FFT.H:0-0`, `src/Integrator/CahnHilliard.cpp`, `src/Integrator/PFC.cpp`.
