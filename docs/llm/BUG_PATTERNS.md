# Device-Side Bug Patterns

Raw material mined from `docs/archive/SUCCESS_BOOK.md` and `docs/llm/changelog/`
during REMEDIATION_PLAN Phase 1. Source of truth for `benchmark/lint_device_patterns.sh`
(Phase 3). Do not re-derive rules from the prose descriptions below alone — Phase 3
must validate each rule against the actual pre-fix code shown here.

---

## 1. Elixir / use-after-free (local FArrayBox freed while an async kernel still reads it)

**Site:** `src/Operator/Operator.cpp:728`, `Operator<Grid::Node>::interpolation()`.

**Pre-fix pattern:** a local `FArrayBox tmpfab` created inside an `MFIter` loop,
populated by an async `ParallelFor`, goes out of scope at the end of the loop body.
AMReX's arena frees the device memory immediately, but a kernel on a different CUDA
stream (MFIter cycles a stream *pool*) may still be reading it. Single-box runs never
hit it (only one stream); multi-box MLMG did, as silently corrupted interpolation.

**Fix:**
```cpp
FArrayBox tmpfab(...);
// ... async ParallelFor populates tmpfab ...
tmpfab.elixir();   // <-- keeps device memory alive until all touching streams sync
```

**Rule of thumb:** any locally-created `FArrayBox`/`MultiFab` temporary fed to an async
device kernel inside an `MFIter` loop needs `.elixir()`. A fab obtained via
`multifab[mfi]` (not locally created) is already owned by the `MultiFab` and is safe.

**Related instance (same class, different surface):** `Util::DeviceErrorFlag::value()`
(`src/Util/Util.H:68-78`) read back the error flag without a full stream sync, racing
the same MFIter stream pool — a silent false-negative rather than a crash. Fixed by a
`streamSynchronizeAll()` before the read, centralized in `Util::AbortIfDeviceError`
(`Util.H:148-151`); same idiom reused for a `dsol_mf` UAF at `Newton.H:459`.
Cite: `docs/llm/changelog/2026-07-02-gpu-audit.md`.

---

## 2. Host member / implicit `this` capture inside a device lambda

**Site A:** `Base::Mechanics::Integrate` wrote to host members (`trac_hi`, `disp_hi`)
from inside a device kernel — silent on the HMM-backed local A1000, CUDA error 700
(illegal memory access) on real A100 hardware.

**Fix A:** replace the direct host-member write with `amrex::ReduceOps`, accumulating
into a device-local result and writing to the host member only after the kernel
completes (pattern already used in `Flame::Integrate`).

**Site B (dormant, unfixed):** `PhaseFieldMicrostructure::Integrate`
(`src/PhaseFieldMicrostructure.cpp:449-483`) accumulates host members via `+=` inside a
device lambda's implicit `this` capture — the same pre-fix pattern as Site A. Inert
today only because this integrator isn't in the GPU-supported closure; must be ported
to `ReduceOps` before it ever is. Cite: `docs/llm/changelog/2026-07-02-gpu-audit.md`.

**Site C:** `IC::StarAftGrain` (`src/IC/StarAftGrain.H`, `Add()`'s `ParallelFor`) had
three variants of this trap at once:
1. `std::vector<...>.data()` host pointers captured by value and dereferenced on device.
2. `geom[lev]` referenced inside `[=]`, which captures `this` and dereferences a host
   reference member on device.
3. `Set::Constant::Pi`, an `extern const` host global, referenced in device code
   (undefined there — compiles fine on CPU, only surfaces under nvcc).

**Fix C:** stage vector payloads into `amrex::Gpu::DeviceVector` members (capture the
device `.data()` pointer, mirroring `Util::BMP`); hoist `ProbLo`/`CellSize` into local
`Set::Vector` values before the lambda; replace `Set::Constant::Pi` with a local
`constexpr Set::Scalar pi` literal.

**Rule of thumb:** inside `AMREX_GPU_DEVICE`/`[=] AMREX_GPU_DEVICE` lambda bodies,
never reference `this->member`, a reference/pointer member (including implicitly via
a bare member name), or an `extern`/host global. Hoist every value needed into a local
before the lambda, and stage container payloads through `amrex::Gpu::DeviceVector`.

---

## 3. Chained Eigen expression templates evaluated on device

**Site:** `src/Model/Solid/Finite/NeoHookean.H`, 3D branches of `DW` (line 75) and
`DDW` (line 127).

**Pre-fix pattern:**
```cpp
Eigen::Matrix3d FinvT = F.inverse().transpose();   // nested Transpose<Inverse<Matrix3>>
```
This single chained expression faults on device with `CUDA error 719: unspecified
launch failure` on the very first elastic solve — not on CPU, not in 2D (the 2D
branches hand-roll every inverse and never form a chained `.inverse().<op>`
expression). Bisection confirmed each operation is device-safe in isolation
(`F.inverse()`, `.determinant()`, `(F*F^T).trace()`); only the *composed* expression
faults.

**Fix:** materialize the intermediate before applying the next operation:
```cpp
Eigen::Matrix3d Finv = F.inverse();
Eigen::Matrix3d FinvT = Finv.transpose();
```

**Rule of thumb:** any chained Eigen expression inside an `AMREX_GPU_DEVICE` context
(`.inverse().transpose()`, `.cross().normalized()`, etc.) must be broken into separate
statements (or `.eval()`'d) so each intermediate is materialized before the next
operation composes on top of it.
