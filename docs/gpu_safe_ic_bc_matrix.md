# GPU-Safe IC/BC Matrix

This is the Phase 0 support matrix for pure device-arena CUDA runs. "Supported" means the path does not write a device-arena fab from a host-only loop.

| Component | Status | Notes |
|---|---|---|
| `IC::BMP` | Supported | Initializes through `amrex::ParallelFor`; bitmap payload is copied to device memory. |
| `IC::Constant` | Supported | Initializes through `amrex::ParallelFor`. |
| `BC::Constant` | Supported | Constant scalar/vector BC paths are device kernels. |
| `BC::Operator::Elastic::Constant` | Supported | Elastic constant BC was made device-capturable for the chamber GPU path. |
| `IC::Expression` scalar | Supported | Scalar expression initialization uses `amrex::ParallelFor`. |
| `IC::StarAftGrain` | Supported | Star-aft grain (NAWC Motor No. 6). Station arrays are staged to `amrex::Gpu::DeviceVector` (Util::BMP pattern) and geometry is hoisted to locals so the `amrex::ParallelFor` lambda captures only device-valid data. |
| `IC::Expression` vector | Guarded unsupported | Vector expression path still uses `LoopConcurrentOnCpu`; CUDA builds abort before entering it. |
| `IC::PNG` | Guarded unsupported | Uses `LoopConcurrentOnCpu`; CUDA builds abort before entering it. |
| `IC::PSRead` | Guarded unsupported | Uses `LoopConcurrentOnCpu`; CUDA builds abort before entering it. |
| `IC::Trig` | Guarded unsupported | Uses `LoopConcurrentOnCpu`; CUDA builds abort before entering it. |
| `IC::Laminate` | Guarded unsupported | Uses `LoopConcurrentOnCpu`; CUDA builds abort before entering it. |
| `BC::Operator::Elastic::Expression` | Guarded unsupported | Uses `LoopConcurrentOnCpu`; CUDA builds abort before entering it. |

Canonical chamber GPU inputs should use BMP/Constant ICs and Constant elastic BCs until the guarded unsupported paths are ported to device kernels.
