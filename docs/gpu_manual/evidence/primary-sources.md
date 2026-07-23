# Primary sources for invariant device semantics

This file grounds universal CUDA and AMReX claims. It does not prove that an
ALAMO call site satisfies a rule; compiler results, runtime evidence, and
recorded inspection provide that proof. Links were checked on 2026-07-21.

## Function space and device callability

- [CUDA Programming Guide: Intro to CUDA C++](https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/intro-to-cuda-cpp.html) defines host, device, and host/device function execution spaces, including member functions, functors, and lambdas.
- [CUDA Programming Guide: C/C++ language extensions](https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/cpp-language-extensions.html) is the primary reference for execution-space specifiers and CUDA C++ restrictions.

These sources ground the invariant behind GPU-001 through GPU-006, GPU-012,
GPU-017, GPU-022, GPU-024, and GPU-030: code and state reached by a device
execution path must be valid in that execution space. They do not prescribe an
ALAMO model hierarchy, boundary-condition representation, or error message.

## AMReX launches, dimensions, reductions, and explicit execution

- [AMReX: Overview of GPU Support](https://amrex-codes.github.io/amrex/docs_html/GPU.html) documents `ParallelFor`, component-aware launches, GPU reductions, stream synchronization, memory arenas, temporary lifetimes, and profiling guidance.
- [AMReX: Basics and dimensionality](https://amrex-codes.github.io/amrex/docs_html/Basics.html) documents `AMREX_SPACEDIM`, CPU/GPU `ParallelFor` behavior, and explicit `RunOn::Host`/`RunOn::Device` operations for `BaseFab`.

These sources ground GPU-007 through GPU-011, GPU-013, GPU-015, GPU-018
through GPU-021, GPU-023, and the baseline-efficiency contract. A port still
has to prove receiver types, component bounds, reduction identities, and its
chosen dimensional scope.

## Asynchrony, ownership, and synchronization

- [CUDA Programming Guide: Asynchronous Execution](https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/asynchronous-execution.html) explains that launches can return before work completes and that results and lifetimes require an appropriate synchronization relationship.
- [AMReX: GPU memory allocation, Elixir, and streams](https://amrex-codes.github.io/amrex/docs_html/GPU.html#memory-allocation) documents arena-backed storage, `Elixir`, async-safe temporaries, and current/all-stream synchronization.

These sources ground GPU-008, GPU-009, and the lifetime portions of GPU-012
and GPU-013. Synchronization placement remains a port inspection decision;
adding a global fence is not automatically correct or efficient.

## Runtime checking and profiling

- [NVIDIA Compute Sanitizer](https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html) defines memcheck, racecheck, initcheck, and synccheck and their diagnostic scope.
- [AMReX: Profiling with GPUs](https://amrex-codes.github.io/amrex/docs_html/GPU.html#profiling-with-gpus) describes the effect of asynchronous launches on timing and the role of whole-loop timing and GPU profiling tools.

These sources ground the sanitizer and trace categories in `VALIDATION.md`
and `PERFORMANCE.md`. A clean sanitizer run is evidence for the exercised
path only; it is not an analytic oracle, a convergence study, or a performance
result.

## Workload shape and resource evidence

- [CUDA C++ Best Practices Guide](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html) covers coalesced global-memory access, shared memory, register pressure, occupancy, branching/divergence, and allocation reuse.
- [CUDA Programming Guide: Writing SIMT Kernels](https://docs.nvidia.com/cuda/cuda-programming-guide/02-basics/writing-cuda-kernels.html) defines the execution/resource model and compiler resource reporting used to interpret block size, registers, shared memory, and occupancy.
- [AMReX: Overview of GPU Support](https://amrex-codes.github.io/amrex/docs_html/GPU.html) documents box-level versus whole-level launches, `MFIter` stream behavior, `TilingIfNotGPU`, device residency, reductions, arenas, block-size selection, and profiler limitations.
- [AMReX: MultiFab Tutorial](https://amrex-codes.github.io/amrex/tutorials_html/MultiFab.html) documents component-bearing `MultiFab`/`Array4` access and portable iteration forms.

These sources ground `GPU_NATIVE_SHAPE.md` and GPU-031. They establish what
must be measured or reasoned about; they do not prove that SoA, shared memory,
fusion, a particular block size, or any occupancy threshold is optimal for a
port.
