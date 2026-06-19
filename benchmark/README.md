# Flame/chamber GPU build + benchmarking

This directory holds the tooling to build the `Integrator::Flame` (chamber)
solver for NVIDIA GPUs, run it with overlapped async IO, and benchmark CPU vs
GPU wall-clock with performance flame graphs.

## 1. Build

CPU (baseline, with profiling):

```bash
./configure --profile
make
# -> bin/alamo-2d-<comp>
```

GPU (CUDA, with profiling). Requires the **CUDA toolkit** (`nvcc`) and an
NVIDIA GPU; `--cuda` auto-detects the local compute capability via `nvidia-smi`:

```bash
./configure --cuda --profile          # or --cuda 86 to pin sm_86
make
# -> bin/alamo-2d-cuda86-<comp>
```

Both binaries coexist (the POSTFIX differs).

## 2. Run with async IO + profiling

Append `benchmark/async_profile.inputs` to your chamber input, or pass the
params on the command line. Async output (`amrex.async_out=1`) hands plotfile
writes to a background thread so the GPU marches step N+1 while the CPU drains
the step-N write -- the requested compute/IO overlap, with no source changes
(`Integrator::WritePlotFile` already routes through
`amrex::WriteMultiLevelPlotfile`, which honors `async_out`).

```bash
mpiexec -np 1 ./bin/alamo-2d-cuda86-g++ input \
    amrex.async_out=1 tiny_profiler.device_synchronize_around_region=1
```

## 3. Benchmark CPU vs GPU

```bash
benchmark/benchmark_gpu_cpu.sh input 50      # 50 steps on each available binary
```

This runs the same input on the CPU and CUDA binaries, captures the AMReX
TinyProfiler per-region wall-clock tables (the holdup breakdown: `Flame::Advance`
compute vs `WritePlotFile` IO vs `Operator::Elastic` solve vs `Flame::Regrid`),
wraps the GPU run in `nsys` when present, and prints a side-by-side speedup
table. Artifacts land in `benchmark/results_<stamp>/`.

### Flame graphs

* **GPU (best):** the script captures `gpu_trace.nsys-rep`. Open it in NVIDIA
  Nsight Systems for the CPU+GPU timeline, or `nsys stats` it. This is the
  truest call-stack/kernel flame graph for GPU runs.
* **CPU:** for a nested call-stack flame graph, run under `perf` and Brendan
  Gregg's FlameGraph:
  ```bash
  perf record -F 997 -g -- mpiexec -np 1 ./bin/alamo-2d-g++ input max_step=50
  perf script | stackcollapse-perf.pl | flamegraph.pl > cpu.svg
  ```
* **Cross-platform (flat):** `tinyprofiler_to_folded.py` turns a TinyProfiler
  table into folded stacks; `benchmark_gpu_cpu.sh` feeds it to `flamegraph.pl`
  if available. Flat (one bar per `BL_PROFILE` region) because TinyProfiler is
  not hierarchical.

## Holdup / wall-clock tracking

`tiny_profiler.device_synchronize_around_region=1` attributes GPU kernel time to
the correct `BL_PROFILE` region. The TinyProfiler tables then show where time
goes (compute / IO / solver / regrid). The driver script also records total
wall-clock per run in `summary.txt`. For finer GPU/RAM/IO attribution use the
`nsys` trace (kernel time, memcpy/HtoD-DtoH stalls, CPU gaps).

## Optimization knobs (making the GPU win big)

- **Never build with `--debug`.** `--use_fast_math` (fast `exp`/`tanh`/`pow`
  intrinsics — a large win for the Arrhenius mobility and heat-flux kernels)
  lives only in the non-debug nvcc path, alongside `--ptxas-options=-O3`,
  `-O3 -finline-limit`, and `-maxrregcount=255`. `nova_flame_gpu.slurm` builds
  non-debug by default.
- **`MODE=fast` vs `MODE=bench`** (env var in `nova_flame_gpu.slurm`):
  - `bench` (default) — `--profile` + `device_synchronize_around_region=1` +
    managed arena. Accurate per-region holdups, but the syncs serialize the GPU.
  - `fast` — no profiler, no per-region sync, **pure device arena**
    (`the_arena_is_managed=0`, no host page-migration). This is the
    "blow CPU out of the water" production build. Run `MODE=fast sbatch ...`.
- **Problem size matters.** GPU advantage grows with cells-per-kernel. The
  center-bore default is `n_cell = 64 64` / `max_level=3` (small — kernel-launch
  overhead dominates on GPU). To actually showcase the GPU, scale the grid:
  `amr.n_cell="256 256" amr.max_level=4`. Pass on the launch line or in a copy
  of the input.
- **One GPU is usually right for 2D.** Multi-GPU (`--ntasks=N --gres=gpu:N`)
  helps only when the domain is large enough that each GPU still has plenty of
  boxes; for a 64^2 grid a single GPU wins (no halo-exchange overhead).

## Status

Done and CPU-verified on branch `chamber-gpu`:
- CUDA build system (`configure --cuda`, Makefile nvcc indirection) from `origin/gpu`.
- Propellant models `Homogenize` + `FullFeedback` made `AMREX_GPU_HOST_DEVICE`
  and `const`-correct (matching the already-ported `Constant`/`PowerLaw`).
- `Integrator::Flame` (`Flame.cpp`): device-clean kernels, all `this`-captures
  eliminated, `Integrate` chamber accumulation rewritten as an `amrex::ReduceOps`
  sum (the per-cell `+=` into members was a GPU race).
- **Shared dependency tree ported** from `origin/gpu`: `Base/Mechanics.H`,
  `Solver/Nonlocal/Newton.H`, `BC/Operator/Elastic/*`, `IC/*`, plus the
  foundational `Set/Base.H` (device `Position` overload), `Set/Matrix4_*`,
  `Util/{Util,BMP,PNG}`. Every device `ParallelFor` lambda is now `this`-free;
  the remaining `[=]` this-captures are intentional **host** `LoopConcurrentOnCpu`
  loops (ICs / elastic-BC boundary handling), which nvcc accepts. The full binary
  links on CPU and the center-bore run is **bit-for-bit identical** to before the
  port.

**Only remaining step:** final `nvcc` compile + on-GPU runtime validation on a
CUDA-toolkit machine (NOVA). The dev box used here has the driver `libcuda.so`
but no `nvcc`, so the `--cuda` build could not be compiled/run locally.
