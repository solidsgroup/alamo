# Flame/chamber GPU build + benchmarking

This directory holds the tooling to build the `Integrator::Flame` (chamber)
solver for NVIDIA GPUs, run it with overlapped async IO, and benchmark CPU vs
GPU wall-clock with performance flame graphs.

## 0. NOVA quickstart (ISU)

On the NOVA **login node**:

```bash
# 1. download + build both A100 (sm_80) and H200 (sm_90) binaries (submits a
#    CPU build job on the `nova` partition, account brunnels). Builds in the
#    current directory -- cd to where you want it first:
mkdir -p /work/brunnels/jackplum/alamo && cd /work/brunnels/jackplum/alamo
sh /path/to/build_alamo_nova.sh          # HTTPS clone by default (no SSH key)
#    (override: ALAMO_DIR=/abs/path  REPO_URL=git@github.com:solidsgroup/alamo.git ...)

# 2. once the build job finishes (squeue -j <id>), run on a GPU:
sbatch --gres=gpu:a100:1 --ntasks=1            benchmark/nova_flame_gpu.slurm   # A100
GPU_TYPE=h200 sbatch --gres=gpu:h200:1 --ntasks=1 benchmark/nova_flame_gpu.slurm # H200
MODE=fast     sbatch --gres=gpu:a100:1 --ntasks=1 benchmark/nova_flame_gpu.slurm # max speed

# 3. CPU baseline for the comparison:
sbatch benchmark/nova_flame_cpu.slurm
```

Account `brunnels`, partition `nova`, email pre-filled. If GitHub SSH isn't set
up on NOVA, set `REPO_URL=https://github.com/solidsgroup/alamo.git`. Module names
(`cuda`/`gcc`/`openmpi`) may need tweaking to match `module avail`.

## 1. Local workstation build

For quick testing on the local NVIDIA GPU, build a binary for the detected
compute capability. On this workstation that is currently `sm_86`, so the helper
produces `bin/alamo_gpu-2d-cuda86-g++` instead of the NOVA-only `cuda80`/`cuda90`
binaries:

```bash
benchmark/build_alamo_local_gpu.sh          # detects local arch and smoke-tests
PROFILE=1 benchmark/build_alamo_local_gpu.sh # profiler-enabled local build
SMOKE=0 benchmark/build_alamo_local_gpu.sh   # build only
```

Equivalent manual form:

```bash
./configure --comp=g++ --dim 2 --cuda local
make -j$(nproc) bin/alamo_gpu
```

`--cuda`, `--cuda auto`, `--cuda local`, and `--cuda native` all auto-detect via
`nvidia-smi`. You can still pin an architecture explicitly, e.g. `--cuda 86`.
CUDA builds now embed PTX in addition to the native cubin; that helps smoke tests
when a binary is run on a newer compatible GPU, but benchmark numbers should use
a binary built for the exact target (`sm_80` for A100, `sm_90` for H200, `sm_86`
for the local RTX A1000).

## 2. Build (manual / non-NOVA)

CPU (baseline, with profiling):

```bash
./configure --profile
make
# -> bin/alamo-2d-<comp>
```

GPU (CUDA, with profiling). Requires the **CUDA toolkit** (`nvcc`) and an
NVIDIA GPU; `--cuda` auto-detects the local compute capability via `nvidia-smi`.
For local workstation testing prefer the helper above.
CUDA builds only compile `alamo_gpu.cc` (the Flame-only entry point) and its
actual dependency closure -- the other integrators (AllenCahn, CahnHilliard,
Dendrite, ...) pulled in by the general-purpose `alamo.cc` launcher aren't
needed for chamber/Flame runs and were never made nvcc-clean:

```bash
./configure --cuda local --profile    # or --cuda 86 to pin sm_86
make
# -> bin/alamo_gpu-2d-cuda86-<comp>
```

For correctness gates, build a no-fast-math CUDA binary:

```bash
./configure --comp=g++ --dim 2 --cuda local --cuda-fp strict
make
# -> bin/alamo_gpu-2d-nofast-cuda86-g++
```

Both binaries coexist (the POSTFIX differs).

## 3. Run with async IO + profiling

Append `benchmark/async_profile.inputs` to your chamber input, or pass the
params on the command line. Async output (`amrex.async_out=1`) hands plotfile
writes to a background thread so the GPU marches step N+1 while the CPU drains
the step-N write -- the requested compute/IO overlap, with no source changes
(`Integrator::WritePlotFile` already routes through
`amrex::WriteMultiLevelPlotfile`, which honors `async_out`).

```bash
mpiexec -np 1 ./bin/alamo_gpu-2d-cuda86-g++ input \
    amrex.async_out=1 tiny_profiler.device_synchronize_around_region=1
```

## 4. Benchmark CPU vs GPU

```bash
benchmark/benchmark_gpu_cpu.sh input 50      # prefers a CUDA binary matching this GPU
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

## Correctness gate

Before optimizing, compare CPU and GPU scalar diagnostics from the same short
Flame run:

```bash
benchmark/golden_compare_flame.sh input 1
```

The harness writes isolated CPU/GPU plot directories under
`benchmark/golden_<stamp>/` and compares `thermo.dat` with
`benchmark/compare_thermo.py`. Override binaries or tolerances with:

```bash
CPU_BIN=bin/alamo-2d-clang++ \
GPU_BIN=bin/alamo_gpu-2d-nofast-cuda86-g++ \
ABS_TOL=1e-8 REL_TOL=1e-6 \
benchmark/golden_compare_flame.sh input 1
```

Guarded CUDA-unsupported IC/BC paths should fail before entering host-loop
writes:

```bash
GPU_BIN=bin/alamo_gpu-2d-nofast-cuda86-g++ benchmark/test_gpu_guarded_ic.sh input
```

For a standardized multi-case baseline suite, record current CPU/GPU references
and later check future changes against them:

```bash
python3 benchmark/baseline_suite.py list
python3 benchmark/baseline_suite.py record
python3 benchmark/baseline_suite.py check
python3 benchmark/baseline_suite.py report
```

References are small JSON summaries under `benchmark/baseline_references/`;
full run logs and plot directories go under `benchmark/baseline_runs/`.
Set `CPU_NP`, `GPU_FAST_NP`, and `GPU_STRICT_NP` to record fair per-profile
MPI sizes. The local G0 baseline uses:

```bash
CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py record
CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py check
CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py report
```

The G0 baseline of record is summarized in `benchmark/G0_BASELINE_OF_RECORD.md`.
For the required Nsight Compute occupancy/register capture, use:

```bash
benchmark/g0_ncu_capture.sh
```

This requires NVIDIA performance counter access. If `ncu` reports
`ERR_NVGPUCTRPERM`, enable counter access on the host or run the capture on a
permitted NOVA node.

## Optimization knobs (making the GPU win big)

- **Use `--cuda-fp strict` for correctness gates.** This omits `--use_fast_math`
  and adds `--fmad=false`, producing a separate `nofast` binary.
- **Never build performance runs with `--debug`.** `--use_fast_math` (fast `exp`/`tanh`/`pow`
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

Local testing now has a separate helper (`benchmark/build_alamo_local_gpu.sh`)
that builds for the workstation GPU architecture. Production benchmark binaries
for the compute cluster should still be built natively for A100/H200 with
`benchmark/build_alamo_nova.sh`.
