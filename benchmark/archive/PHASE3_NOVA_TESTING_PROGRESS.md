# Phase 3 NOVA Testing Progress

**Date:** 2026-06-20  
**Branch:** `chamber-gpu`  
**Purpose:** Track Phase 3 NOVA one-shot testing, observed failures, fixes, and the current restart point.

## Current Status

The Phase 3 3D GPU path is locally validated and the NOVA one-shot workflow has been iterated through two submitted bundles. All failures observed so far have been in the NOVA launch/workflow layer, not in the local 3D CUDA build or 3D single-GPU smoke path.

Latest pushed fix:

```text
23c924721 Disable async output for NOVA Phase 3 runs
```

Before the next NOVA run, the checkout on NOVA should be updated to this commit or newer:

```bash
cd ~/alamo
git pull origin chamber-gpu
git log --oneline -3
bash benchmark/phase3_nova_oneshot.sh
```

## Local Pre-NOVA Validation

The local 3D GPU branch compiled and completed a 100-step run before NOVA testing:

```bash
source benchmark/local_cuda_env.sh
DIM=3 PROFILE=1 SMOKE=0 ARCH=86 ./benchmark/build_alamo_local_gpu.sh
mpiexec -np 1 ./bin/alamo_gpu-3d-profile-cuda86-g++ input_3d_flame_128 \
  max_step=100 stop_time=1e-2_s \
  plot_file=output_local_cuda_3d_100 \
  amr.plot_int=-1 amr.thermo.plot_int=-1
```

Result:

```text
STEP 100 ends. TIME = 0.01
```

Peak local GPU arena use was approximately 5.9 GB on the RTX A1000.

## NOVA Bundle 1: Initial One-Shot Failure

Bundle:

```text
phase3_nova_bundle_20260620_194349.tar.gz
```

Build result:

```text
alamo_build_3d.11159969.out
DONE bin/alamo_gpu-3d-profile-cuda80-g++
DONE bin/alamo_gpu-3d-profile-cuda90-g++
```

Observed GPU failure:

```text
srun: error: Unable to create step for job <id>: Invalid generic resource (gres) specification
```

Cause:

The sweep submitted GPU jobs with dependencies but without explicit GPU resources. The scripts used `srun --gpus-per-task=1`, but the enclosing `sbatch` allocation did not request `--gres`.

Observed CPU failure:

```text
fatal error: eigen3/Eigen/Core: No such file or directory
fatal error: eigen3/Eigen/Eigenvalues: No such file or directory
```

Additional CPU issue:

Concurrent CPU baseline jobs all attempted to configure/build in the same checkout, which caused AMReX checkout/build races such as `.git/index.lock`.

Fix pushed:

```text
d599bc99a Fix Phase 3 NOVA one-shot resources
```

Changes:

- `benchmark/phase3_scaling_sweep.sh` now submits GPU jobs with explicit `--nodes`, `--gres=gpu:<type>:<count>`, `--ntasks`, and `--ntasks-per-node`.
- `8` GPU jobs default to `2` nodes x `4` GPUs via `GPUS_PER_NODE=4`.
- `benchmark/phase3_nova_oneshot.sh` submits one serialized CPU build job.
- CPU baseline jobs depend on that CPU build job.
- CPU build uses `./configure --get-eigen`.
- CPU baseline run jobs use `BUILD_CPU_IF_MISSING=0` under the one-shot so they fail fast instead of rebuilding concurrently.

## NOVA Bundle 2: Async Output Runtime Abort

Bundle:

```text
phase3_nova_bundle_20260620_200739.tar.gz
```

Checkout:

```text
d599bc99ae9a2242c76622a0f6779caa220fee02
```

Job graph status:

The resource/dependency fix worked. `commands.txt` contained GPU submissions shaped like:

```text
sbatch --nodes=1 --gres=gpu:a100:1 --ntasks=1 --ntasks-per-node=1 ...
sbatch --nodes=1 --gres=gpu:a100:2 --ntasks=2 --ntasks-per-node=2 ...
sbatch --nodes=1 --gres=gpu:a100:4 --ntasks=4 --ntasks-per-node=4 ...
sbatch --nodes=2 --gres=gpu:a100:4 --ntasks=8 --ntasks-per-node=4 ...
```

CPU baselines depended on the serialized CPU build job:

```text
BUILD_CPU_IF_MISSING=0 sbatch --dependency=afterok:11159999 benchmark/nova_flame_cpu_3d.slurm
```

Observed runtime failure across captured GPU and CPU jobs:

```text
amrex::Abort::0::AsyncOut with 0 and 1 processes requires MPI_THREAD_MULTIPLE at runtime, but got MPI_THREAD_SINGLE !!!
```

The same failure appeared for 1 GPU, 2 GPU, and 16-rank CPU baseline jobs. The jobs reached AMReX/MPI/CUDA initialization, then aborted before stepping:

```text
Initializing AMReX (26.06)...
MPI initialized with thread support level 0
Initializing CUDA...
CUDA initialized with 1 device.
```

Cause:

The NOVA `srun --mpi=pmix` launch reports `MPI_THREAD_SINGLE`, while the Phase 3 run scripts forced:

```text
amrex.async_out=1
```

AMReX async output requires `MPI_THREAD_MULTIPLE`, so the job abort is expected in this launch environment.

Fix pushed:

```text
23c924721 Disable async output for NOVA Phase 3 runs
```

Changes:

- `benchmark/nova_flame_gpu_3d.slurm`
- `benchmark/nova_flame_gpu_3d_multi.slurm`
- `benchmark/nova_flame_cpu_3d.slurm`

now default:

```text
ASYNC_OUT=0
amrex.async_out=${ASYNC_OUT}
```

`ASYNC_OUT=1` remains available as an explicit override if NOVA is later launched/configured with `MPI_THREAD_MULTIPLE`.

## Current Known Fixed Issues

- Missing GPU `gres` allocation.
- CPU baseline build race in shared checkout.
- Missing Eigen during NOVA CPU build.
- 8-GPU allocation needing two 4-GPU nodes.
- AMReX async output abort under NOVA `MPI_THREAD_SINGLE`.

## Remaining Risks To Validate

These have not been ruled out yet because jobs have not completed stepping on NOVA after the async-output fix:

- Memory pressure for the 512 case on 1 GPU or CPU node.
- Multi-GPU domain decomposition behavior at 4 and 8 GPUs.
- Walltime or queue/node-specific failures.
- Runtime performance variation from disabling async output.
- Whether H200 runs need a different `GPU_TYPE=h200` and queue/resource shape.

## Next NOVA Run

Cancel any old Phase 3 jobs submitted before commit `23c924721`, because Slurm may have captured the old script contents.

Then:

```bash
cd ~/alamo
git pull origin chamber-gpu
git log --oneline -3
bash benchmark/phase3_nova_oneshot.sh
```

Expected latest commits:

```text
23c924721 Disable async output for NOVA Phase 3 runs
d599bc99a Fix Phase 3 NOVA one-shot resources
d850ce22d Stage Phase 3 NOVA one-shot launcher
```

## Data To Capture After Next Run

Capture the same bundle inputs:

- `commands.txt`
- `alamo_build_3d*.out`
- `alamo_build_3d*.err`
- `alamo_cpu_build_3d*.out`
- `alamo_cpu_build_3d*.err`
- `flame_gpu_3d*.out`
- `flame_gpu_3d*.err`
- `flame_gpu_3d_multi*.out`
- `flame_gpu_3d_multi*.err`
- `flame_cpu_3d*.out`
- `flame_cpu_3d*.err`
- `benchmark/phase3_nova_oneshot.sh`
- `benchmark/phase3_scaling_sweep.sh`
- `benchmark/nova_flame_gpu_3d.slurm`
- `benchmark/nova_flame_gpu_3d_multi.slurm`
- `benchmark/nova_flame_cpu_3d.slurm`

Key success markers to look for:

```text
STEP <n> ends.
TinyProfiler total time / region table
=== done ===
```

Key failure markers to search:

```text
amrex::Abort
CUDA
Out of memory
srun: error
CANCELLED
TIME LIMIT
SIGABRT
SIGKILL
```
