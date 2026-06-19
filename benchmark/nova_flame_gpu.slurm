#!/usr/bin/env bash
#SBATCH --job-name=flame_gpu
#SBATCH --output=flame_gpu_%j.out
#SBATCH --error=flame_gpu_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2                 # MPI ranks; set EQUAL to total #GPUs
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:2              # GPUs per node (NOVA: a100; --gres=gpu:a100:2)
#SBATCH --gpus-per-task=1         # bind 1 GPU per rank (AMReX maps rank->GPU)
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
##### ---- EDIT THESE FOR YOUR NOVA ALLOCATION ----------------------------
##SBATCH --partition=<gpu-partition>    # e.g. `gpu`; run `sinfo` to list
##SBATCH --account=<your-slurm-account> # run `sacctmgr show assoc user=$USER`
#########################################################################
# ===========================================================================
# Multi-GPU (1 MPI rank per GPU) run of the Flame center-bore case on NOVA.
#
# PREREQUISITE: the --cuda binary must compile first. On this branch the
# Flame solver + propellant + build system are GPU-ready, but the shared
# dependency tree (Mechanics/Newton/Elastic-BC/IC) still needs porting from
# origin/gpu (benchmark/README.md, "remaining work"). Until then this job
# will fail at the build step -- use nova_flame_cpu.slurm in the meantime.
#
# Submit:   sbatch benchmark/nova_flame_gpu.slurm
# ===========================================================================
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# --- Modules (adjust to `module avail` on NOVA) ----------------------------
module purge
module load cuda/12       || module load cuda
module load gcc/12        || module load gcc
module load openmpi       || module load openmpi4   # prefer a CUDA-aware MPI

INPUT="input_nova_centre_bore"
CUDA_ARCH="${CUDA_ARCH:-80}"   # A100 = sm_80. Override: CUDA_ARCH=90 for H100.

# --- Build once if the CUDA binary is missing ------------------------------
GPU_BIN="$(ls -t bin/alamo-2d-cuda* 2>/dev/null | head -1 || true)"
if [ -z "${GPU_BIN}" ]; then
  echo "No CUDA binary found -- building (./configure --cuda ${CUDA_ARCH} --profile && make)"
  ./configure --cuda "${CUDA_ARCH}" --profile
  make -j"${SLURM_CPUS_ON_NODE:-8}"
  GPU_BIN="$(ls -t bin/alamo-2d-cuda* | head -1)"
fi
echo "GPU binary: ${GPU_BIN}  (${SLURM_NTASKS} ranks / GPUs)"

# --- Runtime overrides ------------------------------------------------------
# async_out  : background-thread IO so the GPU runs step N+1 during step N write
# managed arena : lets host diagnostics + async writer touch device field data
# tiny_profiler : attribute GPU kernel time to BL_PROFILE regions (--profile build)
OVERRIDES="amrex.async_out=1 amrex.async_out_nfiles=0 \
           amrex.the_arena_is_managed=1 \
           tiny_profiler.device_synchronize_around_region=1"

# Optional: capture an Nsight Systems trace for a GPU/CPU flame-graph timeline.
# Prefix the launch with: nsys profile --force-overwrite true -o flame_gpu_${SLURM_JOB_ID} --
NSYS_PREFIX="${NSYS_PREFIX:-}"

echo "=== launching ==="
time ${NSYS_PREFIX} srun --mpi=pmix -n "${SLURM_NTASKS}" --gpus-per-task=1 \
     "${GPU_BIN}" "${INPUT}" ${OVERRIDES}
echo "=== done ==="
