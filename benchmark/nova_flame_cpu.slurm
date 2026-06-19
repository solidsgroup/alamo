#!/usr/bin/env bash
#SBATCH --job-name=flame_cpu
#SBATCH --output=flame_cpu_%j.out
#SBATCH --error=flame_cpu_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4                 # MPI ranks = domain decomposition
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
##### ---- EDIT THESE FOR YOUR NOVA ALLOCATION ----------------------------
##SBATCH --partition=<cpu-partition>   # e.g. `nova`; run `sinfo` to list
##SBATCH --account=<your-slurm-account> # run `sacctmgr show assoc user=$USER`
#########################################################################
# ===========================================================================
# CPU run of the Flame (chamber) center-bore case on ISU's NOVA cluster.
# Works with the binary that builds today (the GPU build needs the dependency
# tree ported first -- see benchmark/README.md).
#
# Submit:   sbatch benchmark/nova_flame_cpu.slurm
# ===========================================================================
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# --- Modules (adjust to `module avail` on NOVA) ----------------------------
module purge
module load gcc/12        || module load gcc
module load openmpi       || module load openmpi4
# Optional, only if you built with --profile:
# module load <anything the profiler needs>

INPUT="input_nova_centre_bore"

# --- Build once if the CPU binary is missing -------------------------------
CPU_BIN="$(ls -t bin/alamo-2d-*-* 2>/dev/null | grep -v cuda | head -1 || true)"
if [ -z "${CPU_BIN}" ]; then
  echo "No CPU binary found -- building (./configure --profile && make)"
  ./configure --profile
  make -j"${SLURM_CPUS_ON_NODE:-8}"
  CPU_BIN="$(ls -t bin/alamo-2d-*-* | grep -v cuda | head -1)"
fi
echo "CPU binary: ${CPU_BIN}"

# --- Runtime overrides ------------------------------------------------------
# async_out gives compute/IO overlap; NOVA's MPI provides MPI_THREAD_MULTIPLE.
# tiny_profiler.* is only consumed by a --profile build; drop it otherwise.
OVERRIDES="amrex.async_out=1 amrex.async_out_nfiles=0 \
           tiny_profiler.device_synchronize_around_region=1"

echo "=== launching ${SLURM_NTASKS} MPI rank(s) ==="
time srun --mpi=pmix -n "${SLURM_NTASKS}" "${CPU_BIN}" "${INPUT}" ${OVERRIDES}
echo "=== done ==="
