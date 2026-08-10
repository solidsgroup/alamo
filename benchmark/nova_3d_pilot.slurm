#!/usr/bin/env bash
#SBATCH -A brunnels
#SBATCH -J pilot3d
#SBATCH --partition=nova
#SBATCH --time=00:20:00
#SBATCH --output=pilot3d.%j.out
#SBATCH --error=pilot3d.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
#
# 3D pilot cost-per-step measurement. One deck, one hardware config per job.
# Resources (nodes/ntasks/mem/gres) are set on the sbatch command line by
# benchmark/submit_nova_3d_pilot.sh -- see the matrix there.
#
# Env:
#   INPUT    (required)  input_nova_3d_pilot_{rod_and_tube,cross}
#   BACKEND  cpu|gpu     (default cpu)
#   ARCH     80|90       CUDA arch for GPU binary glob (a100=80, h200=90)
#   TAG      label appended to the output dir (e.g. cpu64, a100x2)
#
# Build first (once, on a login/build node):
#   CPU: ./configure --comp=g++ --dim=3 --profile && make -j8
#   GPU: benchmark/build_alamo_nova_3d.sh   (cuda80 + cuda90)
#
# The 0.05 s deck finishes if it can; otherwise the 20-min wall kills it.
# Either way cost/step comes from the step log + thermo.dat + tiny profiler.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

INPUT="${INPUT:?set INPUT to one input_nova_3d_pilot_* deck}"
BACKEND="${BACKEND:-cpu}"
ARCH="${ARCH:-80}"
TAG="${TAG:-${BACKEND}}"
TIME_CMD="${TIME_CMD:-/usr/bin/time -v}"

module purge 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true
if [[ "${BACKEND}" == "gpu" ]]; then
    module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
fi

if [[ "${BACKEND}" == "gpu" ]]; then
    BIN="$(ls -t bin/alamo_gpu-3d*cuda${ARCH}*-g++ 2>/dev/null | head -1 || true)"
    [[ -n "${BIN}" ]] || { echo "No bin/alamo_gpu-3d*cuda${ARCH}* -- run benchmark/build_alamo_nova_3d.sh" >&2; exit 1; }
    # Bigger boxes for GPU (fewer, fatter kernels); pin arena so init never
    # grabs 3/4 of host RAM under managed memory.
    OVERRIDES=(
        "amr.max_grid_size=64"
        "amrex.async_out=0" "amrex.async_out_nfiles=0"
        "amrex.the_arena_is_managed=0"
    )
    LAUNCH=(srun --mpi=pmix -n "${SLURM_NTASKS}" --gpus-per-task=1)
else
    BIN="$(ls -t bin/alamo-3d*-g++ 2>/dev/null | grep -v cuda | head -1 || true)"
    [[ -n "${BIN}" ]] || { echo "No 3D CPU binary -- ./configure --comp=g++ --dim=3 --profile && make -j8" >&2; exit 1; }
    OVERRIDES=(
        "amrex.async_out=0" "amrex.async_out_nfiles=0"
        "tiny_profiler.device_synchronize_around_region=1"
    )
    LAUNCH=(srun --mpi=pmix -n "${SLURM_NTASKS}")
fi

OUTDIR="$(grep -E '^plot_file' "${INPUT}" | awk -F'= *' '{print $2}').${TAG}_${SLURM_JOB_ID}"
OVERRIDES+=("plot_file=${OUTDIR}")
# EXTRA: space-separated additional deck overrides, e.g.
#   EXTRA="elastic.solver.bottom_solver=smoother"
if [[ -n "${EXTRA:-}" ]]; then
    read -r -a EXTRA_ARR <<< "${EXTRA}"
    OVERRIDES+=("${EXTRA_ARR[@]}")
fi

echo "=== 3D pilot ==="
echo "input:   ${INPUT}"
echo "binary:  ${BIN}"
echo "backend: ${BACKEND} tag=${TAG} ranks=${SLURM_NTASKS}"
echo "outdir:  ${OUTDIR}"
echo "git:     $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

rc=0
${TIME_CMD} "${LAUNCH[@]}" "${BIN}" "${INPUT}" "${OVERRIDES[@]}" || rc=$?
echo "=== done rc=${rc} (walltime cap 20 min; nonzero after ~19 min = killed mid-run, still usable) ==="
