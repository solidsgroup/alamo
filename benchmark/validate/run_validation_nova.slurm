#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J validate_phase1
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --output=validate_phase1.%j.out
#SBATCH --error=validate_phase1.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
# ============================================================================
# run_validation_nova.slurm -- Phase 1 task 1.B, NOVA entry point.
#
# NOVA-only: separate from run_validation_local.sh by design (see
# benchmark/validate/README.md, "Two entry points, not one auto-detecting
# script") -- NOVA's Slurm/module environment doesn't share a launch path
# with the local mpiexec runner. Conventions here (account, mail-user,
# GPU_TYPE->arch map, module loads, srun --mpi=pmix, ASYNC_OUT=0) follow
# benchmark/NOVA_SLURM_RUNBOOK.md and the existing nova_flame_*.slurm scripts.
#
# A single submission runs ONE profile (resource request differs: GPU profiles
# need --gres, cpu doesn't) -- mirrors the existing nova_flame_cpu_3d.slurm /
# nova_flame_gpu_3d.slurm split. Submit once per profile you want a bundle for.
#
# Build the binaries first: benchmark/build_alamo_nova.sh (2D) and/or
# benchmark/build_alamo_nova_3d.sh (3D). Those scripts only produce `--profile`
# (fast-math) targets today -- see run_validation_nova.py's module docstring
# for the gap if you ask for PROFILE=gpu_strict before a strict NOVA binary
# exists.
#
# Examples:
#   sbatch benchmark/validate/run_validation_nova.slurm                       # PROFILE=cpu default, see below
#   PROFILE=gpu_fast GPU_TYPE=a100 sbatch --gres=gpu:a100:1 --ntasks=1 \
#       --ntasks-per-node=1 --cpus-per-task=8 --mem=32G \
#       benchmark/validate/run_validation_nova.slurm
#   PROFILE=cpu sbatch --ntasks=64 --cpus-per-task=1 --mem=64G \
#       benchmark/validate/run_validation_nova.slurm
#   CASE=centre_bore_3d_256_a2_tuned PROFILE=gpu_fast GPU_TYPE=h200 \
#       sbatch --gres=gpu:h200:1 --ntasks=1 --cpus-per-task=8 --mem=64G \
#       benchmark/validate/run_validation_nova.slurm
#
# Output bundles land under benchmark/validate/runs/ on NOVA's own
# filesystem (e.g. /work/brunnels/$USER/alamo/benchmark/validate/runs/) --
# copy them back to compare against local bundles per
# benchmark/NOVA_SLURM_RUNBOOK.md's file-transfer guidance (this script does
# not attempt to sync across machines itself).
# ============================================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

PROFILE="${PROFILE:-cpu}"            # cpu | gpu_strict | gpu_fast
GPU_TYPE="${GPU_TYPE:-a100}"         # v100 (sm_70), a100 (sm_80), h200 (sm_90)
CASE="${CASE:-}"                     # optional: restrict to one case id

module purge 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true
if [ "${PROFILE}" != "cpu" ]; then
  module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
fi
module load python 2>/dev/null || module load python3 2>/dev/null || true
module list 2>&1 | sed 's/^/    /' || true

python3 -c "import yaml, numpy" 2>/dev/null || {
  echo "ERROR: python3 needs pyyaml + numpy (yt is also needed unless every case is" >&2
  echo "  run with extract_metrics.py --skip-field-norms). Load/activate the right" >&2
  echo "  python environment before submitting, or extend this script's module loads." >&2
  exit 1
}

echo "=== run_validation_nova: profile=${PROFILE} gpu_type=${GPU_TYPE} case=${CASE:-<all nova cases>} ==="
echo "SLURM_NTASKS=${SLURM_NTASKS:-unset} SLURM_JOB_ID=${SLURM_JOB_ID:-unset}"

ARGS=(--profile "${PROFILE}" --gpu-type "${GPU_TYPE}")
[ -n "${CASE}" ] && ARGS+=(--case "${CASE}")

python3 benchmark/validate/run_validation_nova.py "${ARGS[@]}"
echo "=== done ==="
