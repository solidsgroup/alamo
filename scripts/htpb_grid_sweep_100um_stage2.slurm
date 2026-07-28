#!/usr/bin/env bash
# Slurm job array: one task per (pre_exponential, activation_temperature)
# grid point, each evaluating the 100um group geometry against EVERY fit
# pressure (0.8-4.0 MPa, 7 points) -- i.e. the same objective as
# optimize_htpb_fullfeedback_combined.py's stage 2, not the cheap endpoints
# used by htpb_grid_sweep_100um.slurm. Meant to seed stage 2's GP with real
# full-pressure data before it runs any expensive local evaluations (see
# --stage2-prior-dir on optimize_htpb_fullfeedback_combined.py).
#
# Grid is intentionally narrow (unlike the first, exploratory endpoints-only
# sweep): scoped around the region the endpoints sweep and local stage 1
# have already shown is numerically stable and low-residual
# (pre_exponential ~1e-4 to 1.5e-3, activation_temperature ~3000-5000 K),
# rather than re-covering the wide range that mostly just finds crashes.
#
# Each array task gets 32 MPI ranks total, split across the 7 concurrent
# fit-pressure sims (~4-5 ranks/sim). --time is generous (8h) since
# full-pressure evaluations include the slower high-pressure endpoints.
#
# Usage:
#   sbatch --array=0-34 scripts/htpb_grid_sweep_100um_stage2.slurm /path/to/workdir

#SBATCH --job-name=htpb-grid-100um-s2
#SBATCH --partition=nova
#SBATCH --account=brunnels
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

module load openmpi/4.1.6-jd5jqt4

REPO_ROOT="${SLURM_SUBMIT_DIR}"
SCRIPT_DIR="${REPO_ROOT}/scripts"

source "${REPO_ROOT}/.venv_htpb_sweep/bin/activate"

WORKDIR="${1:?usage: sbatch --array=0-34 scripts/htpb_grid_sweep_100um_stage2.slurm <workdir>}"

export LOWMACH_BIN="${LOWMACH_BIN:-${REPO_ROOT}/bin/lowmach-2d-hdf5-g++}"

python3 "${SCRIPT_DIR}/htpb_grid_sweep.py" \
    --mode eval \
    --workdir "${WORKDIR}" \
    --task-id "${SLURM_ARRAY_TASK_ID}" \
    --lowmach-bin "${LOWMACH_BIN}" \
    --cores 32 \
    --all-pressures \
    --pre-exponential-grid 0.0001 0.0002 0.0003 0.0005 0.0007 0.001 0.0015 \
    --activation-temperature-grid 3000 3500 4000 4500 5000
