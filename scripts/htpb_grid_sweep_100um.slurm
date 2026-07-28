#!/usr/bin/env bash
# Slurm job array: one task per (pre_exponential, activation_temperature)
# grid point, each evaluating the 100um group geometry's cheap min/max
# pressure endpoints only (scripts/htpb_grid_sweep.py --mode eval, default
# endpoints-only). Meant to feed the local GP-based Bayesian optimization
# (scripts/optimize_htpb_fullfeedback_combined.py) with extra prior sample
# points, not to replace it.
#
# Sized for the reworked 100um template: halved/symmetric domain (n_cell.x=16
# instead of 32), rocfire thermo/transport (tpg EOS), reflect_even/odd x BCs.
# Each array task gets 16 MPI ranks total, split across however many
# pressures run concurrently (2 endpoints here -> 8 ranks/sim) via
# scripts/optimize_htpb_fullfeedback_combined.py's cores-based MPI_NP sizing.
#
# --time is generous (6h) so slow/near-crash parameter points still finish
# or hit their chemistry-integration abort cleanly, rather than timing out
# with no result the way the earlier wider/pricier sweep did.
#
# Usage:
#   NGRID=$(python3 scripts/htpb_grid_sweep.py --mode list)
#   sbatch --array=0-$((NGRID - 1)) scripts/htpb_grid_sweep_100um.slurm /path/to/workdir
#
# After all array tasks complete:
#   python3 scripts/htpb_grid_sweep.py --mode aggregate --workdir /path/to/workdir

#SBATCH --job-name=htpb-grid-100um
#SBATCH --partition=nova
#SBATCH --account=brunnels
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -euo pipefail

module load openmpi/4.1.6-jd5jqt4

# Slurm copies this script to a per-job spool directory before running it,
# so $(dirname "${BASH_SOURCE[0]}") would resolve there, not to the repo --
# use SLURM_SUBMIT_DIR (the cwd sbatch was invoked from) instead.
REPO_ROOT="${SLURM_SUBMIT_DIR}"
SCRIPT_DIR="${REPO_ROOT}/scripts"

source "${REPO_ROOT}/.venv_htpb_sweep/bin/activate"

WORKDIR="${1:?usage: sbatch --array=0-N scripts/htpb_grid_sweep_100um.slurm <workdir>}"

export LOWMACH_BIN="${LOWMACH_BIN:-${REPO_ROOT}/bin/lowmach-2d-hdf5-g++}"

python3 "${SCRIPT_DIR}/htpb_grid_sweep.py" \
    --mode eval \
    --workdir "${WORKDIR}" \
    --task-id "${SLURM_ARRAY_TASK_ID}" \
    --lowmach-bin "${LOWMACH_BIN}" \
    --cores 16
