#!/usr/bin/env bash
# Runs the real optimize_htpb_fullfeedback_combined.py least_squares
# optimization as a single long-lived Slurm batch job (not the login node --
# this needs to survive for many hours/iterations).
#
# Sizing: run_htpb_sandwich_pressure_sweep.sh backgrounds ALL of a
# geometry's fit-pressure sims at once within one objective evaluation, up
# to 7 concurrent (group_100um) x MPI_NP ranks each. --ntasks below is sized
# for that (7 * MPI_NP); update it if MPI_NP or the fit-pressure counts in
# optimize_htpb_fullfeedback_combined.py's DEFAULT_GEOMETRIES change.
#
# Usage:
#   sbatch scripts/optimize_htpb_combined.slurm /path/to/workdir [extra optimize_htpb_fullfeedback_combined.py args...]

#SBATCH --job-name=htpb-combined-opt
#SBATCH --partition=nova
#SBATCH --account=brunnels
#SBATCH --nodes=1
#SBATCH --ntasks=14
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -euo pipefail

module load openmpi/4.1.6-jd5jqt4

# Slurm copies this script to a per-job spool dir before running it, so
# $(dirname "${BASH_SOURCE[0]}") would resolve there, not the repo -- use
# SLURM_SUBMIT_DIR (the cwd sbatch was invoked from) instead.
REPO_ROOT="${SLURM_SUBMIT_DIR}"
SCRIPT_DIR="${REPO_ROOT}/scripts"

source "${REPO_ROOT}/.venv_htpb_sweep/bin/activate"

WORKDIR="${1:?usage: sbatch scripts/optimize_htpb_combined.slurm <workdir> [extra args...]}"
shift

export LOWMACH_BIN="${LOWMACH_BIN:-${REPO_ROOT}/bin/lowmach-2d-hdf5-g++}"
export MPI_NP="${MPI_NP:-2}"  # keep <= 8 per your per-sim core cap

python3 "${SCRIPT_DIR}/optimize_htpb_fullfeedback_combined.py" \
    --workdir "${WORKDIR}" \
    --lowmach-bin "${LOWMACH_BIN}" \
    "$@"
