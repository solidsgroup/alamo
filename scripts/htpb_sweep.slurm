#!/usr/bin/env bash
# Slurm job array: one task per (pre_exponential, activation_temperature)
# grid point, each evaluating both HTPB sandwich geometries' fit pressures
# (scripts/htpb_grid_sweep.py --mode eval). Used to find a good initial guess
# for scripts/optimize_htpb_fullfeedback_combined.py before resuming the real
# least_squares optimization, rather than its current untuned defaults.
#
# CLUSTER-SPECIFIC -- fill in before submitting:
#   - #SBATCH --partition / --account below
#   - the `module load` lines (MPI + HDF5 to match how bin/lowmach-2d-hdf5-clang++
#     was built; see FULLFEEDBACK_CALIBRATION.md's build notes for what it
#     links against -- an MPI-enabled HDF5, Eigen, libpng)
#   - LOWMACH_BIN below if the HPC build lands somewhere other than
#     $REPO_ROOT/bin/lowmach-2d-hdf5-clang++
#
# Each array task evaluates ONE grid point against BOTH geometries; within a
# geometry, run_htpb_sandwich_pressure_sweep.sh backgrounds ALL of that
# geometry's fit-pressure sims at once (up to 7, for group_100um). Each sim
# is capped at 8 cores per your instruction (MPI_NP below, default 2 ranks
# -- the value validated locally -- comfortably under that cap), but the
# job's total allocation has to cover the worst-case *concurrent* sims, not
# just one: 7 concurrent sims x MPI_NP ranks each. --ntasks below is sized
# for that (7 * MPI_NP); update it if you raise MPI_NP or the fit-pressure
# counts in optimize_htpb_fullfeedback_combined.py's DEFAULT_GEOMETRIES.
# This is still a modest, easily-schedulable per-task request (not a full
# node), just not literally 8 cores total -- 8 is the per-sim cap, not the
# per-task cap.
#
# Usage:
#   NGRID=$(python3 scripts/htpb_grid_sweep.py --mode list)
#   sbatch --array=0-$((NGRID - 1)) scripts/htpb_sweep.slurm /path/to/htpb_sweep_workdir
#
# After all array tasks complete:
#   python3 scripts/htpb_grid_sweep.py --mode aggregate --workdir /path/to/htpb_sweep_workdir

#SBATCH --job-name=htpb-grid-sweep
#SBATCH --partition=nova
#SBATCH --account=brunnels
#SBATCH --nodes=1
#SBATCH --ntasks=14
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
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

WORKDIR="${1:?usage: sbatch --array=0-N scripts/htpb_sweep.slurm <workdir>}"

export LOWMACH_BIN="${LOWMACH_BIN:-${REPO_ROOT}/bin/lowmach-2d-hdf5-g++}"
export MPI_NP="${MPI_NP:-2}"  # keep <= 8 per your per-sim core cap

python3 "${SCRIPT_DIR}/htpb_grid_sweep.py" \
    --mode eval \
    --workdir "${WORKDIR}" \
    --task-id "${SLURM_ARRAY_TASK_ID}" \
    --lowmach-bin "${LOWMACH_BIN}"
