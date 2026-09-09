#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000
#SBATCH --job-name="alamo_shock_at1_dM0p0"
#SBATCH --output="%x-%j-log.txt"
#SBATCH --mail-user=thoopul@iastate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
module purge
module load openmpi
module load gcc/11.4.1
module load hdf5/1.14.6-openmpi4

### OLD: srun --mpi=pmix ./bin/lagrangian_shock_cp_damage-2d-hdf5-g++ ./tests/shock_cp_damage/input6
# damage.M = 0.0, growth.M = 1e-2, nucleation.p_thresh = 4.7e9
# diag_file -> lagshockcp_damage_diag_at1_dM0p0.dat   (bare path: lands in CWD)
# xt.file   -> lagshockcp_damage_xt_at1_dM0p0.dat     (bare path: lands in CWD)
# plot_file -> tests/shock_cp_damage/output_at1_dM0p0
srun --mpi=pmix ./bin/lagrangian_shock_cp_damage-2d-hdf5-g++ ./tests/shock_cp_damage/input_at1_dM0p0.txt
