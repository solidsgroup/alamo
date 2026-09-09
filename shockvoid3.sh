#!/usr/bin/env bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem-per-cpu=1000
#SBATCH --job-name="alamo_shock_void_growth_3"
#SBATCH --output="%x-%j-log.txt"
#SBATCH --mail-user=thoopul@iastate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
module purge
module load openmpi
module load gcc/11.4.1
module load hdf5/1.14.6-openmpi4
srun --mpi=pmix ./bin/lagrangian_shock_cp_damage-2d-hdf5-g++ ./tests/shock_cp_damage/input9
