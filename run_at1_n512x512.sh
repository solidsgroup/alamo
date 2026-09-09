#!/usr/bin/env bash
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name="alamo_shock_at1_n512x512"
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
# 512x512: 262144 cells, dt~9.4e-11 (dy-limited) -> ~128x the 256x64 baseline.
# 96h at 1 rank is very likely NOT enough. Check the 512x128 timing first.
# damage.M=0.1, growth.M=1e-1, p_thresh=4.75e9 (was 4.5e9), max_voids=20, av.CQ=2 av.CL=0.5, Gc/ell unchanged.
# diag_file -> lagshockcp_damage_diag_at1_n512x512.dat   (bare path: lands in CWD)
# xt.file   -> lagshockcp_damage_xt_at1_n512x512.dat     (bare path: lands in CWD)
# plot_file -> tests/shock_cp_damage/output_at1_n512x512
srun --mpi=pmix ./bin/lagrangian_shock_cp_damage-2d-hdf5-g++ ./tests/shock_cp_damage/input_at1_n512x512.txt
