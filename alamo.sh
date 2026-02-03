#!/usr/bin/env bash
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=15000
#SBATCH --job-name="alamo"
#SBATCH --output="%x-%j-log.txt"
#SBATCH --mail-user=mungerct@iastate.edu
#SBATCH --mail-type=BEGIN,END,FAIL

echo "======================================================"
echo " Job '$SLURM_JOB_NAME' (ID: $SLURM_JOB_ID) is starting"
echo " Submitted by: $SLURM_JOB_USER"
echo " Running on node(s): $SLURM_NODELIST"
echo " Start time: $(date)"
echo "======================================================"

module load openmpi_hpc
# srun --mpi=pmix ./bin/alamo-2d-g++ input.scpspheres_kodga_validation phi.ic.psread.file.name="setC_xyzrs/uni_R240um_AP6499_10.xyzr" plot_file="output.scpsphereselastic_set10_phi_eps_8"
srun --mpi=pmix ./bin/alamo-2d-g++ input.scpspheres_kodga_validation
# srun --mpi=pmix ./bin/alamo-2d-g++ input.scpspheres_arbitary_geometry
