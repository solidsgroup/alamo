#!/usr/bin/env bash
#SBATCH --account=brunnels
#SBATCH --partition=nova
#SBATCH --job-name=fapply_a100_build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=nova-build.%j.out
#SBATCH --error=nova-build.%j.err

set -euo pipefail

REPO=/work/brunnels/jackplum/alamo-fapply-20260722
EXPECTED_HEAD=2e6a8f8f5430a58dd6a85b07b08d472e2bd1d5cd
EXPECTED_AMREX=bf15fdce52093c715a1dd9a9da64ba66f0233be1
BIN="$REPO/bin/alamo_gpu-3d-profile-cuda80-g++"

cd "$REPO"
test "$(git rev-parse HEAD)" = "$EXPECTED_HEAD"
test "$(git -C ext/AMReX-Codes/amrex rev-parse HEAD)" = "$EXPECTED_AMREX"
git diff --exit-code -- src/Operator/Elastic.cpp src/Operator/Elastic.H \
  src/Set/Matrix4_Major.H
test "$(sha256sum ext/AMReX-Codes/amrex/Src/LinearSolvers/MLMG/AMReX_MLCGSolver.H | awk '{print $1}')" = \
  9b96d56e2b97d0181f4d36c4033d29af9b9ee77df4fe3f2977c7ca2b6e9c9ee9
test "$(sha256sum ext/AMReX-Codes/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG.H | awk '{print $1}')" = \
  3b04e11295779d17a0b1e1dc23e6a2305474c4a8ed4997ac13b9a46e3c36b842

module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12
module load gcc 2>/dev/null || module load gcc/12
module load openmpi 2>/dev/null || module load openmpi4
module -t list
nvcc --version
mpicxx --version | head -1

./configure --comp=g++ --dim 3 --cuda 80 --cuda-fp fast --profile \
  --get-eigen --offline
make -j"${SLURM_CPUS_PER_TASK}" bin/alamo_gpu

test -x "$BIN"
sha256sum "$BIN" input_3d_centre_bore_128_a2 src/Operator/Elastic.cpp \
  ext/AMReX-Codes/amrex/Src/LinearSolvers/MLMG/AMReX_MLCGSolver.H \
  ext/AMReX-Codes/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG.H
echo NOVA_A100_BUILD_PASS
