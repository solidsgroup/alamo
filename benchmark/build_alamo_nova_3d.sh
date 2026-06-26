#!/bin/bash
# ============================================================================
# build_alamo_nova_3d.sh  --  download + build the GPU (chamber-gpu) Alamo on NOVA (3D)
#
# Run this ON THE NOVA LOGIN NODE:   sh build_alamo_nova_3d.sh
#
# It (1) clones/updates the chamber-gpu branch, (2) primes the AMReX checkout on
# the login node (the only step that needs internet), then (3) submits a CPU
# build job that compiles the currently supported NOVA GPU targets.
# nvcc cross-compiles for every arch below, so the build needs no GPU -- it runs on the
# fast EPYC nodes. When the job finishes you'll have:
#     bin/alamo_gpu-3d-profile-cuda70-g++   (V100)
#     bin/alamo_gpu-3d-profile-cuda80-g++   (A100)
#     bin/alamo_gpu-3d-profile-cuda90-g++   (H200)
#     bin/alamo-3d-profile-g++              (CPU baseline for the crossover)
# All are --profile builds: fully optimized (--use_fast_math etc. on GPU) AND able
# to emit TinyProfiler tables when you pass the profiler runtime params. The CPU
# baseline is built here too so nova_flame_cpu_3d.slurm never has to build in-job.
# Override arches with e.g. ARCHES=80 (A100 only) for a faster build.
# This is the 3D sibling of build_alamo_nova.sh (configure uses --dim 3).
# ============================================================================
set -euo pipefail

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'

REPO_URL="${REPO_URL:-https://github.com/solidsgroup/alamo.git}"
BRANCH="${BRANCH:-chamber-gpu}"
ALAMO_DIR="${ALAMO_DIR:-$PWD}"
ACCOUNT="${ACCOUNT:-brunnels}"
BUILD_PARTITION="${BUILD_PARTITION:-nova}"
ARCHES="${ARCHES:-70 80 90}"                    # 70=V100, 80=A100, 90=H200
BUILD_JOBS="${BUILD_JOBS:-64}"
CPU_BUILD_JOBS="${CPU_BUILD_JOBS:-8}"           # bounded -j for the CPU (-flto) baseline; a high -j OOMs the node
COMP="${COMP:-g++}"
EMAIL="${EMAIL:-jackplum@iastate.edu}"

echo -e "${BOLD}${BLUE}=== Alamo GPU build setup (NOVA, 3D) ===${NC}"
echo -e "  repo=${REPO_URL}  branch=${BRANCH}"
echo -e "  dir=${ALAMO_DIR}  archs='${ARCHES}'  account=${ACCOUNT}  build-partition=${BUILD_PARTITION}"

echo -e "${YELLOW}Loading modules...${NC}"
module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || echo -e "${RED}  WARN: load a cuda module manually${NC}"
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || echo -e "${RED}  WARN: load a gcc module manually${NC}"
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || echo -e "${RED}  WARN: load an openmpi module manually${NC}"
module list 2>&1 | sed 's/^/    /' || true

# Fail loudly if CUDA didn't load. Without nvcc on PATH, ./configure --cuda
# silently falls back to a CPU-only build (wrong binaries, and the -flto CPU
# path then OOMs at high -j). Refuse rather than waste a build job.
if ! command -v nvcc >/dev/null 2>&1; then
  echo -e "${RED}ERROR: nvcc not found after module load.${NC}" >&2
  echo -e "${RED}  Load a CUDA module first (e.g. 'module avail cuda' then 'module load cuda'), then re-run.${NC}" >&2
  exit 1
fi
echo -e "${GREEN}  nvcc: $(command -v nvcc)${NC}"

mkdir -p "${ALAMO_DIR}"
ALAMO_DIR="$(cd "${ALAMO_DIR}" && pwd)"
echo -e "${YELLOW}Building in ${ALAMO_DIR}${NC}"
if [ -d "${ALAMO_DIR}/.git" ]; then
  echo -e "${YELLOW}Updating existing checkout in place...${NC}"
  git -C "${ALAMO_DIR}" fetch origin "${BRANCH}"
  git -C "${ALAMO_DIR}" checkout "${BRANCH}"
  git -C "${ALAMO_DIR}" pull --ff-only origin "${BRANCH}"
elif [ -f "${ALAMO_DIR}/configure" ] && [ -d "${ALAMO_DIR}/src" ]; then
  echo -e "${GREEN}  existing Alamo source found -- building in place${NC}"
else
  echo -e "${YELLOW}Fetching ${BRANCH} into ${ALAMO_DIR}...${NC}"
  git -C "${ALAMO_DIR}" init -q
  git -C "${ALAMO_DIR}" remote add origin "${REPO_URL}" 2>/dev/null     || git -C "${ALAMO_DIR}" remote set-url origin "${REPO_URL}"
  git -C "${ALAMO_DIR}" fetch --depth 1 origin "${BRANCH}"
  git -C "${ALAMO_DIR}" checkout -b "${BRANCH}" FETCH_HEAD
fi
cd "${ALAMO_DIR}"

echo -e "${YELLOW}Priming AMReX checkout (login node)...${NC}"
./configure --comp="${COMP}" --dim 3 --cuda "$(echo ${ARCHES} | awk '{print $1}')" --profile --get-eigen >/tmp/alamo_prime.log 2>&1 || {
  echo -e "${RED}configure prime failed -- see /tmp/alamo_prime.log${NC}"; tail -20 /tmp/alamo_prime.log; exit 1; }
echo -e "${GREEN}  AMReX present under ext/${NC}"

SB="${ALAMO_DIR}/build_alamo_gpu_3d.sbatch"
cat > "${SB}" <<END_OF_SBATCH
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -J alamo_build_3d
#SBATCH -D ${ALAMO_DIR}
#SBATCH --partition=${BUILD_PARTITION}
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${BUILD_JOBS}
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=alamo_build_3d.%j.out
#SBATCH --error=alamo_build_3d.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${EMAIL}
set -euo pipefail
cd ${ALAMO_DIR}
module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true
if ! command -v nvcc >/dev/null 2>&1; then
    echo "ERROR: nvcc not found in build job after module load -- refusing to silently build a CPU-only binary. Fix the cuda module name and resubmit." >&2
    exit 1
fi
# GPU (CUDA) binaries. The nvcc path is -fno-lto, so a high -j is safe.
for arch in ${ARCHES}; do
    echo "=== building cuda sm_\${arch} (3D) ==="
    ./configure --comp=${COMP} --dim 3 --cuda \${arch} --profile --get-eigen
    make -j\${SLURM_CPUS_PER_TASK:-${BUILD_JOBS}} bin/alamo_gpu
done
# CPU baseline binary for the GPU-vs-CPU crossover (what nova_flame_cpu_3d.slurm
# launches). This is the -flto path, so build ONLY bin/alamo with a bounded -j;
# a high -j here OOMs the node -> "lto1: error: ...file too short".
echo "=== building CPU baseline bin/alamo (3D, -flto, -j${CPU_BUILD_JOBS}) ==="
./configure --comp=${COMP} --dim 3 --profile --get-eigen
make -j${CPU_BUILD_JOBS} bin/alamo
echo "=== build complete ==="
ls -lh bin/alamo_gpu-3d*cuda* bin/alamo-3d*-g++ || true
END_OF_SBATCH

echo -e "${BOLD}${GREEN}Submitting build job...${NC}"
JOBID=$(sbatch --parsable "${SB}")
echo -e "${GREEN}  submitted job ${JOBID} on partition ${BUILD_PARTITION} (account ${ACCOUNT})${NC}"
echo ""
echo -e "${BOLD}Next:${NC}"
echo -e "  watch:   squeue -j ${JOBID}"
echo -e "  log:     tail -f ${ALAMO_DIR}/alamo_build_3d.${JOBID}.out"
echo -e "  when done, the binaries are in ${ALAMO_DIR}/bin/ :"
echo -e "     alamo_gpu-3d-profile-cuda70-g++   (V100)"
echo -e "     alamo_gpu-3d-profile-cuda80-g++   (A100)"
echo -e "     alamo_gpu-3d-profile-cuda90-g++   (H200)"
echo -e "     alamo-3d-profile-g++              (CPU baseline)"
echo -e "  then run:  sbatch ${ALAMO_DIR}/benchmark/nova_flame_gpu_3d.slurm   (edit GPU_TYPE first)"
