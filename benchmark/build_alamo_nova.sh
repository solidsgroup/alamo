#!/bin/bash
# ============================================================================
# build_alamo_nova.sh  --  download + build the GPU (chamber-gpu) Alamo on NOVA
#
# Run this ON THE NOVA LOGIN NODE:   sh build_alamo_nova.sh
#
# It (1) clones/updates the chamber-gpu branch, (2) primes the AMReX checkout on
# the login node (the only step that needs internet), then (3) submits a CPU
# build job that compiles BOTH the A100 (sm_80) and H200 (sm_90) binaries.
# nvcc cross-compiles for both archs, so the build needs no GPU -- it runs on the
# fast EPYC nodes. When the job finishes you'll have:
#     bin/alamo_gpu-2d-profile-cuda80-g++   (A100)
#     bin/alamo_gpu-2d-profile-cuda90-g++   (H200)
# Both are --profile builds: fully optimized (--use_fast_math etc.) AND able to
# emit TinyProfiler tables when you pass the profiler runtime params.
# ============================================================================
set -euo pipefail

# Colors
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'

# ---- Config (override via env, e.g. ALAMO_DIR=/some/path sh build_alamo_nova.sh) ----
# Defaults to HTTPS (no SSH key needed). If the repo is private, cache a GitHub
# token first (e.g. `git config --global credential.helper store` then one auth'd
# clone) or set REPO_URL to the SSH form.
REPO_URL="${REPO_URL:-https://github.com/solidsgroup/alamo.git}"
BRANCH="${BRANCH:-chamber-gpu}"
# Build in the directory the script is run from (e.g. cd /work/brunnels/jackplum/alamo
# then run this). Override with ALAMO_DIR=/abs/path.
ALAMO_DIR="${ALAMO_DIR:-$PWD}"
ACCOUNT="${ACCOUNT:-brunnels}"
BUILD_PARTITION="${BUILD_PARTITION:-nova}"   # CPU EPYC nodes; build needs no GPU
ARCHES="${ARCHES:-80 90}"                    # 80=A100, 90=H200
BUILD_JOBS="${BUILD_JOBS:-32}"
COMP="${COMP:-g++}"                          # nvcc host compiler (gcc is the safe choice)
EMAIL="${EMAIL:-jackplum@iastate.edu}"

echo -e "${BOLD}${BLUE}=== Alamo GPU build setup (NOVA) ===${NC}"
echo -e "  repo=${REPO_URL}  branch=${BRANCH}"
echo -e "  dir=${ALAMO_DIR}  archs='${ARCHES}'  account=${ACCOUNT}  build-partition=${BUILD_PARTITION}"

# ---- Modules (EDIT to match `module avail` on NOVA if names differ) ----------
echo -e "${YELLOW}Loading modules...${NC}"
module purge 2>/dev/null || true
module load cuda    2>/dev/null || module load cuda/12   2>/dev/null || echo -e "${RED}  WARN: load a cuda module manually${NC}"
module load gcc     2>/dev/null || module load gcc/12    2>/dev/null || echo -e "${RED}  WARN: load a gcc module manually${NC}"
module load openmpi 2>/dev/null || module load openmpi4  2>/dev/null || echo -e "${RED}  WARN: load an openmpi module manually${NC}"
module list 2>&1 | sed 's/^/    /' || true

# ---- 1. Get the source into ALAMO_DIR (the current directory) ---------------
mkdir -p "${ALAMO_DIR}"
ALAMO_DIR="$(cd "${ALAMO_DIR}" && pwd)"   # absolutize
echo -e "${YELLOW}Building in ${ALAMO_DIR}${NC}"
if [ -d "${ALAMO_DIR}/.git" ]; then
  echo -e "${YELLOW}Updating existing checkout in place...${NC}"
  git -C "${ALAMO_DIR}" fetch origin "${BRANCH}"
  git -C "${ALAMO_DIR}" checkout "${BRANCH}"
  git -C "${ALAMO_DIR}" pull --ff-only origin "${BRANCH}"
elif [ -f "${ALAMO_DIR}/configure" ] && [ -d "${ALAMO_DIR}/src" ]; then
  echo -e "${GREEN}  existing Alamo source found -- building in place${NC}"
else
  # Populate the current directory as a chamber-gpu checkout. Using init+fetch
  # (not `git clone`) so it works even though this dir already holds the script.
  echo -e "${YELLOW}Fetching ${BRANCH} into ${ALAMO_DIR}...${NC}"
  git -C "${ALAMO_DIR}" init -q
  git -C "${ALAMO_DIR}" remote add origin "${REPO_URL}" 2>/dev/null \
    || git -C "${ALAMO_DIR}" remote set-url origin "${REPO_URL}"
  git -C "${ALAMO_DIR}" fetch --depth 1 origin "${BRANCH}"
  git -C "${ALAMO_DIR}" checkout -b "${BRANCH}" FETCH_HEAD
fi
cd "${ALAMO_DIR}"

# ---- 2. Prime the AMReX checkout on the login node (needs internet) ---------
# configure clones AMReX into ext/ if missing; doing it here keeps the build
# job network-free (it only compiles).
echo -e "${YELLOW}Priming AMReX checkout (login node)...${NC}"
./configure --comp="${COMP}" --dim 2 --cuda "$(echo ${ARCHES} | awk '{print $1}')" --profile --get-eigen >/tmp/alamo_prime.log 2>&1 || {
  echo -e "${RED}configure prime failed -- see /tmp/alamo_prime.log${NC}"; tail -20 /tmp/alamo_prime.log; exit 1; }
echo -e "${GREEN}  AMReX present under ext/${NC}"

# ---- 3. Generate + submit the build job -------------------------------------
SB="${ALAMO_DIR}/build_alamo_gpu.sbatch"
cat > "${SB}" <<END_OF_SBATCH
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -J alamo_build
#SBATCH -D ${ALAMO_DIR}
#SBATCH --partition=${BUILD_PARTITION}
#SBATCH -N 1
#SBATCH -n ${BUILD_JOBS}
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=alamo_build.%j.out
#SBATCH --error=alamo_build.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${EMAIL}
set -euo pipefail
cd ${ALAMO_DIR}
module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true
for arch in ${ARCHES}; do
    echo "=== building cuda sm_\${arch} ==="
    ./configure --comp=${COMP} --dim 2 --cuda \${arch} --profile --get-eigen
    make -j${BUILD_JOBS} bin/alamo_gpu
done
echo "=== build complete ==="
ls -lh bin/alamo_gpu-2d*cuda* || true
END_OF_SBATCH

echo -e "${BOLD}${GREEN}Submitting build job...${NC}"
JOBID=$(sbatch --parsable "${SB}")
echo -e "${GREEN}  submitted job ${JOBID} on partition ${BUILD_PARTITION} (account ${ACCOUNT})${NC}"
echo ""
echo -e "${BOLD}Next:${NC}"
echo -e "  watch:   squeue -j ${JOBID}"
echo -e "  log:     tail -f ${ALAMO_DIR}/alamo_build.${JOBID}.out"
echo -e "  when done, the binaries are in ${ALAMO_DIR}/bin/ :"
echo -e "     alamo_gpu-2d-profile-cuda80-g++   (A100)"
echo -e "     alamo_gpu-2d-profile-cuda90-g++   (H200)"
echo -e "  then run:  sbatch ${ALAMO_DIR}/benchmark/nova_flame_gpu.slurm   (edit GPU_TYPE first)"
