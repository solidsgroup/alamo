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
#
# VARIANTS (default "profile") selects which builds to produce:
#   profile -> bin/alamo_gpu-<dim>d-profile-cuda<arch>-g++   diagnostic binary.
#              TinyProfiler is on, and on CUDA every BL_PROFILE region also emits
#              an NVTX range (vendored AMReX 26.06 AMReX_TinyProfiler.cpp:134,211),
#              which is what makes an nsys timeline readable.
#   plain   -> bin/alamo_gpu-<dim>d-cuda<arch>-g++            timing binary.
#              No profiler instrumentation. Wall-clock claims must come from this
#              one -- a profile build pays a push/pop per region (16k+ Fapply
#              ranges in a 2-step 2D run) and MODE=bench additionally sets
#              tiny_profiler.device_synchronize_around_region=1, which serializes
#              the asynchrony a memory-strategy campaign is trying to measure.
#   Build both:  VARIANTS="profile plain" sh build_alamo_nova.sh
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
VARIANTS="${VARIANTS:-profile}"              # "profile", "plain", or "profile plain"
DIMS="${DIMS:-2}"                            # 2, 3, or "2 3"
BUILD_JOBS="${BUILD_JOBS:-64}"
BUILD_MEM="${BUILD_MEM:-64G}"
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
if [ "${SKIP_GIT:-0}" = 1 ]; then
  # The working tree was placed here by an rsync push (benchmark/phase0_capture.sh
  # push) and IS the thing to build. Pulling would fight it: a dirty local branch
  # makes `pull --ff-only` abort, and if it did succeed it would silently discard
  # the pushed source. Provenance for what is actually here lives in
  # benchmark/_pushed_rev.txt, not in this checkout's git metadata.
  echo -e "${YELLOW}SKIP_GIT=1 -- building the tree as-is, no fetch/checkout/pull${NC}"
  if [ -f "${ALAMO_DIR}/benchmark/_pushed_rev.txt" ]; then
    sed 's/^/    /' "${ALAMO_DIR}/benchmark/_pushed_rev.txt"
  else
    echo -e "${RED}  WARNING: no benchmark/_pushed_rev.txt -- no trustworthy source provenance${NC}"
  fi
elif [ -d "${ALAMO_DIR}/.git" ]; then
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
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${BUILD_JOBS}
#SBATCH --mem=${BUILD_MEM}
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
for dim in ${DIMS}; do
for arch in ${ARCHES}; do
for variant in ${VARIANTS}; do
    case "\${variant}" in
      profile) VFLAG="--profile" ;;
      plain)   VFLAG="" ;;
      *) echo "unknown VARIANT '\${variant}' (use profile|plain)"; exit 1 ;;
    esac
    echo "=== building dim=\${dim} cuda sm_\${arch} variant=\${variant} ==="
    ./configure --comp=${COMP} --dim \${dim} --cuda \${arch} \${VFLAG} --get-eigen
    make -j\${SLURM_CPUS_PER_TASK:-${BUILD_JOBS}} bin/alamo_gpu
done
done
done
echo "=== build complete ==="
ls -lh bin/alamo_gpu-*d*cuda* || true
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
