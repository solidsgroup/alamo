#!/usr/bin/env bash
# ============================================================================
# run_3d_cone_a100.sh  --  partial launcher: kick the 3D cone-grain Flame case
#                          off on a single A100 (sm_80) on NOVA.
#
# Run on the NOVA login node from the repo root:
#     sh benchmark/run_3d_cone_a100.sh
#
# This is the run-only sibling of 3d_cone_grain_test.sh, retargeted to A100
# (more plentiful than H200). It does NOT re-run the full build orchestration:
#   - if an A100 (cuda80) binary already exists, it just submits the run;
#   - if not (your one-shot built ARCHES=90 / H200 only), it submits an sm_80
#     build first and chains the run afterok.
#
# Plot cadence is taken from the deck (input_3d_cone now uses amr.plot_dt=0.01,
# plot_int commented out) -- nothing is overridden here.
#
# NOTE: a sm_90 (H200) binary will NOT run on an A100; A100 needs its own sm_80
# build. NOTE: the run uses the input_3d_cone present in this working tree, so
# make sure it is the updated (plot_dt) deck on NOVA before launching.
# ============================================================================
set -euo pipefail

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

ALAMO_DIR="${ALAMO_DIR:-${ROOT_DIR}}"
BRANCH="${BRANCH:-chamber-gpu}"
ACCOUNT="${ACCOUNT:-brunnels}"
PARTITION="${PARTITION:-nova}"
EMAIL="${EMAIL:-jackplum@iastate.edu}"
GPU_TYPE="${GPU_TYPE:-a100}"          # this script targets A100 (sm_80)
ARCH=80
GPUS="${GPUS:-1}"
CPUS_PER_TASK="${CPUS_PER_TASK:-32}"
RUN_TIME="${RUN_TIME:-16:00:00}"
INPUT="${INPUT:-input_3d_cone}"
PLOT_FILE="${PLOT_FILE:-output_3d_cone}"
MODE="${MODE:-fast}"                  # fast (lean) | bench (TinyProfiler-friendly)
NO_BUILD="${NO_BUILD:-0}"             # 1 = never build; error if the A100 binary is missing
DRY_RUN="${DRY_RUN:-0}"

ALAMO_DIR="$(cd "${ALAMO_DIR}" && pwd)"
cd "${ALAMO_DIR}"

echo -e "${BOLD}${BLUE}=== 3D cone-grain Flame -- A100 (sm_${ARCH}) launcher ===${NC}"
echo -e "  repo  = ${ALAMO_DIR}"
echo -e "  input = ${INPUT}   mode = ${MODE}   gres = gpu:${GPU_TYPE}:${GPUS}"
echo -e "  NO_BUILD=${NO_BUILD}  DRY_RUN=${DRY_RUN}"

if [[ ! -f "${ALAMO_DIR}/${INPUT}" ]]; then
    echo -e "${RED}error: ${INPUT} not found in ${ALAMO_DIR}. Sync the updated deck first.${NC}" >&2
    exit 1
fi

# ---- 1) Ensure an A100 (sm_80) binary exists, building only if necessary ----
A100_BIN="$(ls -t bin/alamo_gpu-3d*cuda${ARCH}*-g++ 2>/dev/null | head -1 || true)"
BUILD_JOBID=""
if [[ -n "${A100_BIN}" ]]; then
    echo -e "${GREEN}  found A100 binary: ${A100_BIN} (no build needed)${NC}"
else
    echo -e "${YELLOW}  no bin/alamo_gpu-3d*cuda${ARCH}*-g++ found (your build was H200-only).${NC}"
    if [[ "${NO_BUILD}" -eq 1 ]]; then
        echo -e "${RED}error: NO_BUILD=1 and no A100 binary present. Build sm_80 first:${NC}" >&2
        echo "   ARCHES=80 sh benchmark/build_alamo_nova_3d.sh" >&2
        exit 1
    fi
    echo -e "${YELLOW}  submitting an sm_80 build (ARCHES=80) and chaining the run afterok...${NC}"
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        echo "   (dry run) ARCHES=80 ACCOUNT='${ACCOUNT}' BRANCH='${BRANCH}' sh benchmark/build_alamo_nova_3d.sh"
    else
        BUILD_LOG="$(mktemp -t cone_a100_build.XXXXXX.log)"
        trap 'rm -f "${BUILD_LOG}"' EXIT
        ARCHES=80 ACCOUNT="${ACCOUNT}" BRANCH="${BRANCH}" \
        BUILD_PARTITION="${PARTITION}" EMAIL="${EMAIL}" ALAMO_DIR="${ALAMO_DIR}" \
            sh "${ROOT_DIR}/benchmark/build_alamo_nova_3d.sh" 2>&1 | tee "${BUILD_LOG}"
        BUILD_JOBID="$(grep -oE 'submitted job [0-9]+' "${BUILD_LOG}" | tail -1 | awk '{print $3}')"
        [[ -z "${BUILD_JOBID}" ]] && { echo -e "${RED}error: could not parse sm_80 build job id${NC}" >&2; exit 1; }
        echo -e "${GREEN}  sm_80 build job id: ${BUILD_JOBID}${NC}"
    fi
fi

# ---- 2) Write the A100 run sbatch (sibling of nova_flame_gpu_3d.slurm) -------
RUN_SB="${ALAMO_DIR}/run_3d_cone_a100.sbatch"
if [[ "${DRY_RUN}" -ne 1 ]]; then
cat > "${RUN_SB}" <<EOF_SB
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -J cone_a100
#SBATCH -D ${ALAMO_DIR}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:${GPU_TYPE}:${GPUS}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --mem=128G
#SBATCH --time=${RUN_TIME}
#SBATCH --output=cone_a100.%j.out
#SBATCH --error=cone_a100.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${EMAIL}
# Single-A100 (1 MPI rank / GPU) 3D conical-grain Flame run.
set -euo pipefail
cd "\${SLURM_SUBMIT_DIR:-${ALAMO_DIR}}"

ARCH=${ARCH}
INPUT=${INPUT}
MODE=${MODE}
NRANKS="\${SLURM_NTASKS:-1}"

module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true

GPU_BIN="\$(ls -t bin/alamo_gpu-3d*cuda\${ARCH}*-g++ 2>/dev/null | head -1 || true)"
if [ -z "\${GPU_BIN}" ]; then
  echo "No bin/alamo_gpu-3d*cuda\${ARCH}*-g++ found -- sm_80 build did not produce an A100 binary."
  exit 1
fi
echo "GPU binary: \${GPU_BIN}  (A100 sm_\${ARCH}, \${NRANKS} rank/GPU, mode=\${MODE})"

if [ "\${MODE}" = "fast" ]; then
  OVERRIDES="amrex.async_out=0 amrex.async_out_nfiles=0 amrex.the_arena_is_managed=0"
else
  OVERRIDES="amrex.async_out=0 amrex.async_out_nfiles=0 amrex.the_arena_is_managed=1 tiny_profiler.device_synchronize_around_region=1"
fi

# Plot cadence comes from the deck (amr.plot_dt=0.01). Only the plotfile prefix
# is overridden, to tag it with the job id and avoid clobbering.
echo "=== launching A100 cone run ==="
/usr/bin/time -v srun --mpi=pmix -n "\${NRANKS}" --gpus-per-task=1 \\
    "\${GPU_BIN}" "\${INPUT}" plot_file=${PLOT_FILE}_\${SLURM_JOB_ID} \${OVERRIDES}
echo "=== done -- plotfiles in ${PLOT_FILE}_\${SLURM_JOB_ID} (celloutput.visit / nodeoutput.visit) ==="
EOF_SB
    echo -e "${GREEN}  wrote ${RUN_SB}${NC}"
fi

# ---- 3) Submit the run (afterok the build, if we built one) -----------------
DEP_ARG=()
[[ -n "${BUILD_JOBID}" ]] && DEP_ARG=(--dependency="afterok:${BUILD_JOBID}")

if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo -e "${GREEN}(dry run) would: sbatch ${DEP_ARG[*]:-} ${RUN_SB}${NC}"
    exit 0
fi

RUN_JOBID="$(sbatch --parsable "${DEP_ARG[@]}" "${RUN_SB}")"
RUN_JOBID="${RUN_JOBID%%;*}"

echo
echo -e "${BOLD}${GREEN}=== Staged ===${NC}"
[[ -n "${BUILD_JOBID}" ]] && echo -e "  sm_80 build job : ${BUILD_JOBID}"
echo -e "  A100 run job    : ${RUN_JOBID}"
echo -e "  watch           : squeue -u \$USER"
echo -e "  run log         : tail -f ${ALAMO_DIR}/cone_a100.${RUN_JOBID}.out"
echo -e "  output          : ${ALAMO_DIR}/${PLOT_FILE}_${RUN_JOBID}/"
echo
echo -e "  (tip: cancel the stuck H200 job with  scancel <jobid>  if you no longer want it.)"
