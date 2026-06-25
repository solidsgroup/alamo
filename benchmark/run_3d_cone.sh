#!/usr/bin/env bash
# ============================================================================
# run_3d_cone.sh  --  NOVA launcher for the 3D conical-grain Flame test that
#                      dynamically picks its own hardware.
#
# Run this ON THE NOVA LOGIN NODE from the repo root:
#     bash benchmark/run_3d_cone.sh
#
# Replaces run_3d_cone_h200.sh / run_3d_cone_a100.sh / run_3d_cone_any_gpu.sh,
# which all requested a fixed, hand-picked hardware profile. This script
# instead asks the cluster what's actually free right now (via
# select_nova_resources.sh) and submits the closest-to-preferred profile that
# can start immediately:
#   - GPU:  h200 preferred, a100 fallback (a100-pcie is never selected -- its
#     CUDA arch on NOVA is unconfirmed)
#   - CPUs: 32 preferred, 16 minimum
#   - RAM:  128G preferred, 64G minimum
# If nothing matches any profile right now, it falls back to submitting the
# most-preferred profile (h200/32/128G) and lets it queue normally -- it will
# still run, just not necessarily immediately.
#
# It then:
#   1) writes the input deck (input_3d_cone),
#   2) builds ONLY the GPU binary needed for the selected GPU type,
#   3) submits the run with the selected GPU/CPU/RAM, afterok the build.
# ============================================================================
set -euo pipefail

GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; RED='\033[0;31m'; BOLD='\033[1m'; NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ---- knobs ----------------------------------------------------------------
ALAMO_DIR="${ALAMO_DIR:-${ROOT_DIR}}"
BRANCH="${BRANCH:-chamber-gpu}"
ACCOUNT="${ACCOUNT:-brunnels}"
PARTITION="${PARTITION:-nova}"
EMAIL="${EMAIL:-jackplum@iastate.edu}"
GPUS="${GPUS:-1}"
RUN_TIME="${RUN_TIME:-16:00:00}"
INPUT="${INPUT:-input_3d_cone}"
PLOT_FILE="${PLOT_FILE:-output_3d_cone}"
MODE="${MODE:-fast}"
SKIP_BUILD="${SKIP_BUILD:-0}"
BUILD_ONLY="${BUILD_ONLY:-0}"
DRY_RUN="${DRY_RUN:-0}"
EXPLAIN="${EXPLAIN:-0}"

# Hardware cutoffs -- consumed by select_nova_resources.sh.
export MIN_CPUS="${MIN_CPUS:-16}"
export PREF_CPUS="${PREF_CPUS:-32}"
export MIN_MEM_GB="${MIN_MEM_GB:-64}"
export PREF_MEM_GB="${PREF_MEM_GB:-128}"
export PARTITION

usage() {
    cat <<'EOF'
Usage: bash benchmark/run_3d_cone.sh [--dry-run] [--build-only] [--skip-build] [--explain]

One-shot NOVA launcher for the 3D conical-grain Flame test. Dynamically
selects GPU type (h200 preferred, a100 fallback), CPU count (32 preferred,
16 minimum), and RAM (128G preferred, 64G minimum) based on what's free on
the cluster right now, then builds and submits.

Flags:
  --explain     print the full per-node availability matrix before deciding
  --dry-run     print actions, submit nothing
  --build-only  submit the build and stop
  --skip-build  reuse an existing binary, no build job

Env knobs (defaults in brackets):
  ALAMO_DIR     repo root / working tree            [repo root from script path]
  BRANCH        branch the build pulls/builds        [chamber-gpu]
  ACCOUNT       slurm account                        [brunnels]
  PARTITION     slurm partition                      [nova]
  EMAIL         slurm mail-user                      [jackplum@iastate.edu]
  GPUS          number of GPUs per rank              [1]
  RUN_TIME      walltime for the run job             [16:00:00]
  INPUT         input deck name                      [input_3d_cone]
  PLOT_FILE     plotfile prefix                      [output_3d_cone]
  MODE          fast | bench                         [fast]
  MIN_CPUS / PREF_CPUS      hard floor / preferred CPU count   [16 / 32]
  MIN_MEM_GB / PREF_MEM_GB  hard floor / preferred RAM in GB   [64 / 128]
  FORCE_GPU_TYPE / FORCE_CPUS / FORCE_MEM_GB
                manual override -- skip cluster probing entirely and use
                these values (e.g. FORCE_GPU_TYPE=a100 FORCE_CPUS=16
                FORCE_MEM_GB=64)
  SKIP_BUILD    1 = reuse existing binary, no build   [0]
  BUILD_ONLY    1 = submit build and stop             [0]
  DRY_RUN       1 = print actions, submit nothing     [0]
EOF
}

for arg in "$@"; do
    case "$arg" in
        --dry-run)    DRY_RUN=1 ;;
        --build-only) BUILD_ONLY=1 ;;
        --skip-build) SKIP_BUILD=1 ;;
        --explain)    EXPLAIN=1 ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "error: unknown argument '$arg' (try --help)" >&2; exit 2 ;;
    esac
done

ALAMO_DIR="$(cd "${ALAMO_DIR}" && pwd)"

echo -e "${BOLD}${BLUE}============================================================${NC}"
echo -e "${BOLD}${BLUE} 3D conical-grain Flame test -- NOVA dynamic launcher${NC}"
echo -e "  repo        = ${ALAMO_DIR}"
echo -e "  branch      = ${BRANCH}"
echo -e "  account     = ${ACCOUNT}  partition=${PARTITION}"
echo -e "  CPU cutoffs = min=${MIN_CPUS} pref=${PREF_CPUS}"
echo -e "  RAM cutoffs = min=${MIN_MEM_GB}G pref=${PREF_MEM_GB}G"
echo -e "  SKIP_BUILD=${SKIP_BUILD}  BUILD_ONLY=${BUILD_ONLY}  DRY_RUN=${DRY_RUN}"
echo -e "${BOLD}${BLUE}============================================================${NC}"

# ---------------------------------------------------------------------------
# 0) Ask the cluster what's free right now and pick GPU/CPU/RAM.
# ---------------------------------------------------------------------------
echo
echo -e "${YELLOW}=== Selecting hardware profile ===${NC}"
# shellcheck source=./select_nova_resources.sh
source "${SCRIPT_DIR}/select_nova_resources.sh"

if [[ "${EXPLAIN}" -eq 1 ]]; then
    echo -e "${BLUE}--- node availability ---${NC}"
    _nova_explain || true
    echo
fi

select_nova_resources

case "${SELECTED_GPU_TYPE}" in
    h200) ARCH=90 ;;
    a100) ARCH=80 ;;
    *) echo -e "${RED}error: unexpected SELECTED_GPU_TYPE='${SELECTED_GPU_TYPE}'${NC}" >&2; exit 1 ;;
esac
CPUS_PER_TASK="${SELECTED_CPUS}"
MEM_GB="${SELECTED_MEM_GB}"

if [[ "${SELECTED_IMMEDIATE}" -eq 1 ]]; then
    AVAIL_NOTE="immediately available now on ${SELECTED_NODE}"
elif [[ "${SELECTED_IMMEDIATE}" -eq -1 ]]; then
    AVAIL_NOTE="manual FORCE_* override"
else
    AVAIL_NOTE="nothing free at any profile right now -- will queue normally"
fi

echo -e "${GREEN}  selected: GPU=${SELECTED_GPU_TYPE} (sm_${ARCH})  CPUs=${CPUS_PER_TASK}  RAM=${MEM_GB}G${NC}"
echo -e "${GREEN}  reason:   ${AVAIL_NOTE}${NC}"

# ---------------------------------------------------------------------------
# 1) Materialize the exact input deck into the repo root.
# ---------------------------------------------------------------------------
INPUT_PATH="${ALAMO_DIR}/${INPUT}"
echo
echo -e "${YELLOW}=== Writing input deck ${INPUT_PATH} ===${NC}"
if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "  (dry run) would write ${INPUT_PATH}"
else
    cat > "${INPUT_PATH}" <<'EOF_INPUT'
####################################################################################
# input_3d_cone — 3D GPU Flame test: conical propellant grain (NEW geometry)
#
# Geometry (NEW, never run before):
#   A SOLID CONE of propellant (eta=1) standing on the nozzle end, surrounded by
#   combustion gas (eta=0). The cone:
#     - base (large end) sits on the nozzle face, z = 0
#     - apex (point) terminates at z = 0.75 * H  (cone height = 75% of motor height)
#     - base radius = 0.25 * W = 50% of the motor-grain width (diameter)
#   The cone is generated analytically (Expression IC) — no BMP (BMP is 2D-only).
#
#   Motor grain = cylinder of radius R_motor (phi=1) inscribed in the square domain;
#   the square corners are casing (phi=0). The cone is well inside this cylinder.
#
# Coordinate convention: z is the motor axis. Nozzle/base at z=0 (bottom),
#   apex points up. x-y is the circular cross-section, centred at (cx,cy).
#
# Physics: PHASE-FIELD REGRESSION ONLY. Elastic is DISABLED (elastic.type=disable,
#   per D1: elastic is CPU-resident and the device elastic solve hits CUDA-719).
#   With type=disable the elastic model/bc blocks MUST be omitted or the strict
#   parser aborts.
#
# Grid: n_cell = 128 128 64, max_level = 2. Cubic cells:
#   dx = dy = 0.1754/128 = 1.370e-3 m ;  dz = 0.0877/64 = 1.370e-3 m.
#
# Purpose: LOCAL GPU correctness smoke (10 steps, plot every step) to verify the
#   cone IC is correct and regresses before scaling up on NOVA.
####################################################################################

alamo.program = flame
plot_file = output_3d_cone

amr.plot_dt = 0.01
amr.max_level = 2
amr.blocking_factor = 16
amr.max_grid_size = 128
amr.grid_eff = 0.7
amr.node.all = 1
amr.n_cell = 128 128 64
amr.thermo.int = 1
amr.thermo.plot_int = 1

# Domain: wide-shallow box. x-y is the circular cross-section (motor width),
#   z is the motor axis (motor height). Cubic cells.
geometry.prob_lo = 0.0_m 0.0_m 0.0_m
geometry.prob_hi = 0.1754_m 0.1754_m 0.0877_m
geometry.is_periodic = 0 0 0

system.length = m
system.time = s

timestep = 1.0e-4_s
stop_time = 6.5_s

# Phase field parameters
pf.eps = 5.0e-5_m
pf.lambda = 0.001_J/m^2
pf.kappa = 1.0_J/m^2
pf.relax_steps = 0
pf.w1 = 1.0_1
pf.w12 = 2.0_1
pf.w0 = 0.0_1
small = 1E-4

####################################################################################
# eta IC — SOLID CONE of propellant (eta=1 solid, eta=0 gas), smooth tanh interface.
#
#   cx, cy : motor axis (centre of the square cross-section)
#   rb     : cone base radius = 0.25 * W = 0.0425 m  (base diameter = 0.085 m = 50% W)
#   zap    : apex height = 0.75 * H = 0.0658 m       (cone height = 75% of motor)
#   w      : interface half-width (~ equilibrium width eps*sqrt(kappa/lambda)=1.6e-3)
#
#   Local cone radius at height z:  rc(z) = rb * (1 - z/zap).
#   eta = 0.5 + 0.5*tanh( (rc(z) - rho) / w ),  rho = sqrt((x-cx)^2 + (y-cy)^2).
#   For z > zap, rc(z) < 0 so eta -> 0 automatically (gas above the apex), i.e. the
#   tanh naturally closes the cone tip; no separate top cutoff is needed.
####################################################################################
pf.eta.ic.type = expression
pf.eta.ic.expression.constant.cx  = 0.0877
pf.eta.ic.expression.constant.cy  = 0.0877
pf.eta.ic.expression.constant.rb  = 0.0425
pf.eta.ic.expression.constant.zap = 0.0658
pf.eta.ic.expression.constant.w   = 0.002
pf.eta.ic.expression.region0 = "0.5 + 0.5*tanh((rb*(1.0 - z/zap) - sqrt((x-cx)^2 + (y-cy)^2))/w)"

# eta BCs — Neumann (zero-flux) on all six faces.
#   The cone is interior and surrounded by gas, so the walls must NOT pin eta to
#   solid. Zero-flux at z=0 also makes the cone base inhibited (does not burn
#   through the nozzle face), which is the physical case for a base-bonded grain.
pf.eta.bc.type = constant
pf.eta.bc.constant.type.xlo = neumann
pf.eta.bc.constant.type.xhi = neumann
pf.eta.bc.constant.type.ylo = neumann
pf.eta.bc.constant.type.yhi = neumann
pf.eta.bc.constant.type.zlo = neumann
pf.eta.bc.constant.type.zhi = neumann
pf.eta.bc.constant.val.xlo = 0.0
pf.eta.bc.constant.val.xhi = 0.0
pf.eta.bc.constant.val.ylo = 0.0
pf.eta.bc.constant.val.yhi = 0.0
pf.eta.bc.constant.val.zlo = 0.0
pf.eta.bc.constant.val.zhi = 0.0

####################################################################################
# phi IC — motor-grain mask (phi=1 inside the cylindrical motor, phi=0 casing).
#   Cylinder radius R = 0.085 m (motor width ~ domain width); corners are casing.
#   phi gates the combustion physics (get_L, get_qdot all scale by phi).
####################################################################################
phi.ic.type = expression
phi.ic.expression.constant.cx = 0.0877
phi.ic.expression.constant.cy = 0.0877
phi.ic.expression.constant.R  = 0.085
phi.ic.expression.constant.w  = 0.003
phi.ic.expression.region0 = "0.5 + 0.5*tanh((R - sqrt((x-cx)^2 + (y-cy)^2))/w)"

# Refinement criteria
amr.refinement_criterion = 0.001
amr.refinement_criterion_temp = 0.001_K

# Propellant model (homogenized; unchanged from the 2D/3D reference)
propellant.type = homogenize
propellant.homogenize.dispersion1 = 0.025_W/m/K
propellant.homogenize.dispersion2 = 1.2_kg/m^3
propellant.homogenize.dispersion3 = 1000_J/kg/K
propellant.homogenize.h1 = 1.5e6_W/m^2
propellant.homogenize.h2 = 8.0e6_W/m^2
propellant.homogenize.P_reference = 1_MPa
propellant.homogenize.rho_prop = 1745.0_kg/m^3
propellant.homogenize.k_prop = 1.5_W/m/K
propellant.homogenize.cp_prop = 1476.0_J/kg/K
propellant.homogenize.m_prop = 8000_cm/s
propellant.homogenize.mlocal_prop = 0.1_kg/m^2/s
propellant.homogenize.E_prop = 3500.0_K
propellant.homogenize.mob_prop = true
propellant.homogenize.pressure_exponent = 0.372531
propellant.homogenize.bound = 0_K
propellant.homogenize.bound_width = 50_K

# Thermal
thermal.on = 1
thermal.Tref = 300
thermal.Tfluid = 300
thermal.temp.bc.type = constant
thermal.temp.bc.constant.type.xlo = dirichlet
thermal.temp.bc.constant.type.xhi = dirichlet
thermal.temp.bc.constant.type.ylo = dirichlet
thermal.temp.bc.constant.type.yhi = dirichlet
thermal.temp.bc.constant.type.zlo = dirichlet
thermal.temp.bc.constant.type.zhi = dirichlet
thermal.temp.bc.constant.val.xlo = 300_K
thermal.temp.bc.constant.val.xhi = 300_K
thermal.temp.bc.constant.val.ylo = 300_K
thermal.temp.bc.constant.val.yhi = 300_K
thermal.temp.bc.constant.val.zlo = 300_K
thermal.temp.bc.constant.val.zhi = 300_K

temp.ic.type = constant
temp.ic.constant.value = 300_K

# Laser ignition: uniform volumetric heat for t < 0.025 s to light the grain.
laser.ic.type = expression
laser.ic.expression.region0 = "(t < 0.025) * (150000000)"
thermal.hc = 1.0

# Pressure / chamber
variable_pressure = 1
chamber.pressure = 0.101325_MPa
chamber.ballistic.At = 3.69565e-5_m^2
chamber.ballistic.T0 = 3200_K
chamber.ballistic.R = 287_J/kg/K
chamber.ballistic.gamma = 1.25
chamber.ballistic.pressure = 4.0e6_Pa

# Elastic DISABLED for the GPU path (D1). With type=disable the elastic model/bc
# blocks must be OMITTED or ALAMO's strict parser aborts on unused entries.
elastic.type = disable
EOF_INPUT
    echo -e "${GREEN}  wrote $(wc -l < "${INPUT_PATH}") lines${NC}"
fi

# ---------------------------------------------------------------------------
# 2) Submit the build job for ONLY the selected arch (unless SKIP_BUILD=1).
# ---------------------------------------------------------------------------
BUILD_JOBID=""
if [[ "${SKIP_BUILD}" -eq 1 ]]; then
    echo
    echo -e "${YELLOW}=== SKIP_BUILD=1: not submitting a build; expecting an existing sm_${ARCH} binary ===${NC}"
else
    echo
    echo -e "${YELLOW}=== Submitting GPU build job (ARCHES=${ARCH}) ===${NC}"
    if [[ "${DRY_RUN}" -eq 1 ]]; then
        echo "  (dry run) would run: ARCHES='${ARCH}' ACCOUNT='${ACCOUNT}' BRANCH='${BRANCH}' \\"
        echo "                       sh ${ROOT_DIR}/benchmark/build_alamo_nova_3d.sh"
    else
        BUILD_LOG="$(mktemp -t cone_nova_build.XXXXXX.log)"
        trap 'rm -f "${BUILD_LOG}"' EXIT
        ARCHES="${ARCH}" ACCOUNT="${ACCOUNT}" BRANCH="${BRANCH}" \
        BUILD_PARTITION="${PARTITION}" EMAIL="${EMAIL}" ALAMO_DIR="${ALAMO_DIR}" \
            sh "${ROOT_DIR}/benchmark/build_alamo_nova_3d.sh" 2>&1 | tee "${BUILD_LOG}"
        BUILD_JOBID="$(grep -oE 'submitted job [0-9]+' "${BUILD_LOG}" | tail -1 | awk '{print $3}')"
        if [[ -z "${BUILD_JOBID}" ]]; then
            echo -e "${RED}error: could not parse build job id from ${BUILD_LOG}${NC}" >&2
            exit 1
        fi
        echo -e "${GREEN}  build job id: ${BUILD_JOBID}${NC}"
    fi
fi

if [[ "${BUILD_ONLY}" -eq 1 ]]; then
    echo
    echo -e "${GREEN}BUILD_ONLY=1: stopping after build submission.${NC}"
    exit 0
fi

# ---------------------------------------------------------------------------
# 3) Write + submit the run job on the selected GPU/CPU/RAM (afterok build).
# ---------------------------------------------------------------------------
RUN_SB="${ALAMO_DIR}/run_3d_cone.sbatch"
echo
echo -e "${YELLOW}=== Writing run job ${RUN_SB} ===${NC}"
if [[ "${DRY_RUN}" -ne 1 ]]; then
cat > "${RUN_SB}" <<EOF_SB
#!/bin/bash
#SBATCH -A ${ACCOUNT}
#SBATCH -J cone_${SELECTED_GPU_TYPE}
#SBATCH -D ${ALAMO_DIR}
#SBATCH --partition=${PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:${SELECTED_GPU_TYPE}:${GPUS}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --mem=${MEM_GB}G
#SBATCH --time=${RUN_TIME}
#SBATCH --output=cone_${SELECTED_GPU_TYPE}.%j.out
#SBATCH --error=cone_${SELECTED_GPU_TYPE}.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${EMAIL}
# Single-GPU (1 MPI rank/GPU) 3D conical-grain Flame run.
# Hardware (GPU=${SELECTED_GPU_TYPE}, cpus=${CPUS_PER_TASK}, mem=${MEM_GB}G) was
# chosen dynamically by select_nova_resources.sh based on cluster state at
# submit time (${AVAIL_NOTE}).
set -euo pipefail
cd "\${SLURM_SUBMIT_DIR:-${ALAMO_DIR}}"

ARCH=${ARCH}
GPU_TYPE=${SELECTED_GPU_TYPE}
INPUT=${INPUT}
MODE=${MODE}
NRANKS="\${SLURM_NTASKS:-1}"

module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || module load gcc/12 2>/dev/null || true
module load openmpi 2>/dev/null || module load openmpi4 2>/dev/null || true

GPU_BIN="\$(ls -t bin/alamo_gpu-3d*cuda\${ARCH}*-g++ 2>/dev/null | head -1 || true)"
if [ -z "\${GPU_BIN}" ]; then
  echo "No bin/alamo_gpu-3d*cuda\${ARCH}*-g++ found -- build did not produce a \${GPU_TYPE} binary."
  exit 1
fi
echo "GPU binary: \${GPU_BIN}  (GPU_TYPE=\${GPU_TYPE} sm_\${ARCH}, \${NRANKS} rank/GPU, mode=\${MODE})"

if [ "\${MODE}" = "fast" ]; then
  OVERRIDES="amrex.async_out=0 amrex.async_out_nfiles=0 amrex.the_arena_is_managed=0"
else
  OVERRIDES="amrex.async_out=0 amrex.async_out_nfiles=0 amrex.the_arena_is_managed=1 tiny_profiler.device_synchronize_around_region=1"
fi

echo "=== launching cone run ==="
/usr/bin/time -v srun --mpi=pmix -n "\${NRANKS}" --gpus-per-task=1 \\
    "\${GPU_BIN}" "\${INPUT}" \\
    plot_file=${PLOT_FILE}_\${SLURM_JOB_ID} \\
    \${OVERRIDES}
echo "=== done -- plotfiles in ${PLOT_FILE}_\${SLURM_JOB_ID} (open celloutput.visit / nodeoutput.visit) ==="
EOF_SB
    echo -e "${GREEN}  wrote ${RUN_SB}${NC}"
fi

echo
echo -e "${YELLOW}=== Submitting run job ===${NC}"
DEP_ARG=()
if [[ -n "${BUILD_JOBID}" ]]; then
    DEP_ARG=(--dependency="afterok:${BUILD_JOBID}")
    echo "  (afterok:${BUILD_JOBID})"
fi
if [[ "${DRY_RUN}" -eq 1 ]]; then
    echo "  (dry run) would run: sbatch ${DEP_ARG[*]:-} ${RUN_SB}"
    echo
    echo -e "${GREEN}Dry run complete -- nothing submitted.${NC}"
    exit 0
fi
RUN_JOBID="$(sbatch --parsable "${DEP_ARG[@]}" "${RUN_SB}")"
RUN_JOBID="${RUN_JOBID%%;*}"
echo -e "${GREEN}  run job id: ${RUN_JOBID}${NC}"

echo
echo -e "${BOLD}${GREEN}=== Staged ===${NC}"
echo -e "  hardware  : GPU=${SELECTED_GPU_TYPE} (sm_${ARCH})  CPUs=${CPUS_PER_TASK}  RAM=${MEM_GB}G  (${AVAIL_NOTE})"
[[ -n "${BUILD_JOBID}" ]] && echo -e "  build job : ${BUILD_JOBID}"
echo -e "  run job   : ${RUN_JOBID}"
echo -e "  watch     : squeue -u \$USER"
echo -e "  run log   : tail -f ${ALAMO_DIR}/cone_${SELECTED_GPU_TYPE}.${RUN_JOBID}.out"
echo -e "  output    : ${ALAMO_DIR}/${PLOT_FILE}_${RUN_JOBID}/  (celloutput.visit, nodeoutput.visit)"
