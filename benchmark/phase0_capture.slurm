#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J phase0_capture
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=phase0_capture.%j.out
#SBATCH --error=phase0_capture.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jackplum@iastate.edu
# ===========================================================================
# phase0_capture.slurm -- one submission captures the whole Phase 0 figure set.
#
# Campaign: docs/agent_plans/20260727-gpu-memory-strategy/PLAN.md
# Task:     docs/agent_plans/20260727-phase0-baseline/PLAN.md  (steps N1-N4)
#
# WHY ONE JOB
# -----------
# GPU access here is queued, not interactive (campaign PLAN §2). One job per
# figure would spend the campaign in the queue. This captures F1-F10 for one
# dimension in a single allocation; submit it twice (DIM=2, DIM=3).
#
# LEGS (select with LEGS="...", default all, in this order)
#   env      Node inventory: GPU model, driver, toolkit, ncu counter permission,
#            and the GPU-aware MPI probe carried over from Phase 0.5 (N4).
#   timing   F10. NON-profile binary, managed arena vs device arena, MAX_STEP
#            steps each, external wall + wall/step. The only admissible timing
#            numbers in this job.
#   flip     N2. Device arena + abort_on_out_of_gpu_memory=1. INVENTORY ONLY --
#            a crash here is the deliverable, not a failure. Never fix in-job.
#   nsys     F1,F2,F3,F8,F9. Profile binary (BL_PROFILE regions emit NVTX under
#            CUDA), CUDA+NVTX+OSRT+MPI trace, CUDA memory usage, and unified-
#            memory page-fault counters. Exports stats CSVs so collection is
#            small files, not just a .nsys-rep.
#   ncu      F4,F5,F6,F7. Top kernels only, a few launches each, selected by
#            NVTX range rather than kernel name -- AMReX kernels all present as
#            `launch_global`, which is why a name-based filter returns nothing.
#            Skipped automatically if counters are barred (ERR_NVGPUCTRPERM).
#   smooth   MLMG pre/post-smoothing A/B (4/4 vs 2/2) at production cadence,
#            reporting MLMG ITERATION COUNTS alongside wall. Not in the default
#            leg set -- run it on the deck whose smoothing question is open.
#
# USAGE
#   DECK=input_copy sbatch --gres=gpu:a100:1 benchmark/phase0_capture.slurm
#   DECK=input_3d_centre_bore_128_a2 sbatch --gres=gpu:a100:1 benchmark/phase0_capture.slurm
#   DRYRUN=1 bash benchmark/phase0_capture.slurm      # print every command, run none
#   LEGS="env nsys" DECK=input sbatch --gres=gpu:a100:1 benchmark/phase0_capture.slurm
#   LEGS=smooth SMOOTH_STEP=800 DECK=input_copy sbatch --gres=gpu:a100:1 ...
#
# Build first (BOTH variants -- timing needs the non-profile binary):
#   VARIANTS="profile plain" DIMS="2 3" sh benchmark/build_alamo_nova.sh
# ===========================================================================
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"

# --- Selectors --------------------------------------------------------------
# DECK is the primary selector; DIM is derived from it. One job per deck.
DECK="${DECK:-input_copy}"
GPU_TYPE="${GPU_TYPE:-a100}"
LEGS="${LEGS:-env timing flip nsys ncu}"
DRYRUN="${DRYRUN:-0}"
NCU_LAUNCHES="${NCU_LAUNCHES:-3}"   # ncu serializes; a handful of launches is plenty
MGS="${MGS:-64}"                    # NOVA-like multi-box layout
ARENA_FRAC="${ARENA_FRAC:-0.5}"     # fraction of device memory claimed up front, per rank

# --- Per-deck horizons ------------------------------------------------------
#
# CRITICAL: elastic is ~95% of GPU wall, and every production deck fires it on an
# interval (40-50 steps) after a tstart delay. A 3- or 10-step capture therefore
# contains ZERO elastic work and measures nothing that matters. Horizons below are
# set so each capture spans at least TWO elastic solves at PRODUCTION cadence --
# no elastic.interval override, because overriding it measures a configuration
# nobody runs.
#
#   deck                        dt      tstart   interval  first solve  2nd solve
#   input_copy                  2.5e-4  0.01 s   40        step ~40     step ~80
#   input                       1.0e-4  0.5e-4   50        step ~50     step ~100
#   input_3d_centre_bore_128_a2 1.0e-4  0.5e-4   50        step ~50     step ~100
case "${DECK}" in
  input_copy)
    DIM=2; DEF_MAX_STEP=90;  DEF_TRACE_STEP=90 ;;
  input)
    DIM=2; DEF_MAX_STEP=110; DEF_TRACE_STEP=110 ;;
  input_3d_centre_bore_128_a2)
    DIM=3; DEF_MAX_STEP=110; DEF_TRACE_STEP=110 ;;
  *)
    echo "Unknown DECK='${DECK}'. Add its dt/tstart/interval to the table above"
    echo "and set a horizon that spans two elastic solves -- do not guess."
    exit 2 ;;
esac
MAX_STEP="${MAX_STEP:-${DEF_MAX_STEP}}"
TRACE_STEP="${TRACE_STEP:-${DEF_TRACE_STEP}}"

case "${GPU_TYPE}" in
  a100) ARCH=80 ;;
  h200) ARCH=90 ;;
  *) echo "Unknown GPU_TYPE=${GPU_TYPE} (a100|h200)"; exit 2 ;;
esac

INPUT="${DECK}"

JOB="${SLURM_JOB_ID:-dryrun}"
OUT="${OUT:-$PWD/benchmark/_phase0_${DECK}_${GPU_TYPE}_${JOB}}"

# --- Modules: PINNED. Campaign §2 -- verify spellings once per cluster change.
# Verified on nova-login-2, 2026-07-27: cuda/12.8.1 ships nsys 2024.6.2.225 and
# ncu 2025.1.1.0. Do not switch to bare `module load cuda` -- the version that
# resolves to decides every profiler flag spelling below.
CUDA_MOD="${CUDA_MOD:-cuda/12.8.1-hm53cro}"
GCC_MOD="${GCC_MOD:-gcc/11.4.1-l3zl3n5}"
MPI_MOD="${MPI_MOD:-openmpi/4.1.8-lctrfpx}"

module purge 2>/dev/null || true
module load "${CUDA_MOD}" 2>/dev/null || module load cuda 2>/dev/null || true
module load "${GCC_MOD}"  2>/dev/null || module load gcc  2>/dev/null || true
module load "${MPI_MOD}"  2>/dev/null || module load openmpi 2>/dev/null || true

# --- Binaries ---------------------------------------------------------------
PROF_BIN="${PROF_BIN:-bin/alamo_gpu-${DIM}d-profile-cuda${ARCH}-g++}"
PLAIN_BIN="${PLAIN_BIN:-bin/alamo_gpu-${DIM}d-cuda${ARCH}-g++}"

# --- Deck overrides ---------------------------------------------------------
# amr.plot_int=-1: plotfile IO is a legitimate transfer (campaign §3) but it is
# not what F2 is asking about, and it would dominate a 3-step trace.
COMMON=( "${INPUT}"
         "amr.plot_int=-1"
         "amr.thermo.plot_int=1"
         "amr.thermo.int=1"
         "amr.max_grid_size=${MGS}"
         "elastic.solver.verbose=0"
         "elastic.print_model=0"
         "amrex.abort_on_out_of_gpu_memory=1" )
MANAGED=( "amrex.the_arena_is_managed=1" )
DEVICE=(  "amrex.the_arena_is_managed=0" )

mkdir -p "${OUT}"
say () { printf '\n=== %s ===\n' "$*"; }
run () {  # echo then run (or just echo under DRYRUN)
  printf '+ %s\n' "$*"
  [ "${DRYRUN}" = 1 ] && return 0
  "$@"
}
have_leg () { case " ${LEGS} " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

say "phase0_capture  DECK=${DECK} (${DIM}D) GPU=${GPU_TYPE}(sm_${ARCH}) job=${JOB}"
echo "  legs='${LEGS}'  out=${OUT}"
echo "  horizons: timing=${MAX_STEP} trace=${TRACE_STEP} steps (>=2 elastic solves, production cadence)"
echo "  profile binary=${PROF_BIN}"
echo "  plain   binary=${PLAIN_BIN}"

if [ "${DRYRUN}" != 1 ]; then
  for b in "${PROF_BIN}" "${PLAIN_BIN}"; do
    [ -x "${b}" ] || { echo "MISSING ${b} -- run: VARIANTS=\"profile plain\" DIMS=\"${DIM}\" sh benchmark/build_alamo_nova.sh"; exit 2; }
  done
fi

# ---------------------------------------------------------------- leg: env ---
if have_leg env; then
  say "LEG env -- node inventory, ncu counter permission, GPU-aware MPI (N4)"
  L="${OUT}/env"; mkdir -p "${L}"
  {
    echo "host=$(hostname)  date=$(date -Is)"
    echo "--- source provenance (written by phase0_capture.sh push; the remote"
    echo "    checkout's own git metadata does NOT describe the pushed tree)"
    cat benchmark/_pushed_rev.txt 2>/dev/null || echo "  MISSING benchmark/_pushed_rev.txt -- ledger row will have no trustworthy sha"
    echo "--- nvidia-smi"; nvidia-smi 2>&1 | head -20
    echo "--- compute cap / mem"
    nvidia-smi --query-gpu=name,compute_cap,memory.total,driver_version --format=csv 2>&1
    echo "--- toolkit"
    which nvcc nsys ncu 2>&1; nvcc --version 2>&1 | tail -2
    nsys --version 2>&1 | head -1; ncu --version 2>&1 | grep -i version | head -1
    echo "--- MPI"
    which mpiexec srun 2>&1
    echo "MPICH_GPU_SUPPORT_ENABLED=${MPICH_GPU_SUPPORT_ENABLED:-<unset>}"
    ompi_info --version 2>&1 | head -1
    echo "--- GPU-aware MPI verdict (Phase 0.5 §2 carry-over):"
    ompi_info --parsable --all 2>/dev/null | grep -E 'mpi_built_with_cuda_support' || echo "  (no ompi_info cuda line; not Open MPI?)"
    ompi_info --parsable --all 2>/dev/null | grep -iE 'accelerator|smcuda' | head -5
    echo "--- arena sizing input"
    echo "device_total_MiB=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | head -1)"
  } > "${L}/inventory.txt" 2>&1
  [ "${DRYRUN}" = 1 ] || cat "${L}/inventory.txt"

  # ncu counter permission: a real capture on a trivial run is the only honest test.
  say "ncu counter permission probe"
  NCU_OK=0
  if [ "${DRYRUN}" = 1 ]; then
    echo "+ ncu --target-processes all --launch-count 1 --section SpeedOfLight -o ${L}/ncu_probe ${PLAIN_BIN} ${INPUT} max_step=1 ..."
  else
    ncu --target-processes all --launch-count 1 --section SpeedOfLight \
        -o "${L}/ncu_probe" --force-overwrite \
        "./${PLAIN_BIN}" "${COMMON[@]}" "${MANAGED[@]}" max_step=1 stop_time=1e99_s \
        plot_file="${L}/probe_plot" > "${L}/ncu_probe.log" 2>&1
    if grep -qi "ERR_NVGPUCTRPERM\|The user does not have permission" "${L}/ncu_probe.log"; then
      echo "  ncu counters BARRED on this node (ERR_NVGPUCTRPERM) -- ncu leg will be skipped."
      echo "  Substitute: benchmark/res_usage.sh (root-free register/occupancy proxy)."
    else
      NCU_OK=1; echo "  ncu counters permitted."
    fi
    echo "${NCU_OK}" > "${L}/ncu_ok"
  fi
fi

# ------------------------------------------------------------- leg: timing ---
# F10. Non-profile binary only. See the header note on why.
if have_leg timing; then
  say "LEG timing -- F10 wall/step, managed vs device arena (NON-profile binary)"
  for mode in managed device; do
    L="${OUT}/timing/${mode}"; mkdir -p "${L}"
    if [ "${mode}" = managed ]; then ARENA=( "${MANAGED[@]}" ); else ARENA=( "${DEVICE[@]}" ); fi
    say "timing/${mode}"
    if [ "${DRYRUN}" = 1 ]; then
      echo "+ srun --mpi=pmix -n1 --gpus-per-task=1 ./${PLAIN_BIN} ${COMMON[*]} ${ARENA[*]} max_step=${MAX_STEP} stop_time=1e99_s plot_file=${L}/plot"
    else
      T0=$(date +%s.%N)
      srun --mpi=pmix -n1 --gpus-per-task=1 \
        "./${PLAIN_BIN}" "${COMMON[@]}" "${ARENA[@]}" max_step="${MAX_STEP}" \
        stop_time=1e99_s plot_file="${L}/plot" > "${L}/run.log" 2>&1
      RC=$?; T1=$(date +%s.%N)
      WALL=$(awk -v a="${T0}" -v b="${T1}" 'BEGIN{printf "%.3f", b-a}')
      PER=$(awk -v w="${WALL}" -v n="${MAX_STEP}" 'BEGIN{printf "%.4f", w/n}')
      echo "  rc=${RC} wall=${WALL}s steps=${MAX_STEP} wall_per_step=${PER}s"
      printf 'mode=%s rc=%s wall_s=%s steps=%s wall_per_step_s=%s\n' \
        "${mode}" "${RC}" "${WALL}" "${MAX_STEP}" "${PER}" > "${L}/timing.txt"
    fi
  done
  echo "NOTE: wall includes startup. Campaign requires a startup calibration alongside"
  echo "      wall/step; derive it from the two horizons or a max_step=1 run."
fi

# --------------------------------------------------------------- leg: flip ---
# N2. The deliverable is the crash set, classified. Fix NOTHING in this job.
if have_leg flip; then
  say "LEG flip -- N2 device-arena inventory (a crash here is the RESULT, not a failure)"
  L="${OUT}/flip"; mkdir -p "${L}"
  if [ "${DRYRUN}" = 1 ]; then
    echo "+ srun --mpi=pmix -n1 --gpus-per-task=1 ./${PLAIN_BIN} ${COMMON[*]} ${DEVICE[*]} amrex.the_arena_init_size=<frac ${ARENA_FRAC}> max_step=${MAX_STEP} ..."
  else
    TOT=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits 2>/dev/null | head -1)
    INIT=$(awk -v t="${TOT:-40960}" -v f="${ARENA_FRAC}" 'BEGIN{printf "%d", t*1048576*f}')
    echo "  device_total=${TOT} MiB  arena_init_size=${INIT} B (frac ${ARENA_FRAC}, per rank)"
    srun --mpi=pmix -n1 --gpus-per-task=1 \
      "./${PLAIN_BIN}" "${COMMON[@]}" "${DEVICE[@]}" \
      amrex.the_arena_init_size="${INIT}" max_step="${MAX_STEP}" \
      stop_time=1e99_s plot_file="${L}/plot" > "${L}/run.log" 2>&1
    echo "  rc=$? -> ${L}/run.log"
    grep -inE "illegal|CUDA error|Abort|out of memory|segmentation" "${L}/run.log" | head -20 \
      | tee "${L}/failures.txt"
    echo "  classify each hit per campaign PLAN §7: timestep loop | diagnostic | checkpoint | debug leftover"
  fi
fi

# --------------------------------------------------------------- leg: nsys ---
if have_leg nsys; then
  say "LEG nsys -- F1,F2,F3,F8,F9 (profile binary: BL_PROFILE regions emit NVTX)"
  L="${OUT}/nsys"; mkdir -p "${L}"
  NSYS_ARGS=( profile
    -t cuda,nvtx,osrt,mpi
    --cuda-memory-usage=true
    --cuda-um-cpu-page-faults=true
    --cuda-um-gpu-page-faults=true
    --force-overwrite true
    -o "${L}/trace" )
  if [ "${DRYRUN}" = 1 ]; then
    echo "+ srun --mpi=pmix -n1 --gpus-per-task=1 nsys ${NSYS_ARGS[*]} ./${PROF_BIN} ${COMMON[*]} ${MANAGED[*]} max_step=${TRACE_STEP} ..."
    echo "+ nsys stats --report nvtx_sum,cuda_gpu_kern_sum,cuda_gpu_mem_size_sum,cuda_gpu_mem_time_sum,cuda_api_sum --format csv ${L}/trace.nsys-rep"
  else
    srun --mpi=pmix -n1 --gpus-per-task=1 \
      nsys "${NSYS_ARGS[@]}" \
      "./${PROF_BIN}" "${COMMON[@]}" "${MANAGED[@]}" max_step="${TRACE_STEP}" \
      stop_time=1e99_s plot_file="${L}/plot" > "${L}/run.log" 2>&1
    echo "  nsys rc=$?"
    # Export CSVs so the collect step moves kilobytes, not gigabytes. Report names
    # vary across nsys versions -- try each, record which ones existed.
    for rep in nvtx_sum cuda_gpu_kern_sum cuda_gpu_mem_size_sum cuda_gpu_mem_time_sum \
               cuda_api_sum cuda_gpu_trace um_total_sum um_cpu_page_faults_sum; do
      if nsys stats --report "${rep}" --format csv --output "${L}/${rep}" \
           "${L}/trace.nsys-rep" >> "${L}/stats.log" 2>&1; then
        echo "  stats ok: ${rep}"
      else
        echo "  stats MISSING in this nsys version: ${rep}" | tee -a "${L}/stats_missing.txt"
      fi
    done
  fi
  echo "  F1 is read first at every phase boundary (campaign §12). Open trace.nsys-rep,"
  echo "  collapse to one timestep, and look at the gaps before reading any aggregate."
fi

# ---------------------------------------------------------------- leg: ncu ---
if have_leg ncu; then
  say "LEG ncu -- F4,F5,F6,F7 top kernels only"
  L="${OUT}/ncu"; mkdir -p "${L}"
  NCU_OK_FILE="${OUT}/env/ncu_ok"
  if [ "${DRYRUN}" != 1 ] && [ -f "${NCU_OK_FILE}" ] && [ "$(cat "${NCU_OK_FILE}")" != "1" ]; then
    echo "  SKIPPED: counters barred on this node (see env leg). Use benchmark/res_usage.sh."
  else
    # Select by NVTX range, not kernel name: AMReX launches all appear as
    # `launch_global`, which is why a --kernel-name filter comes back empty
    # (see the ncu-empty root-cause finding).
    SECTIONS=( --section SpeedOfLight
               --section SpeedOfLight_RooflineChart
               --section MemoryWorkloadAnalysis
               --section Occupancy
               --section LaunchStats )
    for range in "Operator::Elastic::Fapply()" "Operator::Fsmooth()" "MLMG::mgVcycle()"; do
      tag=$(echo "${range}" | tr -cd '[:alnum:]')
      if [ "${DRYRUN}" = 1 ]; then
        echo "+ srun --mpi=pmix -n1 --gpus-per-task=1 ncu --nvtx --nvtx-include \"${range}/\" --launch-count ${NCU_LAUNCHES} ${SECTIONS[*]} -o ${L}/${tag} ./${PROF_BIN} ${COMMON[*]} ${MANAGED[*]} max_step=${TRACE_STEP} ..."
      else
        srun --mpi=pmix -n1 --gpus-per-task=1 \
          ncu --nvtx --nvtx-include "${range}/" --launch-count "${NCU_LAUNCHES}" \
              "${SECTIONS[@]}" -o "${L}/${tag}" --force-overwrite \
          "./${PROF_BIN}" "${COMMON[@]}" "${MANAGED[@]}" max_step="${TRACE_STEP}" \
          stop_time=1e99_s plot_file="${L}/plot_${tag}" > "${L}/${tag}.log" 2>&1
        echo "  ${range} rc=$? -> ${L}/${tag}.ncu-rep"
        ncu --import "${L}/${tag}.ncu-rep" --csv --page details > "${L}/${tag}.csv" 2>>"${L}/${tag}.log" \
          && echo "    csv exported"
      fi
    done
  fi
fi

# ------------------------------------------------------------- leg: smooth ---
# The open question both prior tasks refused to close.
#
#   docs/agent_plans/20260721-fapply-runtime-optimization/results/RESULT.md
#     retains 2/2: -12.79% wall 2D, -21.05% 3D, A100 -16.56% wall / -22.54% MLMG.
#     Says twice: "not a longer-evolution stability claim". 3D final residual
#     9.72181e-09 against a 1e-08 gate -- clearing convergence by 3% on a
#     TWO-STEP deck.
#   docs/agent_plans/20260727-gpu-optimization-investigation/results/RESULT.md
#     R1: 2/2 "costs 5-41% more MLMG iterations", declines to apply it because
#     "the longest validation run here was 125 steps, and production runs to
#     stop_time = 6.5_s ... that margin is exactly what gets consumed."
#
# Both asked for a long-horizon run on one geometry. This is that run. It is
# NOT a timing leg: the deliverable is ITERATION COUNT and whether the solve
# survives, with wall as secondary.
#
# Note elastic.solver.verbose is deliberately NOT zeroed here -- MLMG's
# iteration/residual lines are the measurement. Everywhere else verbose=0.
if have_leg smooth; then
  SMOOTH_STEP="${SMOOTH_STEP:-800}"   # >> the 125-step longest prior validation
  say "LEG smooth -- 4/4 vs 2/2 at production cadence, ${SMOOTH_STEP} steps"
  echo "  deliverable = MLMG iteration counts + survival, NOT wall"
  for arm in 4 2; do
    L="${OUT}/smooth/${arm}x${arm}"; mkdir -p "${L}"
    ARM=( "elastic.solver.pre_smooth=${arm}" "elastic.solver.post_smooth=${arm}"
          "elastic.solver.verbose=4" )
    if [ "${DRYRUN}" = 1 ]; then
      echo "+ srun --mpi=pmix -n1 --gpus-per-task=1 ./${PLAIN_BIN} ${INPUT} ... ${ARM[*]} max_step=${SMOOTH_STEP}"
      continue
    fi
    T0=$(date +%s.%N)
    srun --mpi=pmix -n1 --gpus-per-task=1 \
      "./${PLAIN_BIN}" "${INPUT}" \
      amr.plot_int=-1 amr.thermo.plot_int=1 amr.thermo.int=1 \
      "amr.max_grid_size=${MGS}" elastic.print_model=0 \
      amrex.abort_on_out_of_gpu_memory=1 "${MANAGED[@]}" "${ARM[@]}" \
      max_step="${SMOOTH_STEP}" stop_time=1e99_s \
      plot_file="${L}/plot" > "${L}/run.log" 2>&1
    RC=$?; T1=$(date +%s.%N)
    WALL=$(awk -v a="${T0}" -v b="${T1}" 'BEGIN{printf "%.3f", b-a}')
    # MLMG prints one "Final Iter." line per solve at verbose>=1; count them and
    # summarise, rather than trusting a single end-of-run number.
    ITERS=$(grep -oE "Final Iter\.? +[0-9]+" "${L}/run.log" 2>/dev/null | grep -oE "[0-9]+$" | tr '\n' ' ')
    NSOLVE=$(grep -c "Final Iter" "${L}/run.log" 2>/dev/null)
    FAILED=$(grep -ciE "MLMG failed|Abort|did not converge|illegal" "${L}/run.log" 2>/dev/null)
    {
      echo "arm=${arm}/${arm}"
      echo "rc=${RC} wall_s=${WALL} steps=${SMOOTH_STEP}"
      echo "elastic_solves=${NSOLVE}"
      echo "failure_lines=${FAILED}"
      echo "mlmg_iters_per_solve=${ITERS}"
      echo "${ITERS}" | tr ' ' '\n' | awk 'NF{n++;s+=$1;if($1>m)m=$1} END{if(n)printf "mlmg_iters_mean=%.2f\nmlmg_iters_max=%d\n",s/n,m}'
    } > "${L}/summary.txt"
    cat "${L}/summary.txt"
    if [ "${FAILED}" -gt 0 ] 2>/dev/null; then
      echo "  --- first failure line:"
      grep -inE "MLMG failed|Abort|did not converge|illegal" "${L}/run.log" | head -3
    fi
  done
  echo "VERDICT INPUT: compare mlmg_iters_mean/max and survival between arms."
  echo "  2/2 winning wall while iters stay bounded  -> the prior claim extends."
  echo "  2/2 iters climbing toward max_iter, or failing -> it does not, and"
  echo "  the two-step measurements were sampling an easy conditioning state."
fi

say "phase0_capture done -- artifacts under ${OUT}"
echo "Collect with: bash benchmark/phase0_capture.sh collect ${OUT}"
