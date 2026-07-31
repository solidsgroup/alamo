#!/bin/bash
#SBATCH -A brunnels
#SBATCH -J phase0_capture
#SBATCH --partition=nova
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
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
#   arena    F8/T1/T4/T5. AMReX request-layer allocation counts, live bytes,
#            and high-water at a startup and full-horizon endpoint.
#   nsys     F1,F2,F3,F9. Profile binary (BL_PROFILE regions emit NVTX under
#            CUDA), CUDA+NVTX trace (NOT mpi/OSRT -- see the leg), CUDA memory
#            usage, and unified-memory page-fault counters. Exports stats CSVs
#            so collection is small files, not just a .nsys-rep.
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
LEGS="${LEGS:-env timing flip arena nsys ncu}"
DRYRUN="${DRYRUN:-0}"
NCU_LAUNCHES="${NCU_LAUNCHES:-3}"   # ncu serializes; a handful of launches is plenty
NCU_TARGETS="${NCU_TARGETS:-4}"      # distinct functions from nsys's top-10 table
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
ELASTIC_INTERVAL="$(awk -F= '
  {
    key=$1
    gsub(/[[:space:]]/, "", key)
    if (key == "elastic.interval") {
      value=$2
      sub(/#.*/, "", value)
      gsub(/[[:space:]]/, "", value)
      print value
      exit
    }
  }' "${DECK}")"
if ! [[ "${ELASTIC_INTERVAL}" =~ ^[1-9][0-9]*$ ]]; then
  echo "Could not read a positive elastic.interval from ${DECK}; refusing a capture with an unknown elastic horizon." >&2
  exit 2
fi
MIN_ELASTIC_STEPS=$((2 * ELASTIC_INTERVAL))
if [ "${MAX_STEP}" -lt "${MIN_ELASTIC_STEPS}" ] || [ "${TRACE_STEP}" -lt "${MIN_ELASTIC_STEPS}" ]; then
  echo "Capture horizons do not span two elastic intervals for ${DECK}:" >&2
  echo "  elastic.interval=${ELASTIC_INTERVAL} requires MAX_STEP and TRACE_STEP >= ${MIN_ELASTIC_STEPS}" >&2
  echo "  got MAX_STEP=${MAX_STEP} TRACE_STEP=${TRACE_STEP}" >&2
  echo "Update the per-deck horizon or set explicit admissible overrides; do not profile an elastic-free prefix." >&2
  exit 2
fi

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
CAPTURE_ERRORS=0
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
echo "  deck sha256=$(sha256sum "${INPUT}" | awk '{print $1}')  elastic.interval=${ELASTIC_INTERVAL}"
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
    echo "--- MPI (module/PATH view -- NOT necessarily what the binary runs)"
    which mpiexec srun 2>&1
    echo "MPICH_GPU_SUPPORT_ENABLED=${MPICH_GPU_SUPPORT_ENABLED:-<unset>}"
    echo "MPI_MOD=${MPI_MOD}"
    ompi_info --version 2>&1 | head -1
    # The N4 verdict MUST come from the MPI the binary actually loads. Job
    # 11772152 loaded MPI_MOD=openmpi/4.1.8 and reported its ompi_info, while
    # ldd showed the binary resolving libmpi.so.40 to openmpi-5.0.8 via RPATH
    # -- two different implementations behind the same soname. Ask ldd, then
    # ask THAT tree's ompi_info.
    echo "--- MPI actually linked by ${PLAIN_BIN} (authoritative for N4)"
    LINKED_MPI_LIB=$(ldd "${PLAIN_BIN}" 2>/dev/null | awk '/libmpi\.so/ {print $3; exit}')
    echo "linked_libmpi=${LINKED_MPI_LIB:-<none>}"
    LINKED_OMPI_INFO=""
    if [ -n "${LINKED_MPI_LIB}" ]; then
      LINKED_OMPI_INFO="$(dirname "$(dirname "${LINKED_MPI_LIB}")")/bin/ompi_info"
      echo "linked_ompi_info=${LINKED_OMPI_INFO}"
      [ -x "${LINKED_OMPI_INFO}" ] || echo "  WARNING: no ompi_info beside the linked libmpi; verdict below is the PATH one"
    fi
    [ -x "${LINKED_OMPI_INFO}" ] || LINKED_OMPI_INFO="ompi_info"
    "${LINKED_OMPI_INFO}" --version 2>&1 | head -1
    if ! "${LINKED_OMPI_INFO}" --version 2>/dev/null | head -1 | grep -qF "$(ompi_info --version 2>/dev/null | head -1 | grep -oE 'v[0-9.]+')"; then
      echo "  MISMATCH: the loaded module and the linked library are different Open MPI builds."
      echo "  Trust the linked one. Consider setting MPI_MOD to match the build."
    fi
    echo "--- GPU-aware MPI verdict (Phase 0.5 §2 carry-over), from the LINKED MPI:"
    "${LINKED_OMPI_INFO}" --parsable --all 2>/dev/null | grep -E 'mpi_built_with_cuda_support' || echo "  (no ompi_info cuda line; not Open MPI?)"
    "${LINKED_OMPI_INFO}" --parsable --all 2>/dev/null | grep -iE 'accelerator|smcuda' | head -5
    echo "--- GPU visibility on the BATCH step (srun steps may differ)"
    echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"
    echo "SLURM_JOB_GPUS=${SLURM_JOB_GPUS:-<unset>}  SLURM_STEP_GPUS=${SLURM_STEP_GPUS:-<unset>}"
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
    # Run under srun like every other GPU step. The batch step does not reliably
    # hold the job's GRES -- job 11772151 aborted here with `CUDA error 3 ...
    # initialization error` while nvidia-smi still listed the A100, because
    # nvidia-smi queries without allocating and CUDA init does not.
    srun --mpi=pmix -n1 --gpus-per-task=1 \
      ncu --target-processes all --launch-count 1 --section SpeedOfLight \
        -o "${L}/ncu_probe" --force-overwrite \
        "./${PLAIN_BIN}" "${COMMON[@]}" "${MANAGED[@]}" max_step=1 stop_time=1e99_s \
        plot_file="${L}/probe_plot" > "${L}/ncu_probe.log" 2>&1
    PROBE_RC=$?
    # Three outcomes, not two. The old logic treated "no permission string" as
    # success, so a CUDA-init abort was recorded as `counters permitted` and the
    # ncu leg then ran anyway and failed.
    if grep -qi "ERR_NVGPUCTRPERM\|The user does not have permission" "${L}/ncu_probe.log"; then
      echo "  ncu counters BARRED on this node (ERR_NVGPUCTRPERM) -- ncu leg will be skipped."
      echo "  Substitute: benchmark/res_usage.sh (root-free register/occupancy proxy)."
      NCU_OK=0
    elif [ "${PROBE_RC}" -ne 0 ]; then
      echo "  ncu probe FAILED rc=${PROBE_RC} for a reason other than permissions."
      echo "  Counters are UNVERIFIED; skipping the ncu leg rather than guessing."
      grep -iE "error|abort" "${L}/ncu_probe.log" | head -3
      NCU_OK=0
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
  # PLAN v2 §12.1: startup EXCLUDED, calibration derived up front rather than
  # deferred to analysis. Every prior F10 number here included CUDA context
  # creation and InitData in a 12-33 s total, which is why none of them were
  # admissible. Five repetitions satisfy the campaign protocol; randomized,
  # recorded arm order makes machine drift land on both configurations.
  REPS="${REPS:-5}"
  for mode in managed device; do
    L="${OUT}/timing/${mode}"; mkdir -p "${L}"
    if [ "${mode}" = managed ]; then ARENA=( "${MANAGED[@]}" ); else ARENA=( "${DEVICE[@]}" ); fi
    if [ "${DRYRUN}" = 1 ]; then
      echo "+ [startup calib] srun ... ./${PLAIN_BIN} ${COMMON[*]} ${ARENA[*]} max_step=1 ..."
      continue
    fi
    S0=$(date +%s.%N)
    srun --mpi=pmix -n1 --gpus-per-task=1 \
      "./${PLAIN_BIN}" "${COMMON[@]}" "${ARENA[@]}" max_step=1 \
      stop_time=1e99_s plot_file="${L}/plot_calib" > "${L}/startup.log" 2>&1
    STARTUP_RC=$?
    if [ "${STARTUP_RC}" -ne 0 ]; then
      echo "  startup calibration (${mode}) FAILED rc=${STARTUP_RC}; see ${L}/startup.log"
      CAPTURE_ERRORS=$((CAPTURE_ERRORS + 1))
    fi
    S1=$(date +%s.%N)
    STARTUP=$(awk -v a="${S0}" -v b="${S1}" 'BEGIN{printf "%.3f", b-a}')
    echo "  startup calibration (${mode}, max_step=1): ${STARTUP}s"
    printf 'startup_s=%s\n' "${STARTUP}" > "${L}/startup.txt"
  done
  if [ "${DRYRUN}" = 1 ]; then
    echo "+ [5 randomized/interleaved reps] srun ... managed ... max_step=${MAX_STEP}"
    echo "+ [5 randomized/interleaved reps] srun ... device  ... max_step=${MAX_STEP}"
  else
    for mode in managed device; do
      : > "${OUT}/timing/${mode}/reps.txt"
    done
    : > "${OUT}/timing/order.tsv"
    for rep in $(seq 1 "${REPS}"); do
      RUN_ORDER=$(printf '%s\n' managed device | shuf | paste -sd' ' -)
      printf '%s\t%s\n' "${rep}" "${RUN_ORDER}" >> "${OUT}/timing/order.tsv"
      for mode in ${RUN_ORDER}; do
        L="${OUT}/timing/${mode}"
        if [ "${mode}" = managed ]; then
          ARENA=( "${MANAGED[@]}" )
        else
          ARENA=( "${DEVICE[@]}" )
        fi
        STARTUP=$(awk -F= '/startup_s/{print $2}' "${L}/startup.txt" 2>/dev/null)
        STARTUP="${STARTUP:-0}"
        T0=$(date +%s.%N)
        srun --mpi=pmix -n1 --gpus-per-task=1 \
          "./${PLAIN_BIN}" "${COMMON[@]}" "${ARENA[@]}" max_step="${MAX_STEP}" \
          stop_time=1e99_s plot_file="${L}/plot_rep${rep}" \
          > "${L}/run_rep${rep}.log" 2>&1
        RC=$?; T1=$(date +%s.%N)
        if [ "${RC}" -ne 0 ]; then
          CAPTURE_ERRORS=$((CAPTURE_ERRORS + 1))
        fi
        WALL=$(awk -v a="${T0}" -v b="${T1}" 'BEGIN{printf "%.3f", b-a}')
        # Steady-state per step: startup removed, so this is comparable across
        # decks and across phases. The raw wall is kept for audit.
        PER=$(awk -v w="${WALL}" -v s="${STARTUP}" -v n="${MAX_STEP}" \
              'BEGIN{printf "%.5f", (w-s)/n}')
        echo "  ${mode} rep${rep}: rc=${RC} wall=${WALL}s steady_per_step=${PER}s"
        printf '%s\t%s\t%s\t%s\n' "${rep}" "${RC}" "${WALL}" "${PER}" \
          >> "${L}/reps.txt"
      done
    done
    for mode in managed device; do
      L="${OUT}/timing/${mode}"
      STARTUP=$(awk -F= '/startup_s/{print $2}' "${L}/startup.txt" 2>/dev/null)
      STARTUP="${STARTUP:-0}"
      awk -v m="${mode}" -v n="${MAX_STEP}" -v s="${STARTUP}" '
        {
          if ($2 != 0) {fail++}
          w[NR]=$3; p[NR]=$4; sw+=$3; sp+=$4
        }
        END{
          k=NR; if(k==0){print "no reps"; exit}
          mw=sw/k; mp=sp/k;
          for(i=1;i<=k;i++){vw+=(w[i]-mw)^2; vp+=(p[i]-mp)^2}
          sdw=(k>1)?sqrt(vw/(k-1)):0; sdp=(k>1)?sqrt(vp/(k-1)):0;
          for(i=1;i<=k;i++){ws[i]=w[i];ps[i]=p[i]}
          for(i=2;i<=k;i++){
            x=ws[i];j=i-1;while(j>=1&&ws[j]>x){ws[j+1]=ws[j];j--}ws[j+1]=x
            x=ps[i];j=i-1;while(j>=1&&ps[j]>x){ps[j+1]=ps[j];j--}ps[j+1]=x
          }
          if(k%2){wm=ws[(k+1)/2];pm=ps[(k+1)/2]}
          else{wm=(ws[k/2]+ws[k/2+1])/2;pm=(ps[k/2]+ps[k/2+1])/2}
          printf "mode=%s reps=%d failed_reps=%d steps=%d startup_s=%s\n",
                 m, k, fail, n, s;
          printf "wall_median_s=%.3f steady_per_step_median_s=%.5f\n", wm, pm;
          printf "wall_mean_s=%.3f wall_sd_s=%.3f\n", mw, sdw;
          printf "steady_per_step_mean_s=%.5f steady_per_step_sd_s=%.5f\n", mp, sdp;
        }' "${L}/reps.txt" | tee "${L}/timing.txt"
    done
  fi
  echo "T0 rule: report steady_per_step_median_s; retain mean +/- sd for the"
  echo "  uncertainty gate. A difference inside +/- 2 sd is not a regression"
  echo "  and not an improvement -- PLAN v2 §12.1."
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

# -------------------------------------------------------------- leg: arena ---
# F8 and the allocation-side targets. Nsight's CUDA memory report sees backing
# cudaMalloc/cudaFree calls; AMReX's arenas can recycle those while application
# alloc/free requests churn every step. The profile build's TinyProfiler hooks
# sit at CArena::alloc/free, so their Nalloc/Nfree/CurrentMem/MaxMem tables see
# the requests T4/T5 actually define. A one-step and a full-horizon endpoint
# distinguish initialization from timestep/elastic-solve activity.
if have_leg arena; then
  say "LEG arena -- F8/T1/T4/T5 AMReX request-layer memory profile"
  L="${OUT}/arena"; mkdir -p "${L}"
  ARENA_STEPS="${ARENA_STEPS:-1 ${MAX_STEP}}"
  for arena_step in ${ARENA_STEPS}; do
    A="${L}/step${arena_step}"; mkdir -p "${A}"
    if [ "${DRYRUN}" = 1 ]; then
      echo "+ srun --mpi=pmix -n1 --gpus-per-task=1 ./${PROF_BIN} ${COMMON[*]} ${DEVICE[*]} tiny_profiler.memprof_enabled=1 max_step=${arena_step} ..."
      continue
    fi
    srun --mpi=pmix -n1 --gpus-per-task=1 \
      "./${PROF_BIN}" "${COMMON[@]}" "${DEVICE[@]}" \
      tiny_profiler.memprof_enabled=1 max_step="${arena_step}" \
      stop_time=1e99_s plot_file="${A}/plot" > "${A}/run.log" 2>&1
    ARENA_RC=$?
    echo "  max_step=${arena_step} rc=${ARENA_RC} -> ${A}/run.log"
    if [ "${ARENA_RC}" -ne 0 ] || ! grep -q "Device Memory Usage:" "${A}/run.log"; then
      echo "  ARENA PROFILE FAILED: request-layer table missing" \
        | tee "${A}/CAPTURE_FAILED"
      CAPTURE_ERRORS=$((CAPTURE_ERRORS + 1))
    fi
  done
  echo "  Compare endpoint Nalloc-Nfree and CurrentMem for T5b; normalize the"
  echo "  Nalloc delta by added steps/elastic solves for T4. Treat only"
  echo "  non-initialization Managed Memory rows as T1 application allocations."
fi

# --------------------------------------------------------------- leg: nsys ---
if have_leg nsys; then
  say "LEG nsys -- F1,F2,F3,F9 (profile binary: BL_PROFILE regions emit NVTX)"
  L="${OUT}/nsys"; mkdir -p "${L}"
  # NO `mpi` OR `osrt` IN -t. Nsight Systems 2024.6.2's OSRT injection aborts
  # inside ompi_mpi_init -> opal_common_ucx_mca_register against the Open MPI
  # 5.0.8 / UCX 1.18.1 this binary links. Isolation job 11825250 proved that
  # `cuda`, `nvtx`, and `cuda,nvtx` all complete, while `osrt` alone and
  # `cuda,nvtx,osrt` both abort with "Expected shared object name, found a path
  # delimiter". Dropping `srun` did not help (job 11825224). Every run here is
  # -n1, and CUDA tracing still supplies the API calls needed by F2/T6a.
  # Do not re-add mpi/OSRT without first re-testing the linked MPI and nsys.
  NSYS_ARGS=( profile
    -t cuda,nvtx
    --cuda-memory-usage=true
    --cuda-um-cpu-page-faults=true
    --cuda-um-gpu-page-faults=true
    --force-overwrite true
    -o "${L}/trace" )
  if [ "${DRYRUN}" = 1 ]; then
    echo "+ srun --mpi=pmix -n1 --gpus-per-task=1 nsys ${NSYS_ARGS[*]} ./${PROF_BIN} ${COMMON[*]} ${MANAGED[*]} max_step=${TRACE_STEP} ..."
    echo "+ nsys stats --report nvtx_sum,nvtx_gpu_proj_sum,nvtx_kern_sum,cuda_gpu_kern_sum,cuda_gpu_mem_size_sum,cuda_gpu_mem_time_sum,cuda_api_sum --format csv ${L}/trace.nsys-rep"
  else
    srun --mpi=pmix -n1 --gpus-per-task=1 \
      nsys "${NSYS_ARGS[@]}" \
      "./${PROF_BIN}" "${COMMON[@]}" "${MANAGED[@]}" max_step="${TRACE_STEP}" \
      stop_time=1e99_s plot_file="${L}/plot" > "${L}/run.log" 2>&1
    NSYS_RC=$?
    echo "  nsys rc=${NSYS_RC}"
    if [ "${NSYS_RC}" -ne 0 ]; then
      CAPTURE_ERRORS=$((CAPTURE_ERRORS + 1))
    fi
    # Distinguish "nsys produced nothing" from "this nsys lacks that report".
    # The old loop reported eight MISSING reports when the real fault was an
    # absent trace, which reads as a version quirk instead of a dead capture.
    # Test SIZE, not existence. nsys creates trace.nsys-rep up front and starts
    # streaming into it, so a capture killed mid-write leaves a 0-byte file that
    # passes -f. That is exactly what happened on the 3D deck in 11775090: the
    # guard stayed quiet, `ls` showed a full set of report files, and every one
    # of them was empty. A dead capture that looks alive is the failure this
    # guard exists to prevent.
    if [ ! -s "${L}/trace.nsys-rep" ]; then
      echo "  CAPTURE FAILED: no ${L}/trace.nsys-rep (nsys rc=${NSYS_RC})."
      echo "  F1,F2,F3,F9 are NOT captured for this deck. Read ${L}/run.log;"
      echo "  An MPI-init abort with a libToolsInjection64.so frame means a trace"
      echo "  feature is incompatible with the linked MPI, not that the deck ran."
      echo "  Do not read the stats CSVs -- there are none."
      echo "no trace.nsys-rep produced; nsys rc=${NSYS_RC}" > "${L}/CAPTURE_FAILED"
    else
      # Export CSVs so the collect step moves kilobytes, not gigabytes. Report names
      # vary across nsys versions -- try each, record which ones existed.
      for rep in nvtx_sum nvtx_gpu_proj_sum nvtx_kern_sum cuda_gpu_kern_sum \
                 cuda_gpu_mem_size_sum cuda_gpu_mem_time_sum cuda_api_sum \
                 cuda_gpu_trace um_total_sum um_cpu_page_faults_sum; do
        nsys stats --report "${rep}" --format csv --output "${L}/${rep}" \
             "${L}/trace.nsys-rep" >> "${L}/stats.log" 2>&1
        # nsys stats exits 0 having written a 0-byte CSV when the sqlite export
        # fails ("Exportation error: Cannot read from stream"). Judge by the
        # output, not the exit code.
        OUTCSV=$(ls "${L}/${rep}"*.csv 2>/dev/null | head -1)
        if [ -s "${OUTCSV:-/nonexistent}" ]; then
          echo "  stats ok: ${rep}"
      else
        echo "  stats EMPTY/UNAVAILABLE: ${rep}" | tee -a "${L}/stats_missing.txt"
        fi
      done
      KERNEL_SUMMARY=$(ls "${L}/cuda_gpu_kern_sum"*.csv 2>/dev/null | head -1)
      if [ -s "${KERNEL_SUMMARY:-/nonexistent}" ]; then
        PYTHONDONTWRITEBYTECODE=1 python3 benchmark/nsys_discover_targets.py \
          "${KERNEL_SUMMARY}" --top 10 --ncu-targets "${NCU_TARGETS}" \
          --top-output "${L}/discovered_top10.tsv" \
          --target-output "${L}/discovered_ncu_targets.tsv" \
          | tee "${L}/discovery.log"
      else
        echo "kernel discovery unavailable: cuda_gpu_kern_sum is empty" \
          | tee "${L}/DISCOVERY_FAILED"
        CAPTURE_ERRORS=$((CAPTURE_ERRORS + 1))
      fi
    fi
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
    echo "  SKIPPED: counters unavailable or unverified on this node (see env leg)."
    echo "  Use benchmark/res_usage.sh as the root-free substitute."
  else
    # The nsys leg ranks raw kernels by total GPU time and derives demangled-name
    # selectors for the first distinct functions. This is load-bearing:
    # Fapply is nested under Fsmooth and mgVcycle, so NVTX range filters profile
    # Fapply three times while never measuring the actual Fsmooth kernel.
    SECTIONS=( --section SpeedOfLight
               --section SpeedOfLight_RooflineChart
               --section MemoryWorkloadAnalysis
               --section Occupancy
               --section LaunchStats )
    TARGETS_FILE="${OUT}/nsys/discovered_ncu_targets.tsv"
    if [ "${DRYRUN}" = 1 ]; then
      echo "+ python3 benchmark/nsys_discover_targets.py <cuda_gpu_kern_sum.csv> --top 10 --ncu-targets ${NCU_TARGETS} ..."
      echo "+ for each discovered selector: srun ... ncu --kernel-name-base demangled --kernel-name <regex> --launch-count ${NCU_LAUNCHES} --kill yes ..."
    elif [ ! -s "${TARGETS_FILE}" ]; then
      echo "  SKIPPED: discovered targets missing (${TARGETS_FILE}); NCU cannot claim F4." \
        | tee "${L}/DISCOVERY_FAILED"
      CAPTURE_ERRORS=$((CAPTURE_ERRORS + 1))
    else
      while IFS=$'\t' read -r rank time_percent total_time_ms instances label selector; do
        [ "${rank}" = "rank" ] && continue
        tag=$(printf 'rank%s_%s' "${rank}" "${label}" | tr -cd '[:alnum:]_')
        # --kill yes is load-bearing, not tidiness. Without it ncu collects its
        # NCU_LAUNCHES matches at the FIRST elastic solve and then lets the app
        # run the whole horizon with NVTX collection still attached, retaining
        # range state for every BL_PROFILE push/pop (~16k Fapply ranges per 2D
        # step). Job 11772161 reached step 81 that way and was OOM-killed at
        # 64 GB host RAM on MLMG::mgVcycle(); the 3D job 11772153 lost all
        # three ranges the same way. The two survivors still wrote 173 MB
        # reports for 3 profiled launches.
        srun --mpi=pmix -n1 --gpus-per-task=1 \
          ncu --kernel-name-base demangled --kernel-name "${selector}" \
              --launch-count "${NCU_LAUNCHES}" \
              --kill yes \
              "${SECTIONS[@]}" -o "${L}/${tag}" --force-overwrite \
          "./${PROF_BIN}" "${COMMON[@]}" "${MANAGED[@]}" max_step="${TRACE_STEP}" \
          stop_time=1e99_s plot_file="${L}/plot_${tag}" > "${L}/${tag}.log" 2>&1
        # rc is EXPECTED to be nonzero under --kill yes: ncu terminates the app
        # once it has its launches, so judge this leg by the .ncu-rep, not rc.
        NCU_RC=$?
        echo "  rank=${rank} ${label} (${time_percent}% GPU, ${selector}) rc=${NCU_RC} -> ${L}/${tag}.ncu-rep"
        if [ -s "${L}/${tag}.ncu-rep" ]; then
          ncu --import "${L}/${tag}.ncu-rep" --csv --page details \
            > "${L}/${tag}.csv" 2>>"${L}/${tag}.log" && echo "    csv exported"
        else
          echo "rank=${rank} label=${label} selector=${selector} produced no report" \
            | tee "${L}/${tag}.CAPTURE_FAILED"
          CAPTURE_ERRORS=$((CAPTURE_ERRORS + 1))
        fi
      done < "${TARGETS_FILE}"
    fi
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
echo "required_leg_failures=${CAPTURE_ERRORS}"
if [ "${CAPTURE_ERRORS}" -ne 0 ]; then
  echo "phase0_capture FAILED: ${CAPTURE_ERRORS} required simulation/profile legs failed."
  exit 4
fi
