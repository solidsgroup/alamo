#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

DIMS="${DIMS:-2 3}"
COMP="${COMP:-clang++}"
JOBS="${JOBS:-8}"
TIMEOUT="${TIMEOUT:-10000}"
TESTS="${TESTS:-}"
OUT_DIR="${OUT_DIR:-benchmark/phase4_cpu_semantics_$(date +%Y%m%d_%H%M%S)}"
CONFIGURE_EXTRA="${CONFIGURE_EXTRA:-}"
RUNTESTS_EXTRA="${RUNTESTS_EXTRA:-}"
MPIRUN_FLAGS="${MPIRUN_FLAGS:-}"
NO_BUILD="${NO_BUILD:-0}"
DRY_RUN="${DRY_RUN:-0}"
PYTHONNOUSERSITE="${PYTHONNOUSERSITE:-1}"
MPLBACKEND="${MPLBACKEND:-Agg}"

export PYTHONNOUSERSITE
export MPLBACKEND

usage() {
    cat <<'USAGE'
Run the Phase 4.1 CPU-semantics regression.

Environment overrides:
  DIMS="2 3"                 dimensions to configure/build/test
  COMP=clang++               compiler name passed to configure and runtests.py
  JOBS=8                     make parallelism
  TIMEOUT=10000              per-test timeout for scripts/runtests.py
  TESTS=""                   tests path/filter; empty runs full ./tests/* inventory
  OUT_DIR=benchmark/...      log directory
  CONFIGURE_EXTRA="..."      extra ./configure arguments
  RUNTESTS_EXTRA="..."       extra scripts/runtests.py arguments
  MPIRUN_FLAGS="..."         flags passed through --mpirun-flags
  NO_BUILD=1                 skip configure and make
  DRY_RUN=1                  list tests without running them
  PYTHONNOUSERSITE=1         hide user-site packages for Python checks
  MPLBACKEND=Agg             non-interactive matplotlib backend

Example:
  DIMS="2 3" MPIRUN_FLAGS="--oversubscribe" benchmark/phase4_cpu_semantics_regression.sh
  TESTS="tests/Unit" benchmark/phase4_cpu_semantics_regression.sh
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

mkdir -p "${OUT_DIR}"

echo "============================================================"
echo " Phase 4.1 CPU semantics regression"
echo " repo            = ${ROOT_DIR}"
echo " dims            = ${DIMS}"
echo " compiler        = ${COMP}"
echo " tests           = ${TESTS:-<all>}"
echo " out_dir         = ${OUT_DIR}"
echo " configure_extra = ${CONFIGURE_EXTRA}"
echo " runtests_extra  = ${RUNTESTS_EXTRA}"
echo " mpirun_flags    = ${MPIRUN_FLAGS}"
echo " no_build        = ${NO_BUILD}"
echo " dry_run         = ${DRY_RUN}"
echo " PYTHONNOUSERSITE= ${PYTHONNOUSERSITE}"
echo " MPLBACKEND      = ${MPLBACKEND}"
echo "============================================================"
echo

for dim in ${DIMS}; do
    log="${OUT_DIR}/dim${dim}.log"
    echo "=== Phase 4.1 dim=${dim} ===" | tee "${log}"

    if [[ "${NO_BUILD}" != "1" ]]; then
        echo "+ ./configure --dim=${dim} --comp=${COMP} ${CONFIGURE_EXTRA}" | tee -a "${log}"
        # shellcheck disable=SC2086
        ./configure --dim="${dim}" --comp="${COMP}" ${CONFIGURE_EXTRA} 2>&1 | tee -a "${log}"

        echo "+ make -j${JOBS}" | tee -a "${log}"
        make -j"${JOBS}" 2>&1 | tee -a "${log}"
    fi

    cmd=(scripts/runtests.py "--dim=${dim}" "--comp=${COMP}" "--timeout=${TIMEOUT}" "--no-backspace")
    if [[ -n "${MPIRUN_FLAGS}" ]]; then
        cmd+=("--mpirun-flags=${MPIRUN_FLAGS}")
    fi
    if [[ "${DRY_RUN}" == "1" ]]; then
        cmd+=("--dryrun" "--cmd")
    fi
    if [[ -n "${RUNTESTS_EXTRA}" ]]; then
        # Keep this as env-driven shell splitting so callers can pass native
        # runtests.py flags without this wrapper mirroring every option.
        # shellcheck disable=SC2206
        extra_args=(${RUNTESTS_EXTRA})
        cmd+=("${extra_args[@]}")
    fi
    if [[ -n "${TESTS}" ]]; then
        # shellcheck disable=SC2206
        test_args=(${TESTS})
        cmd+=("${test_args[@]}")
    fi

    printf '+ %q' "${cmd[@]}" | tee -a "${log}"
    echo | tee -a "${log}"
    "${cmd[@]}" 2>&1 | tee -a "${log}"
done
