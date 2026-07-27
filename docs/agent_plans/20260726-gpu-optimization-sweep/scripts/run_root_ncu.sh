#!/usr/bin/env bash
# Nsight Compute capture of the hot Fapply kernel.  Needs root: this driver has
# NVreg_RestrictProfilingToAdminUsers at its default, so an unprivileged ncu
# fails with ERR_NVGPUCTRPERM.
#
#   sudo bash docs/agent_plans/20260726-gpu-optimization-sweep/scripts/run_root_ncu.sh [2d|3d]
#
# Requires an otherwise idle GPU (the 3D arm needs the whole card).
set -euo pipefail

if [[ "${EUID}" -ne 0 ]]; then
    echo "ERROR: run with sudo" >&2
    exit 1
fi

REGIME="${1:-2d}"
REPO=/home/jackplum/Projects/alamo
TASK="${REPO}/docs/agent_plans/20260726-gpu-optimization-sweep"
OUT="${TASK}/artifacts/ncu"
NCU="${REPO}/.local/nsight-compute/usr/bin/ncu"

cd "${REPO}"
mkdir -p "${OUT}/work"

restore_owner() {
    if [[ -n "${SUDO_UID:-}" && -n "${SUDO_GID:-}" ]]; then
        chown -R "${SUDO_UID}:${SUDO_GID}" "${OUT}" 2>/dev/null || true
    fi
}
trap restore_owner EXIT

if [[ "${REGIME}" == 2d ]]; then
    BIN="${REPO}/bin/alamo_gpu-2d-profile-cuda86-g++"
    INPUT=tests/ElasticSoftVoid/input
    ARGS=(
        max_step=2 amr.max_level=2 amr.n_cell=64 64 8
        explicitmesh.lo1=32 32 0 explicitmesh.hi1=95 95 0
        explicitmesh.lo2=72 72 0 explicitmesh.hi2=183 183 0
        pf.eta.ic.expression.constant.w=0.002
        model_void.kappa=0.2_MPa model_void.mu=0.2_MPa
    )
    SKIP=400
else
    BIN="${REPO}/bin/alamo_gpu-3d-profile-cuda86-g++"
    INPUT=input_3d_centre_bore_128_a2
    ARGS=(
        max_step=2 stop_time=1e99_s elastic.interval=1
        amr.max_grid_size=64 amrex.the_arena_is_managed=1
        elastic.solver.nriters=2
    )
    SKIP=200
fi

test -x "${NCU}"
test -x "${BIN}"

set -x
"${NCU}" --target-processes all \
    --launch-count 1 --launch-skip "${SKIP}" --kill yes \
    --replay-mode kernel --kernel-name-base mangled --kernel-name 'regex:6Fapply' \
    --section SpeedOfLight --section Occupancy --section LaunchStats \
    --section MemoryWorkloadAnalysis --section WarpStateStats \
    --section SchedulerStats --section InstructionStats \
    --force-overwrite --export "${OUT}/${REGIME}_fapply" \
    "${BIN}" "${INPUT}" "${ARGS[@]}" \
    elastic.print_model=0 elastic.solver.verbose=0 \
    elastic.solver.nr_diagnostics=0 amr.plot_int=-1 amr.thermo.plot_int=-1 \
    plot_file="${OUT}/work/${REGIME}_plot"
set +x

"${NCU}" --import "${OUT}/${REGIME}_fapply.ncu-rep" --page details \
    > "${OUT}/${REGIME}_fapply_details.txt"
echo "wrote ${OUT}/${REGIME}_fapply.ncu-rep and _details.txt"
