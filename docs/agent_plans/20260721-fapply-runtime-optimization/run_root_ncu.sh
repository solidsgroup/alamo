#!/usr/bin/env bash
# One-shot privileged 3D application-replay capture for the frozen Step 1 baseline.
# Run manually as:
#   sudo bash docs/agent_plans/20260721-fapply-runtime-optimization/run_root_ncu.sh
set -euo pipefail

if [[ "${EUID}" -ne 0 ]]; then
    echo "ERROR: this script must be run with sudo/root" >&2
    exit 1
fi

REPO=/home/jackplum/Projects/alamo
TASK_ROOT="${REPO}/docs/agent_plans/20260721-fapply-runtime-optimization"
BASE="${TASK_ROOT}/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline"
OUT="${BASE}/ncu_retry3"
RAW_LOG="${BASE}/raw/ncu_fapply_root_retry3.log"
NCU="${REPO}/.local/nsight-compute/usr/lib/nsight-compute/ncu"
BIN_3D="${BASE}/bin/alamo_gpu-3d-fast-sm86-baseline"
KERNEL_3D='_ZN5amrex13launch_globalILi256EZNS_11ParallelForILi256EZNK8Operator7ElasticILi1EE6FapplyEiiRNS_8MultiFabERKS5_EUliiiE_Li3EEENSt9enable_ifIXsr5amrex19MaybeDeviceRunnableIT0_vEE5valueEvE4typeERKNS_3Gpu10KernelInfoERKNS_5BoxNDIXT1_EEERKSB_EUlvE_EEvSB_'

cd "${REPO}"
test -x "${NCU}"
test -x "${BIN_3D}"
test ! -e "${OUT}"
test ! -e "${RAW_LOG}"
mkdir -p "${OUT}/work/3d"

restore_owner() {
    if [[ -n "${SUDO_UID:-}" && -n "${SUDO_GID:-}" ]]; then
        chown -R "${SUDO_UID}:${SUDO_GID}" "${OUT}" "${RAW_LOG}" 2>/dev/null || true
    fi
}
trap restore_owner EXIT

exec > >(tee "${RAW_LOG}") 2>&1
set -x

"${NCU}" --version
"${NCU}" --query-metrics > "${OUT}/query_metrics.txt"
"${NCU}" --list-sections > "${OUT}/sections.txt"
sha256sum "${NCU}" "${BIN_3D}" input_3d_centre_bore_128_a2
nvidia-smi \
    --query-gpu=name,uuid,driver_version,pstate,clocks.current.sm,clocks.current.memory,power.limit,memory.used \
    --format=csv,noheader
nvidia-smi --query-compute-apps=pid,process_name,used_gpu_memory \
    --format=csv,noheader

COMMON_NCU=(
    --target-processes all
    --launch-count 1
    --kill yes
    --replay-mode application
    --section LaunchStats
    --section Occupancy
    --metrics gpu__time_duration.sum,smsp__inst_executed_op_local_ld.sum,smsp__inst_executed_op_local_st.sum,l1tex__t_bytes_pipe_lsu_mem_local_op_ld.sum,l1tex__t_bytes_pipe_lsu_mem_local_op_st.sum
    --force-overwrite
)

"${NCU}" "${COMMON_NCU[@]}" --kernel-name-base mangled \
    --kernel-name "${KERNEL_3D}" --export "${OUT}/3d_fapply" \
    "${BIN_3D}" input_3d_centre_bore_128_a2 \
    max_step=2 stop_time=1e99_s elastic.interval=1 \
    amr.max_grid_size=64 amrex.the_arena_is_managed=1 \
    elastic.solver.nriters=2 elastic.print_model=0 \
    elastic.solver.verbose=0 elastic.solver.nr_diagnostics=0 \
    amr.plot_int=-1 amr.thermo.plot_int=-1 \
    plot_file="${OUT}/work/3d/plot"
test -s "${OUT}/3d_fapply.ncu-rep"

"${NCU}" --import "${OUT}/3d_fapply.ncu-rep" --page raw --csv \
    > "${OUT}/3d_fapply_raw.csv"

nvidia-smi \
    --query-gpu=name,uuid,driver_version,pstate,clocks.current.sm,clocks.current.memory,power.limit,memory.used \
    --format=csv,noheader
set +x
echo "Root Nsight Compute captures completed under ${OUT}"
