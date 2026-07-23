#!/usr/bin/env bash
#SBATCH --account=brunnels
#SBATCH --partition=interactive
#SBATCH --job-name=fapply_a100_ncu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:a100:1
#SBATCH --time=00:20:00
#SBATCH --output=nova-a100-ncu.%j.out
#SBATCH --error=nova-a100-ncu.%j.err

set -euo pipefail

REPO=/work/brunnels/jackplum/alamo-fapply-20260722
BIN="$REPO/bin/alamo_gpu-3d-profile-cuda80-g++"
INPUT="$REPO/input_3d_centre_bore_128_a2"
OUT="$REPO/nova-results-11751561/ncu"
EXPECTED_BIN_SHA=036c7ec520c6ae301b31d1e75a107b5dff2e343f02729c880f7be3ceeda98c84

cd "$REPO"
test "$(sha256sum "$BIN" | awk '{print $1}')" = "$EXPECTED_BIN_SHA"
test ! -e "$OUT"
mkdir -p "$OUT"

module purge 2>/dev/null || true
module load cuda 2>/dev/null || module load cuda/12
module load gcc 2>/dev/null || module load gcc/12
module load openmpi 2>/dev/null || module load openmpi4
nvidia-smi -L
ncu --version

COMMON=(
  max_step=2 stop_time=1e99_s elastic.interval=1
  amr.max_grid_size=64 amrex.the_arena_is_managed=1
  elastic.solver.nriters=2 elastic.print_model=0
  elastic.solver.verbose=0 elastic.solver.nr_diagnostics=0
  amr.plot_int=-1 amr.thermo.plot_int=-1
  elastic.solver.pre_smooth=2 elastic.solver.post_smooth=2
)

set +e
ncu --section LaunchStats --section Occupancy --replay-mode application \
  --kernel-name-base demangled \
  --kernel-name 'regex:.*Elastic.*Fapply.*' \
  --launch-skip 100 --launch-count 1 --force-overwrite \
  --export "$OUT/fapply_2x2" \
  "$BIN" "$INPUT" "${COMMON[@]}" plot_file="$OUT/plot" \
  >"$OUT/run.stdout" 2>"$OUT/run.stderr"
ncu_rc=$?
set -e

if (( ncu_rc != 0 )); then
  if grep -q ERR_NVGPUCTRPERM "$OUT/run.stderr"; then
    echo "NOVA_A100_NCU_COUNTERS_BLOCKED ERR_NVGPUCTRPERM"
    exit 0
  fi
  echo "NOVA_A100_NCU_FAILED rc=$ncu_rc" >&2
  exit "$ncu_rc"
fi

ncu --import "$OUT/fapply_2x2.ncu-rep" --page details --csv \
  --log-file "$OUT/fapply_2x2_details.csv"
sha256sum "$OUT/fapply_2x2.ncu-rep" "$OUT/fapply_2x2_details.csv" \
  >"$OUT/SHA256SUMS"
echo NOVA_A100_NCU_PASS
