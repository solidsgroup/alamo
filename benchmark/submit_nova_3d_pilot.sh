#!/usr/bin/env bash
# Submit the 3D pilot matrix: 2 decks x 5 hardware configs = 10 jobs.
#
#   config   ranks  resources
#   cpu64      64   1 node, 64 cores
#   cpu128    128   2 nodes, 64 cores each
#   a100x1      1   1x A100
#   a100x2      2   2x A100 (1 rank/GPU; NOTE: multi-GPU was a LOSS in the
#                   Phase-3 crossover study -- this arm re-tests that on the
#                   real decks, expect it may lose to a100x1)
#   h200x1      1   1x H200
#
# Dry-run by default; pass --submit to actually sbatch.
# Build both CPU (dim=3) and GPU (cuda80 + cuda90) 3D binaries first.

set -euo pipefail
cd "$(dirname "$0")/.."

DECKS=(input_nova_3d_pilot_rod_and_tube input_nova_3d_pilot_cross)
SLURM=benchmark/nova_3d_pilot.slurm
# Optional: DEP=<jobid> gates every job on a build job (afterok).
DEPFLAG=()
[[ -n "${DEP:-}" ]] && DEPFLAG=(--dependency=afterok:"${DEP}")
DO=echo
[[ "${1:-}" == "--submit" ]] && DO=""

for INPUT in "${DECKS[@]}"; do
    [[ -f "${INPUT}" ]] || { echo "missing deck ${INPUT}" >&2; exit 1; }

    # -- CPU ---------------------------------------------------------------
    ${DO} env INPUT="${INPUT}" BACKEND=cpu TAG=cpu64 \
        sbatch "${DEPFLAG[@]}" --nodes=1 --ntasks=64 --cpus-per-task=1 --mem=180G \
        --export=ALL,INPUT="${INPUT}",BACKEND=cpu,TAG=cpu64 "${SLURM}"

    ${DO} env INPUT="${INPUT}" BACKEND=cpu TAG=cpu128 \
        sbatch "${DEPFLAG[@]}" --nodes=2 --ntasks=128 --ntasks-per-node=64 --cpus-per-task=1 --mem=180G \
        --export=ALL,INPUT="${INPUT}",BACKEND=cpu,TAG=cpu128 "${SLURM}"

    # -- GPU (1 rank per GPU) ---------------------------------------------
    ${DO} env sbatch "${DEPFLAG[@]}" --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=120G \
        --gres=gpu:a100:1 \
        --export=ALL,INPUT="${INPUT}",BACKEND=gpu,ARCH=80,TAG=a100x1 "${SLURM}"

    ${DO} env sbatch "${DEPFLAG[@]}" --nodes=1 --ntasks=2 --cpus-per-task=8 --mem=180G \
        --gres=gpu:a100:2 \
        --export=ALL,INPUT="${INPUT}",BACKEND=gpu,ARCH=80,TAG=a100x2 "${SLURM}"

    ${DO} env sbatch "${DEPFLAG[@]}" --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=120G \
        --gres=gpu:h200:1 \
        --export=ALL,INPUT="${INPUT}",BACKEND=gpu,ARCH=90,TAG=h200x1 "${SLURM}"
done

echo "---"
echo "Analysis after runs finish: for each pilot3d.<jobid>.out, take"
echo "  (1) steady cost/step from the last ~50 step lines (post-ignition),"
echo "  (2) elastic-solve wall from the tiny-profiler / verbose=4 blocks."
echo "Full-length (L=2D) estimate = 4 x slab cost/step x 20000 steps (5 s)."
echo "Cross deck may OOM on a100x1 (16.7M base cells, nodal elastic) --"
echo "an OOM there is itself the answer for single-GPU feasibility."
