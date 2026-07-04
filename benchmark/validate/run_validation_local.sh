#!/usr/bin/env bash
# Thin entry point for the local-A1000 validation runner (roadmap task 1.B).
# The orchestration lives in run_validation_local.py (manifest parsing,
# subprocess/file handling -- safer in Python than in bash string-soup); this
# wrapper just gives it a stable `bash benchmark/validate/run_validation_local.sh`
# invocation, matching the naming convention of the NOVA counterpart
# (run_validation_nova.slurm). All env vars / flags pass straight through, e.g.:
#
#   bash benchmark/validate/run_validation_local.sh --profiles gpu_strict --case canonical_2d_elastic
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
exec python3 "$ROOT/benchmark/validate/run_validation_local.py" "$@"
