#!/usr/bin/env bash
set -euo pipefail

if [ -z "${ALAMO_BINARY:-}" ]; then
    echo "ALAMO_BINARY must name the real Alamo executable" >&2
    exit 2
fi

exec "$ALAMO_BINARY" "$@" \
    elastic.solver.no_gpu_sync=1 \
    amrex.the_arena_init_size="${ARENA_INIT_SIZE:-1073741824}"
