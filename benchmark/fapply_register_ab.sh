#!/usr/bin/env bash
# Static register/stack A/B for the elastic hot kernels (Fapply / Diagonal).
#
# Usage:
#   bash benchmark/fapply_register_ab.sh <binary-or-object> [<binary-or-object-B>]
#
# Examples:
#   bash benchmark/fapply_register_ab.sh bin/alamo_gpu-3d-cuda86-g++
#   bash benchmark/fapply_register_ab.sh \
#       obj/obj-3d-cuda86-g++/Operator/Elastic.cpp.o \
#       /home/jackplum/Projects/alamo-elastic-opt/obj/obj-3d-cuda86-g++/Operator/Elastic.cpp.o
#
# Why this exists (see benchmark/GPU_STRUCTURAL_PLAN_20260703.md §2.1):
# the 3D Sym::Major Fapply kernel's register pressure is visible *statically*
# and *locally* (sm_86 cuobjdump matches the A100 ncu 255-reg finding), so
# register-pressure edits can be A/B'd in minutes on the A1000 workstation —
# no NOVA round-trip. Wall-clock claims still require an A100 run.
#
# Sym enum mapping (src/Set/Matrix4.H): 0=None 1=Major 2=Minor 3=MajorMinor
# 4=Diagonal 5=Full 6=Isotropic. The chamber sims (NeoHookeanPredeformed)
# instantiate Elastic<1> (Sym::Major) — that is the row that matters.
set -euo pipefail
cd "$(dirname "$0")/.."

CUOBJDUMP="${CUOBJDUMP:-}"
if [ -z "${CUOBJDUMP}" ]; then
    for c in cuobjdump .local/cuda/bin/cuobjdump .local/cuda-12.6.3-redist/bin/cuobjdump; do
        if command -v "$c" >/dev/null 2>&1; then CUOBJDUMP="$c"; break; fi
        if [ -x "$c" ]; then CUOBJDUMP="$(pwd)/$c"; break; fi
    done
fi
if [ -z "${CUOBJDUMP}" ]; then
    echo "cuobjdump not found; set CUOBJDUMP=<path>" >&2
    exit 1
fi

report() {
    local target="$1"
    echo "== ${target}"
    # Main ParallelFor kernels only (the hot per-node loops); reduce/probe
    # kernels are filtered out. One line per (class<SYM>, method).
    "${CUOBJDUMP}" --dump-resource-usage "${target}" 2>/dev/null \
    | awk '
        /Function _ZN5amrex13launch_globalILi256EZNS_11ParallelForILi256E/ {
            name=$2
            sym=""; meth=""
            if (match(name, /ElasticILi[0-9]+EE/)) {
                sym=substr(name, RSTART+10, RLENGTH-12)
                rest=substr(name, RSTART+RLENGTH)
                if (rest ~ /^6Fapply/)   meth="Fapply"
                if (rest ~ /^8Diagonal/) meth="Diagonal"
                if (rest ~ /^6Stress/)   meth="Stress"
                if (rest ~ /^6Energy/)   meth="Energy"
            }
            if (meth != "") { getline; printf "  Elastic<%s>::%s  %s %s\n", sym, meth, $1, $2 }
        }' \
    | sort -u
}

report "$1"
if [ $# -ge 2 ]; then
    echo
    report "$2"
    echo
    echo "(compare the Elastic<1>::Fapply rows: REG and STACK; STACK>0 on the"
    echo " main ParallelFor means spill-adjacent local-memory frame)"
fi
