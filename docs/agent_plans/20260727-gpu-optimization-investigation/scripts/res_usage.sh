#!/usr/bin/env bash
# Static device-resource report for the elastic kernels of an alamo GPU binary.
#
#   res_usage.sh <binary> [label]
#
# Uses `cuobjdump -res-usage`, which needs no elevated privileges -- unlike
# Nsight Compute, which fails with ERR_NVGPUCTRPERM on this workstation.
# Reports registers, per-thread stack (local-memory spill) and the implied
# register-limited occupancy on sm_86 for the production Sym::Major (=1)
# instantiation of Operator::Elastic.
#
# Set::Sym enum (src/Set/Matrix4.H:12):
#   None=0 Major=1 Minor=2 MajorMinor=3 Diagonal=4 Full=5 Isotropic=6
# Production models (Flame and tests/ElasticSoftVoid) use
# NeoHookeanPredeformed : NeoHookean : Solid<Set::Sym::Major>, i.e. Elastic<1>.
set -euo pipefail

BIN="${1:?usage: res_usage.sh <binary> [label]}"
LABEL="${2:-$(basename "$BIN")}"

REPO=/home/jackplum/Projects/alamo
# shellcheck source=/dev/null
source "$REPO/benchmark/local_cuda_env.sh"

# sm_86: 65536 32-bit registers per SM, 1536 max resident threads per SM.
REGS_PER_SM=65536
MAX_THREADS_PER_SM=1536

echo "== elastic kernel resources: $LABEL =="
printf '%-28s %-10s %6s %8s %10s %12s\n' kernel sym regs stack_B thr_per_sm occupancy

cuobjdump -res-usage "$BIN" \
| awk -v regs_per_sm="$REGS_PER_SM" -v max_thr="$MAX_THREADS_PER_SM" '
  /^ Function / { fn = $2; next }
  /REG:/ && fn != "" {
    reg = ""; stack = ""
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^REG:/)   { sub(/^REG:/,   "", $i); reg   = $i }
      if ($i ~ /^STACK:/) { sub(/^STACK:/, "", $i); stack = $i }
    }
    # Demangled-name matching on the Operator::Elastic<N>::<Kernel> signature.
    # The mangled form embeds "8Operator7ElasticILi<N>EE<len><name>".
    if (match(fn, /8Operator7ElasticILi[0-9]+EE/)) {
      whole  = substr(fn, RSTART, RLENGTH)
      rest   = substr(fn, RSTART + RLENGTH)
      # Pull N out of the "ILi<N>EE" template-argument token only, so the
      # "8" and "7" name-length prefixes of "8Operator7Elastic" are not
      # swept into the symmetry id.
      match(whole, /ILi[0-9]+EE/)
      symtok = substr(whole, RSTART + 3, RLENGTH - 5)
      # rest starts with <len><name>; pull the identifier out.
      if (match(rest, /^[0-9]+/)) {
        n    = substr(rest, RSTART, RLENGTH) + 0
        name = substr(rest, RSTART + RLENGTH, n)
        if (name == "Fapply" || name == "Diagonal" || name == "Fsmooth" ||
            name == "Energy" || name == "Stress"   || name == "Normalize") {
          thr = int(regs_per_sm / reg)
          if (thr > max_thr) thr = max_thr
          printf "%-28s %-10s %6s %8s %10d %11.1f%%\n", \
                 name, (symtok == "1" ? "Major*" : "Sym=" symtok), reg, stack, \
                 thr, 100.0 * thr / max_thr
        }
      }
    }
    fn = ""
  }
' | sort -k2,2 -k1,1

echo "(* = production instantiation.  occupancy is the register-limited"
echo " theoretical ceiling on sm_86, not an achieved measurement.)"
