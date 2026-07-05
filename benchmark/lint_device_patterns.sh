#!/usr/bin/env bash
# benchmark/lint_device_patterns.sh -- grep-based enforcement of the three
# device-side bug classes documented in docs/llm/BUG_PATTERNS.md (REMEDIATION_PLAN
# Phase 3). Each rule below is derived from, and calibrated against, an actual
# historical pre-fix/post-fix pair in this repo (commit hashes in each section).
# Scope: src/Operator src/Integrator src/Model src/Solver, per the plan.
#
# Exit 0 if every match is either absent or explicitly allowlisted below.
# Exit 1 if any non-allowlisted match is found.
set -uo pipefail
cd "$(git rev-parse --show-toplevel)"

SCAN_DIRS="src/Operator src/Integrator src/Model src/Solver"
VIOLATIONS_FILE="$(mktemp)"
trap 'rm -f "$VIOLATIONS_FILE"' EXIT

# Inline allowlist: exact "file:line" -> reason. Every entry must cite the
# BUG_PATTERNS.md section it corresponds to.
declare -A ALLOWLIST=(
  ["src/Integrator/PhaseFieldMicrostructure.cpp:449"]="dormant landmine: bare 'volume' accumulator inside AMREX_GPU_DEVICE lambda, implicit this-> capture. Inert only because this integrator is not in the GPU-supported closure. Must move to ReduceOps (like Base::Mechanics::Integrate) before GPU-enabling. See docs/llm/BUG_PATTERNS.md #2 Site B."
  ["src/Integrator/PhaseFieldMicrostructure.cpp:459"]="same as :449 (bare 'area' accumulator). BUG_PATTERNS.md #2 Site B."
  ["src/Integrator/PhaseFieldMicrostructure.cpp:463"]="same as :449 (bare 'gbenergy' accumulator). BUG_PATTERNS.md #2 Site B."
  ["src/Integrator/PhaseFieldMicrostructure.cpp:465"]="same as :449 (bare 'realgbenergy' accumulator). BUG_PATTERNS.md #2 Site B."
  ["src/Integrator/PhaseFieldMicrostructure.cpp:466"]="same as :449 (bare 'regenergy' assignment). BUG_PATTERNS.md #2 Site B."
  ["src/Integrator/PhaseFieldMicrostructure.cpp:473"]="same as :449 (bare 'gbenergy' accumulator, anisotropic branch). BUG_PATTERNS.md #2 Site B."
  ["src/Integrator/PhaseFieldMicrostructure.cpp:476"]="same as :449 (bare 'realgbenergy' accumulator, anisotropic branch). BUG_PATTERNS.md #2 Site B."
  ["src/Integrator/PhaseFieldMicrostructure.cpp:481"]="same as :449 (bare 'regenergy' accumulator, anisotropic branch). BUG_PATTERNS.md #2 Site B."
  ["src/Integrator/PhaseFieldMicrostructure.cpp:483"]="same as :449 (bare 'gbenergy' accumulator, 3D branch). BUG_PATTERNS.md #2 Site B."
)

report() {
  local loc="$1" msg="$2"
  local reason="${ALLOWLIST[$loc]:-}"
  if [ -n "$reason" ]; then
    echo "ALLOWLISTED $loc -- $reason"
  else
    echo "VIOLATION   $loc -- $msg"
    echo "$loc" >> "$VIOLATIONS_FILE"
  fi
}

echo "== Rule 1: local FArrayBox temp without a matching .elixir() (UAF) =="
echo "   calibrated against c00f69086 (pre-fix had 0 'elixir' hits in the file)"
while IFS=: read -r file line rest; do
  varname=$(echo "$rest" | sed -E 's/^[[:space:]]*(amrex::)?FArrayBox[[:space:]]+([A-Za-z_][A-Za-z0-9_]*)[[:space:]]*;.*/\2/')
  [ -z "$varname" ] && continue
  if awk -v ln="$line" -v var="$varname" 'NR>ln && $0 ~ (var "\\.elixir\\(\\)") {f=1} END{exit !f}' "$file"; then
    continue
  fi
  report "$file:$line" "local FArrayBox '$varname' has no matching .elixir() call later in the file (BUG_PATTERNS.md #1)"
done < <(grep -rnE '^[[:space:]]*(amrex::)?FArrayBox[[:space:]]+[A-Za-z_][A-Za-z0-9_]*[[:space:]]*;' $SCAN_DIRS)

echo "== Rule 2a: array-typed member hoisted via 'auto x = this->member;' (decays to pointer) =="
echo "   calibrated against c00f69086~1..76ae550e0 (trac_hi/disp_hi in Mechanics.H)"
ARRAY_MEMBERS=$(grep -rn --include='*.H' -E '^\s*[A-Za-z_:]+(<[^>]*>)?\s+[A-Za-z_][A-Za-z0-9_]*\s*\[[A-Za-z0-9_]+\]\s*;' $SCAN_DIRS 2>/dev/null \
  | grep -v 'return ' | grep -oE '[A-Za-z_][A-Za-z0-9_]*\s*\[' | sed -E 's/[[:space:]]*\[//' | sort -u)
for m in $ARRAY_MEMBERS; do
  while IFS=: read -r file line rest; do
    report "$file:$line" "'auto x = this->$m;' hoists an array-typed member; auto deduction decays it to a pointer, so a lambda capturing the local still dereferences host memory on device (BUG_PATTERNS.md #2 Site A)"
  done < <(grep -rnE "auto[[:space:]]+[A-Za-z_][A-Za-z0-9_]*[[:space:]]*=[[:space:]]*this->${m}[[:space:]]*;" $SCAN_DIRS 2>/dev/null)
done

echo "== Rule 2b: bare top-level member mutated inside an AMREX_GPU_DEVICE lambda (implicit this capture) =="
echo "   calibrated against PhaseFieldMicrostructure.cpp (dormant, allowlisted below) -- no fixed"
echo "   before/after pair exists for this exact site, so it is enforced as a live finding, not"
echo "   silently ignored; see docs/llm/BUG_PATTERNS.md #2 Site B."
for cpp in $(grep -rlE 'AMREX_GPU_DEVICE' $SCAN_DIRS --include='*.cpp' 2>/dev/null); do
  hdr="${cpp%.cpp}.H"
  [ -f "$hdr" ] || continue
  # Only class-scope members (4-space indent) to avoid matching fields of
  # nested config structs (e.g. `pf.L`), which are always accessed qualified.
  members=$(grep -oE '^    (Set::Scalar|Set::Vector|Set::Matrix|amrex::Real|double|int)\s+[A-Za-z_][A-Za-z0-9_]*\s*(=|;)' "$hdr" 2>/dev/null \
    | grep -oE '[A-Za-z_][A-Za-z0-9_]*\s*(=|;)$' | sed -E 's/[[:space:]]*[=;]$//' | sort -u)
  for m in $members; do
    while IFS=: read -r file line content; do
      report "$file:$line" "bare member '$m' mutated inside an AMREX_GPU_DEVICE lambda without this-> (implicit this capture; BUG_PATTERNS.md #2 Site B)"
    done < <(awk -v member="$m" -v file="$cpp" '
      /AMREX_GPU_DEVICE/ {inlambda=1}
      inlambda && $0 ~ ("(^|[^A-Za-z0-9_.])" member "[[:space:]]*(\\+=|=[^=])") \
               && $0 !~ ("this->" member) \
               && $0 !~ ("(Set::Scalar|Set::Vector|Set::Matrix|amrex::Real|double|int|auto)[[:space:]]+" member "[[:space:]]*=") {
        print file ":" NR ":" $0
      }
      inlambda && /\}\);/ {inlambda=0}
    ' "$cpp")
  done
done

echo "== Rule 3: chained Eigen expression templates (e.g. .inverse().transpose()) =="
echo "   calibrated against a5c1b2ddf (NeoHookean.H DW/DDW, CUDA error 719)"
while IFS=: read -r file line content; do
  # skip pure-comment lines (leading // after trimming whitespace)
  trimmed="${content#"${content%%[![:space:]]*}"}"
  case "$trimmed" in
    "//"*) continue ;;
  esac
  report "$file:$line" "chained Eigen expression must be materialized into a named intermediate before the next operation composes on it (BUG_PATTERNS.md #3)"
done < <(grep -rnE '\.(inverse|cross|normalized)\(\)\.(transpose|normalized|cross|inverse)\(' $SCAN_DIRS)

echo
n_violations=$(wc -l < "$VIOLATIONS_FILE" | tr -d ' ')
echo "non-allowlisted violations: $n_violations"
[ "$n_violations" -eq 0 ]
