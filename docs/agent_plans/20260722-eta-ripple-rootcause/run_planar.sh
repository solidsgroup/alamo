#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "usage: $0 <x|y|diag> <pressure-pa> <ncell> <width-m> <shift-in-cells>" >&2
    exit 2
fi

orientation=$1
pressure_pa=$2
ncell=$3
width_m=$4
shift_cells=$5

case "$orientation" in
    x) nx=1; ny=0 ;;
    y) nx=0; ny=1 ;;
    diag) nx=0.7071067811865476; ny=0.7071067811865476 ;;
    *) echo "orientation must be x, y, or diag" >&2; exit 2 ;;
esac

repo_root=$(cd "$(dirname "$0")/../../.." && pwd)
input="$repo_root/docs/agent_plans/20260722-eta-ripple-rootcause/input_planar"
result_root="$repo_root/docs/agent_plans/20260722-eta-ripple-rootcause/results/planar"
h=$(awk -v n="$ncell" 'BEGIN { printf "%.17g", 0.04 / n }')
shift_m=$(awk -v h="$h" -v s="$shift_cells" 'BEGIN { printf "%.17g", h * s }')
if [[ "$orientation" == diag ]]; then
    # Put x+y=0.02 m through the interior so the free xhi/yhi boundaries are
    # wholly in eta~=1. shift_cells remains a signed normal displacement.
    shift_m=$(awk -v h="$h" -v s="$shift_cells" \
        'BEGIN { printf "%.17g", -0.01 * sqrt(2.0) + h * s }')
fi
label="${orientation}_p${pressure_pa}_n${ncell}_w${width_m}_s${shift_cells}"
output="$result_root/$label"
log="$result_root/$label.log"

if [[ -e "$output" || -e "$log" ]]; then
    echo "refusing to overwrite existing planar run: $label" >&2
    exit 2
fi
mkdir -p "$result_root"

args=(
    "plot_file=$output"
    "amr.n_cell=$ncell $ncell"
    "amr.max_grid_size=$ncell"
    "pf.eta.ic.expression.constant.nx=$nx"
    "pf.eta.ic.expression.constant.ny=$ny"
    "pf.eta.ic.expression.constant.s0=$shift_m"
    "pf.eta.ic.expression.constant.w=$width_m"
    "elastic.traction=${pressure_pa}_Pa"
)

if [[ "$orientation" == x ]]; then
    args+=(
        "geometry.is_periodic=0 1"
        "pf.eta.bc.constant.type.xlo=neumann"
        "pf.eta.bc.constant.type.xhi=neumann"
        "pf.eta.bc.constant.type.ylo=periodic"
        "pf.eta.bc.constant.type.yhi=periodic"
        "thermal.temp.bc.constant.type.xlo=neumann"
        "thermal.temp.bc.constant.type.xhi=neumann"
        "thermal.temp.bc.constant.type.ylo=periodic"
        "thermal.temp.bc.constant.type.yhi=periodic"
        "elastic.bc.constant.type.xlo=disp disp"
        "elastic.bc.constant.type.xhi=trac trac"
        "elastic.bc.constant.type.ylo=trac trac"
        "elastic.bc.constant.type.yhi=trac trac"
        "elastic.bc.constant.type.xloylo=disp disp"
        "elastic.bc.constant.type.xloyhi=disp disp"
        "elastic.bc.constant.type.xhiylo=trac trac"
        "elastic.bc.constant.type.xhiyhi=trac trac"
    )
elif [[ "$orientation" == y ]]; then
    args+=(
        "geometry.is_periodic=1 0"
        "pf.eta.bc.constant.type.xlo=periodic"
        "pf.eta.bc.constant.type.xhi=periodic"
        "pf.eta.bc.constant.type.ylo=neumann"
        "pf.eta.bc.constant.type.yhi=neumann"
        "thermal.temp.bc.constant.type.xlo=periodic"
        "thermal.temp.bc.constant.type.xhi=periodic"
        "thermal.temp.bc.constant.type.ylo=neumann"
        "thermal.temp.bc.constant.type.yhi=neumann"
        "elastic.bc.constant.type.xlo=trac trac"
        "elastic.bc.constant.type.xhi=trac trac"
        "elastic.bc.constant.type.ylo=disp disp"
        "elastic.bc.constant.type.yhi=trac trac"
        "elastic.bc.constant.type.xloylo=disp disp"
        "elastic.bc.constant.type.xloyhi=trac trac"
        "elastic.bc.constant.type.xhiylo=disp disp"
        "elastic.bc.constant.type.xhiyhi=trac trac"
    )
else
    args+=(
        "geometry.is_periodic=0 0"
        "pf.eta.bc.constant.type.xlo=neumann"
        "pf.eta.bc.constant.type.xhi=neumann"
        "pf.eta.bc.constant.type.ylo=neumann"
        "pf.eta.bc.constant.type.yhi=neumann"
        "thermal.temp.bc.constant.type.xlo=neumann"
        "thermal.temp.bc.constant.type.xhi=neumann"
        "thermal.temp.bc.constant.type.ylo=neumann"
        "thermal.temp.bc.constant.type.yhi=neumann"
        "elastic.bc.constant.type.xlo=disp disp"
        "elastic.bc.constant.type.xhi=trac trac"
        "elastic.bc.constant.type.ylo=disp disp"
        "elastic.bc.constant.type.yhi=trac trac"
        "elastic.bc.constant.type.xloylo=disp disp"
        "elastic.bc.constant.type.xloyhi=disp disp"
        "elastic.bc.constant.type.xhiylo=disp disp"
        "elastic.bc.constant.type.xhiyhi=trac trac"
    )
fi

cd "$repo_root"
mpiexec -n 1 bin/alamo-2d-g++ "$input" allow_unused=1 "${args[@]}" 2>&1 | tee "$log"
