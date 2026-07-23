#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "usage: $0 <pressure-pa> <ncell> <width-m> <normal-shift-in-cells>" >&2
    exit 2
fi

pressure_pa=$1
ncell=$2
width_m=$3
shift_cells=$4

repo_root=$(cd "$(dirname "$0")/../../.." && pwd)
input="$repo_root/docs/agent_plans/20260722-eta-ripple-rootcause/input_planar"
result_root="$repo_root/docs/agent_plans/20260722-eta-ripple-rootcause/results/radial"
h=$(awk -v n="$ncell" 'BEGIN { printf "%.17g", 0.04 / n }')
radius=$(awk -v h="$h" -v s="$shift_cells" \
    'BEGIN { printf "%.17g", 0.02 + h * s }')
label="p${pressure_pa}_n${ncell}_w${width_m}_s${shift_cells}"
output="$result_root/$label"
log="$result_root/$label.log"

if [[ -e "$output" || -e "$log" ]]; then
    echo "refusing to overwrite existing radial run: $label" >&2
    exit 2
fi
mkdir -p "$result_root"

args=(
    "plot_file=$output"
    "amr.n_cell=$ncell $ncell"
    "amr.max_grid_size=$ncell"
    "geometry.is_periodic=0 0"
    "pf.eta.ic.expression.constant.R=$radius"
    "pf.eta.ic.expression.constant.w=$width_m"
    "pf.eta.ic.expression.region0=0.5+0.5*tanh((sqrt((x-0.04)*(x-0.04)+y*y)-R)/w)"
    "elastic.traction=${pressure_pa}_Pa"
    "pf.eta.bc.constant.type.xlo=neumann"
    "pf.eta.bc.constant.type.xhi=neumann"
    "pf.eta.bc.constant.type.ylo=neumann"
    "pf.eta.bc.constant.type.yhi=neumann"
    "thermal.temp.bc.constant.type.xlo=neumann"
    "thermal.temp.bc.constant.type.xhi=neumann"
    "thermal.temp.bc.constant.type.ylo=neumann"
    "thermal.temp.bc.constant.type.yhi=neumann"
    "elastic.bc.constant.type.xlo=trac trac"
    "elastic.bc.constant.type.xhi=disp trac"
    "elastic.bc.constant.type.ylo=trac disp"
    "elastic.bc.constant.type.yhi=trac trac"
    "elastic.bc.constant.type.xloylo=trac disp"
    "elastic.bc.constant.type.xloyhi=trac trac"
    "elastic.bc.constant.type.xhiylo=disp disp"
    "elastic.bc.constant.type.xhiyhi=disp trac"
)

cd "$repo_root"
mpiexec -n 1 bin/alamo-2d-g++ "$input" allow_unused=1 "${args[@]}" 2>&1 | tee "$log"
