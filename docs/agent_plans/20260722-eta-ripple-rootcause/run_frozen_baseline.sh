#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <run-label>" >&2
    exit 2
fi

repo_root=$(cd "$(dirname "$0")/../../.." && pwd)
label=$1
nprocs=${NPROCS:-8}
result_root="$repo_root/docs/agent_plans/20260722-eta-ripple-rootcause/results/runs"
output_dir="$result_root/$label"
log_file="$result_root/$label.log"

if [[ -e "$output_dir" || -e "$log_file" ]]; then
    echo "refusing to overwrite existing run: $label" >&2
    exit 2
fi

mkdir -p "$result_root"
cd "$repo_root"

# Restart both centerings on the already-synchronized final hierarchy.  The
# phase chemical potential is made identically zero, the chamber pressure is
# frozen, the one advance uses a negligible dt, and the forced post-restart
# regrid is deferred.  max_step permits exactly one solve/write cycle.
mpiexec -n "$nprocs" bin/alamo-2d-g++ input_rt1s_ideal \
    "plot_file=$output_dir" \
    "restart_cell=$repo_root/output_ripple_no_postsolve_regrid/05002cell" \
    "restart_node=$repo_root/output_ripple_no_postsolve_regrid/05002node" \
    max_step=5003 \
    stop_time=1.0005_s \
    timestep=1.0e-30_s \
    amr.plot_int=1 \
    amr.plot_dt=-1.0_s \
    amr.base_regrid_int=-1 \
    amr.regrid_int=-1 \
    amr.thermo.plot_int=-1 \
    thermal.on=1 \
    thermal.end_initial_refine_time=2.0_s \
    variable_pressure=0 \
    allow_unused=1 \
    propellant.homogenize.mob_prop=false \
    chamber.pressure=4316386.252462386_Pa \
    pf.lambda=0.0_J/m^2 \
    pf.kappa=0.0_J/m^2 \
    small=0.0 \
    elastic.interval=1 \
    elastic.tstart=0 \
    elastic.frozen_restart_fields=1 \
    elastic.traction=4316386.252462386_Pa \
    elastic.traction_from_chamber=0 \
    elastic.print_residual=1 \
    elastic.max_coarsening_level=0 \
    elastic.solver.bottom_solver=bicgstab \
    elastic.bc.constant.type.xlo="trac trac" \
    elastic.bc.constant.type.xhi="disp trac" \
    elastic.bc.constant.type.ylo="trac disp" \
    elastic.bc.constant.type.yhi="trac trac" \
    elastic.bc.constant.type.xloylo="trac disp" \
    elastic.bc.constant.type.xloyhi="trac trac" \
    elastic.bc.constant.type.xhiylo="disp disp" \
    elastic.bc.constant.type.xhiyhi="disp trac" \
    elastic.output_stress_symmetry.xhi=1 \
    elastic.output_stress_symmetry.ylo=1 \
    elastic.plot_casing_support=1 \
    elastic.casing_support_refinement_criterion=0.2 \
    casing_support.ic.type=expression \
    casing_support.ic.expression.region0="0.5 + 0.5*tanh((a4 - sqrt((x-cx)^2 + (y-cy)^2))/w)" \
    casing_support.ic.expression.constant.a4=0.0877 \
    casing_support.ic.expression.constant.cx=0.0877 \
    casing_support.ic.expression.constant.cy=0.0877 \
    casing_support.ic.expression.constant.w=0.00035 \
    2>&1 | tee "$log_file"
