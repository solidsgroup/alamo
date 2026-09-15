#!/usr/bin/env bash
# Run from the repository root.  The production input is read-only here;
# every plotfile is written below this calibration directory.
set -euo pipefail
base=tests/LMRFMonoAP/input
out=validation_homogeneous_chen2002/calibration/runs

common=(
  stop_time=1.5e-4_s amr.plot_dt=1.0e-5_s
  AP_decomposition.phase_change.rate_multiplier=1.0e5
  AP_decomposition.phase_change.activation_temperature=0.0
  AP_decomposition.phase_change.temperature_cutoff=0.0
  chemistry.model.rocfire.A="0.0 0.0 0.0 0.0"
  heat_source.ic.expression.constant.qflux=0.0
  include_conduction=0 advect_temperature=0
)

bin/lowmach-2d-clang++ "$base" plot_file="$out/T300_N40_R1e5" \
  amr.n_cell="2 40" amr.max_grid_size=40 "${common[@]}"
bin/lowmach-2d-clang++ "$base" plot_file="$out/T300_N80_R1e5" \
  amr.n_cell="2 80" amr.max_grid_size=80 "${common[@]}"

# Keep the gas state isobaric at 600 K and suppress the unrelated Rocfire
# reactions; only the phase-field mechanism is being calibrated.
hot=(
  Final.density.ic.expression.constant.Tambient=600.0
  Final.density.ic.expression.constant.Tignite=600.0
  temperature.ic.constant.value=600.0
  component_density.bc.expression.constant.T=600.0
  chemistry.model.rocfire.A="0.0 0.0 0.0 0.0"
)
bin/lowmach-2d-clang++ "$base" plot_file="$out/T600_N40_R1e5" \
  amr.n_cell="2 40" amr.max_grid_size=40 "${common[@]}" "${hot[@]}"
bin/lowmach-2d-clang++ "$base" plot_file="$out/T600_N80_R1e5" \
  amr.n_cell="2 80" amr.max_grid_size=80 "${common[@]}" "${hot[@]}"
