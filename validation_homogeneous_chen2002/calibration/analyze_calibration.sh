#!/usr/bin/env bash
set -euo pipefail
for run in validation_homogeneous_chen2002/calibration/runs/T300_N40_R1e5 \
           validation_homogeneous_chen2002/calibration/runs/T300_N80_R1e5 \
           validation_homogeneous_chen2002/calibration/runs/T600_N40_R1e5 \
           validation_homogeneous_chen2002/calibration/runs/T600_N80_R1e5; do
  MPLCONFIGDIR=/tmp/mpl-cache /home/esandall/Software/anaconda3/bin/python \
    scripts/lmrf_interface.py "$run"
done
