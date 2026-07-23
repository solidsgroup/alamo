#!/usr/bin/env python3
"""Measure stored stress agreement with reconstructed conservative face stress."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

DIAGNOSTIC = Path(__file__).parents[1] / "20260722-rod-tube-resolution-phi-study"
sys.path.insert(0, str(DIAGNOSTIC))
import analyze_interface_stress as ais  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("plotdir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    header = ais.parse_header(args.plotdir)
    maxima = {"all_finite_nodes_pa": 0.0, "interior_nodes_pa": 0.0,
              "physical_boundary_nodes_pa": 0.0, "levels": {}}
    for level in range(header["finest"] + 1):
        data, _ = ais.dense_level(args.plotdir, level, header)
        stored = np.moveaxis(data[[ais.C_STRESS_XX, ais.C_STRESS_XY,
                                   ais.C_STRESS_YX, ais.C_STRESS_YY]], 0, -1)
        stored = stored.reshape(stored.shape[:2] + (2, 2))
        recovered = ais.face_recovered_stress(data, *header["dx"][level])
        error = np.max(np.abs(stored - recovered), axis=(-2, -1))
        finite = np.isfinite(error)
        boundary = np.zeros(error.shape, dtype=bool)
        boundary[[0, -1], :] = True
        boundary[:, [0, -1]] = True
        for name, mask in (
            ("all_finite_nodes_pa", finite),
            ("interior_nodes_pa", finite & ~boundary),
            ("physical_boundary_nodes_pa", finite & boundary),
        ):
            if np.any(mask):
                maxima[name] = max(maxima[name], float(np.max(error[mask])))
        if np.any(finite & ~boundary):
            masked = np.where(finite & ~boundary, error, -np.inf)
            j, i = np.unravel_index(np.argmax(masked), masked.shape)
            maxima["levels"][str(level)] = {
                "interior_max_pa": float(error[j, i]),
                "index": [int(i), int(j)],
                "position_m": [
                    header["prob_lo"][0] + i * header["dx"][level][0],
                    header["prob_lo"][1] + j * header["dx"][level][1],
                ],
            }

    text = json.dumps(maxima, indent=2) + "\n"
    print(text, end="")
    if args.output:
        args.output.write_text(text)


if __name__ == "__main__":
    main()
