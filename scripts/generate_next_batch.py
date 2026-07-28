#!/usr/bin/env python3
"""Fit a GP (same dimensions/acquisition as optimize_htpb_fullfeedback_combined.py's
stage 2) to a set of already-completed grid_results/task_*.json evaluations, then ask
the optimizer for a batch of new candidate (pre_exponential, activation_temperature)
points -- the GP surrogate's predicted best locations to evaluate next, used to
replace/refine a Nova grid sweep instead of re-scanning the whole original grid.

Usage:
  python3 scripts/generate_next_batch.py --results-dir <dir with task_*.json> \
      --n-points 35 --out points.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from skopt import Optimizer
from skopt.space import Real

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))


def load_results(results_dir: Path) -> tuple[list[list[float]], list[float]]:
    x, y = [], []
    for path in sorted(results_dir.glob("task_*.json")):
        with open(path) as fh:
            r = json.load(fh)
        x.append([float(np.log10(r["pre_exponential"])), float(r["activation_temperature"])])
        y.append(float(r["residual_norm"]))
    return x, y


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--n-points", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--random-state", type=int, default=0)
    args = parser.parse_args()

    x0, y0 = load_results(args.results_dir)
    print(f"loaded {len(x0)} prior evaluations from {args.results_dir}")

    dimensions = [Real(-6.0, 1.0, name="log_pre_exponential"),
                  Real(0.0, 15000.0, name="activation_temperature")]

    optimizer = Optimizer(dimensions, base_estimator="GP", acq_func="LCB",
                           acq_func_kwargs={"kappa": 0.05},
                           n_initial_points=0, random_state=args.random_state)
    optimizer.tell(x0, y0)

    points = optimizer.ask(n_points=args.n_points, strategy="cl_min")
    out = [{"pre_exponential": float(10 ** p[0]), "activation_temperature": float(p[1])}
           for p in points]

    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)

    print(f"wrote {len(out)} candidate points -> {args.out}")
    for p in out:
        print(f"  pre_exponential={p['pre_exponential']:.6g} "
              f"activation_temperature={p['activation_temperature']:.6g}")


if __name__ == "__main__":
    main()
