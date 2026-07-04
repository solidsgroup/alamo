#!/usr/bin/env python3
"""Parse an alamo run log into structured timing data.

Extracts:
  * per-step records (step index, sim TIME, DT)
  * MLMG elastic-solve timers ("MLMG: Timers: Solve=.. Iter=.. Bottom=..")
  * coarse run summary (steps completed, final sim time, MLMG solve count,
    total/avg/max elastic solve seconds)

These give a code-region breakdown of where wall time goes WITHOUT needing the
binary rebuilt with TINY_PROFILE -- the elastic solve is by far the dominant
non-phase-field cost, and the solver prints its own timers.

Usage: parse_log.py <logfile>  ->  JSON on stdout
"""
import json
import re
import statistics
import sys

RE_STEP_END = re.compile(r"^STEP (\d+) ends\. TIME = ([\d.eE+-]+) DT = ([\d.eE+-]+)")
RE_MLMG = re.compile(r"MLMG: Timers: Solve = ([\d.eE+-]+) Iter = ([\d.eE+-]+) Bottom = ([\d.eE+-]+)")
RE_MLMG_ITER = re.compile(r"MLMG: Final Iter\.\s*(\d+)")
RE_ELASTIC = re.compile(r"ELASTIC SOLVE: t=([\d.eE+-]+) step=(\d+)")


def parse(text: str) -> dict:
    steps, solves, mlmg_iters = [], [], []
    for line in text.splitlines():
        m = RE_STEP_END.match(line)
        if m:
            steps.append({"step": int(m.group(1)),
                          "time": float(m.group(2)),
                          "dt": float(m.group(3))})
            continue
        m = RE_MLMG.search(line)
        if m:
            solves.append({"solve_s": float(m.group(1)),
                           "iter_s": float(m.group(2)),
                           "bottom_s": float(m.group(3))})
            continue
        m = RE_MLMG_ITER.search(line)
        if m:
            mlmg_iters.append(int(m.group(1)))

    solve_times = [s["solve_s"] for s in solves]
    summary = {
        "steps_completed": steps[-1]["step"] if steps else 0,
        "final_sim_time": steps[-1]["time"] if steps else 0.0,
        "n_elastic_solves": len(solves),
        "elastic_total_s": round(sum(solve_times), 4),
        "elastic_mean_s": round(statistics.mean(solve_times), 4) if solve_times else 0.0,
        "elastic_max_s": round(max(solve_times), 4) if solve_times else 0.0,
        "mlmg_iters_mean": round(statistics.mean(mlmg_iters), 2) if mlmg_iters else None,
        "mlmg_iters_max": max(mlmg_iters) if mlmg_iters else None,
    }
    return {"summary": summary, "elastic_solves": solves, "steps": steps}


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: parse_log.py <logfile>")
    try:
        with open(sys.argv[1], errors="replace") as fh:
            text = fh.read()
    except FileNotFoundError:
        print(json.dumps({"error": "logfile not found", "path": sys.argv[1]}))
        return
    print(json.dumps(parse(text), indent=2))


if __name__ == "__main__":
    main()
