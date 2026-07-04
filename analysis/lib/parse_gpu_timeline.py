#!/usr/bin/env python3
"""Turn an nvidia-smi --query-gpu CSV time series into stats + an SVG timeline.

Columns (nounits): timestamp, util.gpu, util.mem, mem.used(MiB), power(W),
                   sm_clock(MHz), temp(C)

Reports mean/peak SM utilization (the key "is the GPU actually busy?" metric),
mean power, peak memory, and the fraction of samples with util==0 (idle gaps =
host-bound / launch-latency-bound time, the classic small-kernel killer).

Usage: parse_gpu_timeline.py <csv> <results_dir>
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import svgchart  # noqa: E402


def main():
    csv_path, results = sys.argv[1], sys.argv[2]
    rows = []
    if os.path.exists(csv_path):
        with open(csv_path, errors="replace") as fh:
            for r in csv.reader(fh):
                if not r or "timestamp" in r[0]:
                    continue
                try:
                    rows.append({
                        "util_gpu": float(r[1]), "util_mem": float(r[2]),
                        "mem_used": float(r[3]), "power": float(r[4]),
                        "sm_clk": float(r[5]), "temp": float(r[6]),
                    })
                except (ValueError, IndexError):
                    continue

    def col(k):
        return [x[k] for x in rows] if rows else []

    def stats(k):
        v = col(k)
        return {"mean": round(sum(v) / len(v), 2), "max": max(v),
                "min": min(v)} if v else None

    idle_frac = (round(sum(1 for x in rows if x["util_gpu"] == 0) / len(rows), 3)
                 if rows else None)
    out = {
        "samples": len(rows),
        "util_gpu": stats("util_gpu"),
        "util_mem": stats("util_mem"),
        "mem_used_MiB": stats("mem_used"),
        "power_W": stats("power"),
        "sm_clk_MHz": stats("sm_clk"),
        "idle_sample_fraction": idle_frac,
    }
    with open(os.path.join(results, "gpu_timeline.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    md = ["# GPU Device Timeline (nvidia-smi sampling)", ""]
    if not rows:
        md.append("_No samples captured (nvidia-smi unavailable or run too "
                  "short). For kernel-level flame charts install Nsight Systems "
                  "(`nsys`) and re-run phase 5._")
    else:
        md += [f"- Samples: {len(rows)}",
               f"- **Mean SM utilization: {out['util_gpu']['mean']}%** "
               f"(peak {out['util_gpu']['max']}%)",
               f"- Idle-sample fraction (util==0): {idle_frac} "
               f"(high => host/launch-latency bound)",
               f"- Mean power: {out['power_W']['mean']} W "
               f"(peak {out['power_W']['max']} W)",
               f"- Peak memory used: {out['mem_used_MiB']['max']} MiB",
               f"- Mean SM clock: {out['sm_clk_MHz']['mean']} MHz", ""]
        x = list(range(len(rows)))
        svgchart.timeline(
            os.path.join(results, "chart_gpu_timeline.svg"),
            "GPU utilization & power over run",
            [i * 0.1 for i in x],
            {"SM util %": col("util_gpu"),
             "Mem util %": col("util_mem"),
             "Power W": col("power")},
            xlabel="approx seconds", ylabel="value")
    with open(os.path.join(results, "gpu_timeline.md"), "w") as fh:
        fh.write("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
