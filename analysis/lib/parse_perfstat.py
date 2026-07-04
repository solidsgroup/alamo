#!/usr/bin/env python3
"""Parse `perf stat` text output into derived microarchitecture metrics.

Computes IPC, frontend/backend stall %, cache-miss rate, branch-mispred rate
and host CPU utilization -- the headline numbers an architect reads first.

Usage: parse_perfstat.py <cpu.txt> <gpu.txt> <results_dir>
"""
import json
import os
import re
import sys

NUM = r"[\d,\.]+"


def _f(s):
    return float(s.replace(",", "")) if s and s not in ("<not", "<not>") else None


def parse(path):
    if not path or not os.path.exists(path):
        return None
    txt = open(path, errors="replace").read()
    vals = {}
    for ev in ("task-clock", "context-switches", "cpu-migrations", "page-faults",
               "cycles", "instructions", "branches", "branch-misses",
               "cache-references", "cache-misses",
               "stalled-cycles-frontend", "stalled-cycles-backend"):
        m = re.search(rf"^\s*({NUM})\s+{re.escape(ev)}\b", txt, re.M)
        vals[ev] = _f(m.group(1)) if m else None
    m = re.search(rf"({NUM})\s+seconds time elapsed", txt)
    vals["elapsed_s"] = _f(m.group(1)) if m else None

    cyc, ins = vals.get("cycles"), vals.get("instructions")
    d = {}
    d["ipc"] = round(ins / cyc, 3) if (cyc and ins) else None
    d["cache_miss_pct"] = (round(100 * vals["cache-misses"] / vals["cache-references"], 2)
                           if vals.get("cache-misses") and vals.get("cache-references") else None)
    d["branch_mispred_pct"] = (round(100 * vals["branch-misses"] / vals["branches"], 2)
                               if vals.get("branch-misses") and vals.get("branches") else None)
    d["frontend_stall_pct"] = (round(100 * vals["stalled-cycles-frontend"] / cyc, 2)
                               if vals.get("stalled-cycles-frontend") and cyc else None)
    d["backend_stall_pct"] = (round(100 * vals["stalled-cycles-backend"] / cyc, 2)
                              if vals.get("stalled-cycles-backend") and cyc else None)
    if vals.get("task-clock") and vals.get("elapsed_s"):
        d["host_cpus_utilized"] = round(vals["task-clock"] / 1000.0 / vals["elapsed_s"], 3)
    d["raw"] = vals
    return d


def main():
    cpu, gpu, results = parse(sys.argv[1]), parse(sys.argv[2]), sys.argv[3]
    with open(os.path.join(results, "perfstat.json"), "w") as fh:
        json.dump({"cpu": cpu, "gpu": gpu}, fh, indent=2)

    def r(label, key, fmt="{}"):
        cv = (cpu or {}).get(key)
        gv = (gpu or {}).get(key)
        cs = fmt.format(cv) if cv is not None else "-"
        gs = fmt.format(gv) if gv is not None else "-"
        return f"| {label} | {cs} | {gs} |"

    md = ["# Microarchitecture (perf stat)", "",
          "GPU column = HOST-side counters; low IPC / low host-CPU-utilized is",
          "expected & healthy when work is offloaded to the device.", "",
          "| Metric | CPU | GPU(host) |", "|---|---|---|",
          r("IPC (instructions/cycle)", "ipc", "{:.3f}"),
          r("Frontend stall %", "frontend_stall_pct", "{:.1f}"),
          r("Backend stall %", "backend_stall_pct", "{:.1f}"),
          r("Cache-miss %", "cache_miss_pct", "{:.2f}"),
          r("Branch-mispred %", "branch_mispred_pct", "{:.2f}"),
          r("Host CPUs utilized", "host_cpus_utilized", "{:.2f}"),
          ""]
    if cpu is None and gpu is None:
        md.append("_No perf data captured (counters likely blocked by "
                  "perf_event_paranoid). Run `! sudo sysctl "
                  "kernel.perf_event_paranoid=1` and re-run phase 3._")
    with open(os.path.join(results, "perfstat.md"), "w") as fh:
        fh.write("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
