#!/usr/bin/env python3
"""Phase 1 -- wall-clock & efficiency comparison of the completed CPU/GPU runs.

Consumes the /usr/bin/time -v footer + log timers from the two production logs,
emits:  results/wallclock.json, results/wallclock.md, results/wallclock.csv,
        results/chart_wallclock.svg, results/chart_elastic.svg

This is the "production" view: full 1.5 s runs, plotting off.  I/O detail comes
from phase 2 (strace) and the GPU timeline from phase 5.

Usage: 01_wallclock.py <cpu_log> <gpu_log> <results_dir>
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
import parse_time      # noqa: E402
import parse_log       # noqa: E402
import svgchart        # noqa: E402


def load(logpath):
    if not os.path.exists(logpath):
        return None
    with open(logpath, errors="replace") as fh:
        text = fh.read()
    return {"time": parse_time.parse(text), "log": parse_log.parse(text)}


def row(label, cpu, gpu, fmt="{}"):
    cs = fmt.format(cpu) if cpu is not None else "-"
    gs = fmt.format(gpu) if gpu is not None else "-"
    return f"| {label} | {cs} | {gs} |"


def main():
    cpu_log, gpu_log, results = sys.argv[1], sys.argv[2], sys.argv[3]
    cpu, gpu = load(cpu_log), load(gpu_log)

    data = {"cpu": cpu, "gpu": gpu}
    with open(os.path.join(results, "wallclock.json"), "w") as fh:
        json.dump(data, fh, indent=2)

    ct = (cpu or {}).get("time", {})
    gt = (gpu or {}).get("time", {})
    cl = (cpu or {}).get("log", {}).get("summary", {})
    gl = (gpu or {}).get("log", {}).get("summary", {})

    cw, gw = ct.get("wall_s"), gt.get("wall_s")
    speedup = (cw / gw) if (cw and gw) else None

    md = ["# Wall-clock & Efficiency Report",
          "",
          f"- CPU binary log: `{os.path.basename(cpu_log)}`",
          f"- GPU binary log: `{os.path.basename(gpu_log)}`",
          ""]
    if speedup:
        verdict = "GPU faster" if speedup > 1 else "CPU faster"
        md.append(f"**Headline: CPU/GPU wall-clock ratio = {speedup:.2f}x "
                  f"({verdict}).**\n")

    md += ["## Resource summary", "",
           "| Metric | CPU | GPU |", "|---|---|---|",
           row("Wall clock (s)", cw, gw, "{:.2f}"),
           row("User CPU (s)", ct.get("user_s"), gt.get("user_s"), "{:.2f}"),
           row("System CPU (s)", ct.get("sys_s"), gt.get("sys_s"), "{:.2f}"),
           row("CPU time = user+sys (s)", ct.get("cpu_time_s"), gt.get("cpu_time_s"), "{:.2f}"),
           row("Avg active cores (cpu_time/wall)", ct.get("avg_active_cores"),
               gt.get("avg_active_cores"), "{:.2f}"),
           row("Non-compute / stall (s)", ct.get("non_compute_s"),
               gt.get("non_compute_s"), "{:.2f}"),
           row("Non-compute / stall (%)", ct.get("non_compute_pct"),
               gt.get("non_compute_pct"), "{:.1f}"),
           row("Peak RSS (MB)", ct.get("max_rss_MB"), gt.get("max_rss_MB"), "{:.1f}"),
           row("FS bytes written (MB)", ct.get("fs_outputs_MB"),
               gt.get("fs_outputs_MB"), "{:.2f}"),
           row("Major (I/O) page faults", ct.get("major_faults"),
               gt.get("major_faults"), "{}"),
           row("Voluntary ctx switches", ct.get("vol_ctxsw"),
               gt.get("vol_ctxsw"), "{}"),
           row("Involuntary ctx switches", ct.get("invol_ctxsw"),
               gt.get("invol_ctxsw"), "{}"),
           ""]

    md += ["## Solver work (from MLMG timers in log)", "",
           "| Metric | CPU | GPU |", "|---|---|---|",
           row("Steps completed", cl.get("steps_completed"), gl.get("steps_completed"), "{}"),
           row("Final sim time (s)", cl.get("final_sim_time"), gl.get("final_sim_time"), "{:.4f}"),
           row("Elastic solves", cl.get("n_elastic_solves"), gl.get("n_elastic_solves"), "{}"),
           row("Elastic total (s)", cl.get("elastic_total_s"), gl.get("elastic_total_s"), "{:.2f}"),
           row("Elastic mean/solve (s)", cl.get("elastic_mean_s"), gl.get("elastic_mean_s"), "{:.3f}"),
           row("Elastic max/solve (s)", cl.get("elastic_max_s"), gl.get("elastic_max_s"), "{:.3f}"),
           row("MLMG iters mean", cl.get("mlmg_iters_mean"), gl.get("mlmg_iters_mean"), "{}"),
           ""]

    # Per-step normalized cost (the fair "throughput" metric) ------------------
    def per_step_ms(t, l):
        w, n = t.get("wall_s"), l.get("steps_completed")
        return round(1000.0 * w / n, 3) if (w and n) else None
    md += ["## Throughput", "",
           "| Metric | CPU | GPU |", "|---|---|---|",
           row("Wall per step (ms)", per_step_ms(ct, cl), per_step_ms(gt, gl), "{:.3f}"),
           ""]

    with open(os.path.join(results, "wallclock.md"), "w") as fh:
        fh.write("\n".join(md))

    # CSV --------------------------------------------------------------------
    with open(os.path.join(results, "wallclock.csv"), "w") as fh:
        fh.write("metric,cpu,gpu\n")
        for k in ("wall_s", "user_s", "sys_s", "cpu_time_s", "avg_active_cores",
                  "non_compute_s", "non_compute_pct", "max_rss_MB",
                  "fs_outputs_MB", "major_faults"):
            fh.write(f"{k},{ct.get(k)},{gt.get(k)}\n")

    # Charts ------------------------------------------------------------------
    if cw or gw:
        svgchart.grouped_bars(
            os.path.join(results, "chart_wallclock.svg"),
            "Wall clock vs CPU-time vs stall (s)",
            ["Wall", "CPU time", "Stall/IO"],
            {"CPU": [ct.get("wall_s"), ct.get("cpu_time_s"), ct.get("non_compute_s")],
             "GPU": [gt.get("wall_s"), gt.get("cpu_time_s"), gt.get("non_compute_s")]},
            ylabel="seconds")
    if cl.get("elastic_total_s") or gl.get("elastic_total_s"):
        svgchart.grouped_bars(
            os.path.join(results, "chart_elastic.svg"),
            "Elastic MLMG solve cost",
            ["total (s)", "mean/solve (s)", "max/solve (s)"],
            {"CPU": [cl.get("elastic_total_s"), cl.get("elastic_mean_s"), cl.get("elastic_max_s")],
             "GPU": [gl.get("elastic_total_s"), gl.get("elastic_mean_s"), gl.get("elastic_max_s")]},
            ylabel="seconds")

    print("\n".join(md))
    if speedup:
        print(f"\n[wallclock] CPU/GPU ratio = {speedup:.2f}x")


if __name__ == "__main__":
    main()
