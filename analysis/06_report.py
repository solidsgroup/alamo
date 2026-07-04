#!/usr/bin/env python3
"""Phase 6 -- assemble every artifact into REPORT.md and a self-contained
index.html (SVGs inlined), with an auto-generated architect's analysis.

Usage: 06_report.py <results_dir>
"""
import glob
import html
import json
import os
import sys


def load(results, name):
    p = os.path.join(results, name)
    if os.path.exists(p):
        try:
            return json.load(open(p))
        except json.JSONDecodeError:
            return None
    return None


def g(d, *ks, default=None):
    for k in ks:
        if not isinstance(d, dict):
            return default
        d = d.get(k)
    return d if d is not None else default


def analysis_bullets(wc, io, ps, gpu):
    """Heuristic, metric-driven architect commentary."""
    out = []
    cw = g(wc, "cpu", "time", "wall_s")
    gw = g(wc, "gpu", "time", "wall_s")
    if cw and gw:
        sp = cw / gw
        if sp >= 1.15:
            out.append(f"GPU delivers a **{sp:.2f}x** wall-clock speedup over the "
                       f"single-rank CPU build on this workload.")
        elif sp <= 0.87:
            out.append(f"GPU is **{1/sp:.2f}x slower** than the CPU build here. "
                       f"This workload is still small enough that per-kernel "
                       f"launch latency and host<->device sync dominate, so the "
                       f"device is starved. Expect GPU to win only at larger "
                       f"n_cell / higher max_level / more substeps.")
        else:
            out.append(f"GPU and CPU are within ~15% ({sp:.2f}x): the problem is "
                       f"too small to amortize offload overhead.")
    cstall = g(wc, "cpu", "time", "non_compute_pct")
    if cstall is not None and cstall > 15:
        out.append(f"CPU run spends **{cstall:.0f}%** of wall time NOT on CPU "
                   f"(I/O / sync / sleeps) -- investigate plotfile cadence and "
                   f"MLMG bottom-solve communication.")
    gidle = g(gpu, "idle_sample_fraction")
    gutil = g(gpu, "util_gpu", "mean")
    if gutil is not None:
        if gutil < 35:
            out.append(f"Mean SM utilization is only **{gutil:.0f}%** "
                       f"(idle samples {gidle}). The GPU is host-bound: fuse "
                       f"kernels, batch AMR levels, and cut per-step "
                       f"device synchronizations.")
        elif gutil > 70:
            out.append(f"Mean SM utilization **{gutil:.0f}%** -- the device is "
                       f"well fed; further gains need kernel-level tuning "
                       f"(occupancy, memory coalescing).")
    iocpu = g(io, "cpu", "io_pct_of_syscalls")
    if iocpu is not None and iocpu > 40:
        out.append(f"I/O is **{iocpu:.0f}%** of in-kernel time during plotting; "
                   f"switch to async/HDF5 plotfiles or raise plot_int in "
                   f"production.")
    ipc = g(ps, "cpu", "ipc")
    be = g(ps, "cpu", "backend_stall_pct")
    if ipc is not None:
        tag = ("memory-bound" if (be or 0) > 40 else
               "compute-bound" if ipc > 1.5 else "mixed")
        out.append(f"CPU IPC {ipc} with backend stalls {be}% => **{tag}**.")
    if not out:
        out.append("Insufficient instrumented data captured; run the missing "
                   "phases (perf may be blocked by perf_event_paranoid).")
    return out


def main():
    results = sys.argv[1]
    wc = load(results, "wallclock.json")
    io = load(results, "io_profile.json")
    ps = load(results, "perfstat.json")
    gpu = load(results, "gpu_timeline.json")
    workload_name = os.environ.get("WORKLOAD_NAME", "input_copy")
    workload_note = os.environ.get(
        "WORKLOAD_NOTE",
        "star geometry, void=20 MPa, elastic static MLMG every 50 steps.",
    )

    # ---- REPORT.md (concatenate the per-phase markdown) --------------------
    parts = ["# ALAMO CPU vs GPU -- Performance Analysis", "",
             f"_Workload: `{workload_name}`, {workload_note}_", "",
             "## Executive summary (auto-generated)", ""]
    for b in analysis_bullets(wc, io, ps, gpu):
        parts.append(f"- {b}")
    parts.append("")
    for fn in ("wallclock.md", "io_profile.md", "perfstat.md", "gpu_timeline.md"):
        p = os.path.join(results, fn)
        if os.path.exists(p):
            parts += ["", "---", "", open(p).read()]
    md = "\n".join(parts)
    open(os.path.join(results, "REPORT.md"), "w").write(md)

    # ---- index.html with inlined SVGs -------------------------------------
    def md_to_html(text):
        # tiny markdown: headings, bold, tables, lists, hr
        lines, htmlout, in_tbl = text.splitlines(), [], False
        for ln in lines:
            if ln.startswith("| "):
                cells = [c.strip() for c in ln.strip("|").split("|")]
                if set("".join(cells)) <= set("-: "):
                    continue
                if not in_tbl:
                    htmlout.append("<table>"); in_tbl = True
                tag = "td"
                htmlout.append("<tr>" + "".join(
                    f"<{tag}>{html.escape(c)}</{tag}>" for c in cells) + "</tr>")
                continue
            if in_tbl:
                htmlout.append("</table>"); in_tbl = False
            ln = html.escape(ln)
            ln = ln.replace("**", "")  # bold markers stripped
            if ln.startswith("# "):
                htmlout.append(f"<h1>{ln[2:]}</h1>")
            elif ln.startswith("## "):
                htmlout.append(f"<h2>{ln[3:]}</h2>")
            elif ln.startswith("- "):
                htmlout.append(f"<li>{ln[2:]}</li>")
            elif ln.strip() == "---":
                htmlout.append("<hr>")
            elif ln.strip():
                htmlout.append(f"<p>{ln}</p>")
        if in_tbl:
            htmlout.append("</table>")
        return "\n".join(htmlout)

    svgs = sorted(glob.glob(os.path.join(results, "*.svg")))
    flames = sorted(glob.glob(os.path.join(results, "flamegraph_*.svg")))
    charts = [s for s in svgs if os.path.basename(s).startswith("chart_")]

    body = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>ALAMO CPU vs GPU performance</title><style>",
        "body{font-family:Segoe UI,Helvetica,Arial,sans-serif;margin:40px;"
        "max-width:1200px;color:#111}",
        "table{border-collapse:collapse;margin:12px 0}td{border:1px solid #ddd;"
        "padding:5px 10px}tr:first-child td{background:#f3f4f6;font-weight:600}",
        "h1{border-bottom:3px solid #2563eb}h2{margin-top:34px;color:#1e3a8a}",
        "img,svg,object{max-width:100%}.card{border:1px solid #e5e7eb;"
        "border-radius:10px;padding:16px;margin:18px 0;box-shadow:0 1px 3px "
        "rgba(0,0,0,.06)}li{margin:4px 0}.flame{overflow:auto}</style></head><body>",
        md_to_html(md),
        "<h2>Charts</h2>",
    ]
    for c in charts:
        body.append(f"<div class='card'>{open(c).read()}</div>")
    if flames:
        body.append("<h2>Flame graphs</h2>")
        for f in flames:
            body.append(f"<div class='card flame'><h3>{os.path.basename(f)}</h3>"
                        f"{open(f).read()}</div>")
    body.append("</body></html>")
    open(os.path.join(results, "index.html"), "w").write("\n".join(body))

    print("[report] wrote REPORT.md and index.html")
    print("[report] open: " + os.path.join(results, "index.html"))


if __name__ == "__main__":
    main()
