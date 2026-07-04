#!/usr/bin/env python3
"""Extract the GNU `/usr/bin/time -v` footer embedded in an alamo run log.

The production runs were launched as `/usr/bin/time -v mpiexec ... > run.log 2>&1`,
so the verbose resource report is the tail of the log.  We turn it into JSON with
normalized, typed fields plus a couple of derived efficiency metrics.

Usage: parse_time.py <logfile>  ->  JSON on stdout
"""
import json
import re
import sys


def _to_seconds(hms: str) -> float:
    """'h:mm:ss' or 'm:ss.ss' -> seconds."""
    parts = hms.strip().split(":")
    parts = [float(p) for p in parts]
    sec = 0.0
    for p in parts:
        sec = sec * 60 + p
    return sec


FIELDS = {
    "wall_s":          (r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", _to_seconds),
    "user_s":          (r"User time \(seconds\):\s*([\d.]+)", float),
    "sys_s":           (r"System time \(seconds\):\s*([\d.]+)", float),
    "cpu_pct":         (r"Percent of CPU this job got:\s*([\d.]+)%", float),
    "max_rss_kb":      (r"Maximum resident set size \(kbytes\):\s*(\d+)", int),
    "major_faults":    (r"Major \(requiring I/O\) page faults:\s*(\d+)", int),
    "minor_faults":    (r"Minor \(reclaiming a frame\) page faults:\s*(\d+)", int),
    "vol_ctxsw":       (r"Voluntary context switches:\s*(\d+)", int),
    "invol_ctxsw":     (r"Involuntary context switches:\s*(\d+)", int),
    "fs_inputs":       (r"File system inputs:\s*(\d+)", int),
    "fs_outputs":      (r"File system outputs:\s*(\d+)", int),
    "exit_status":     (r"Exit status:\s*(\d+)", int),
}


def parse(text: str) -> dict:
    out = {}
    for key, (pat, cast) in FIELDS.items():
        m = re.search(pat, text)
        out[key] = cast(m.group(1)) if m else None

    # Derived architect metrics ------------------------------------------------
    wall = out.get("wall_s") or 0.0
    user = out.get("user_s") or 0.0
    sysd = out.get("sys_s") or 0.0
    cpu_time = user + sysd
    out["cpu_time_s"] = cpu_time
    # Threads kept busy on average (>1 => threaded/parallel; <1 => stalled/idle).
    out["avg_active_cores"] = round(cpu_time / wall, 3) if wall else None
    # Time the process was NOT burning CPU: blocked on I/O, sync, GPU, sleeps.
    # For a single-thread workload this is a direct lower bound on stall time.
    out["non_compute_s"] = round(wall - cpu_time, 3) if wall else None
    out["non_compute_pct"] = round(100.0 * (wall - cpu_time) / wall, 2) if wall else None
    # FS bytes (time reports 512-byte blocks).
    out["fs_outputs_MB"] = round((out.get("fs_outputs") or 0) * 512 / 1e6, 3)
    out["fs_inputs_MB"] = round((out.get("fs_inputs") or 0) * 512 / 1e6, 3)
    out["max_rss_MB"] = round((out.get("max_rss_kb") or 0) / 1024.0, 2)
    return out


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: parse_time.py <logfile>")
    try:
        with open(sys.argv[1], errors="replace") as fh:
            text = fh.read()
    except FileNotFoundError:
        print(json.dumps({"error": "logfile not found", "path": sys.argv[1]}))
        return
    print(json.dumps(parse(text), indent=2))


if __name__ == "__main__":
    main()
