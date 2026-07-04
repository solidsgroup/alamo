#!/usr/bin/env python3
"""Parse `strace -c` summaries (CPU & GPU) into a classified I/O report.

strace -c table columns: % time | seconds | usecs/call | calls | errors | syscall

We bucket syscalls so the report answers "how much wall time is locked in I/O?"
broken down by write / read / sync / metadata / gpu-driver / memory / other.

Usage: parse_strace.py <cpu_summary> <gpu_summary> <results_dir>
"""
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
import svgchart  # noqa: E402

BUCKETS = {
    "file_write": {"write", "pwrite64", "writev", "pwritev"},
    "file_read":  {"read", "pread64", "readv", "preadv"},
    "sync":       {"fsync", "fdatasync", "sync", "sync_file_range", "msync"},
    "metadata":   {"openat", "open", "close", "lseek", "stat", "fstat",
                   "newfstatat", "lstat", "statx", "mkdir", "mkdirat",
                   "rename", "renameat", "renameat2", "unlink", "unlinkat",
                   "access", "faccessat", "getdents64", "fcntl", "ftruncate",
                   "readlink", "readlinkat", "dup", "dup2", "dup3", "pipe2"},
    "gpu_driver": {"ioctl"},
    "memory":     {"mmap", "munmap", "mremap", "brk", "mprotect", "madvise"},
}
IO_BUCKETS = ("file_write", "file_read", "sync", "metadata")

ROW = re.compile(
    r"^\s*([\d.]+)\s+([\d.]+)\s+(\d+)?\s+(\d+)\s+(?:(\d+)\s+)?(\w+)\s*$")


def parse_summary(path):
    if not path or not os.path.exists(path):
        return None
    rows = {}
    total = {"seconds": 0.0, "calls": 0}
    with open(path, errors="replace") as fh:
        for line in fh:
            m = ROW.match(line)
            if not m:
                continue
            secs = float(m.group(2))
            calls = int(m.group(4))
            errors = int(m.group(5)) if m.group(5) else 0
            name = m.group(6)
            if name == "total":
                total = {"seconds": secs, "calls": calls}
                continue
            rows[name] = {"seconds": secs, "calls": calls, "errors": errors}

    bucketed = {b: {"seconds": 0.0, "calls": 0} for b in BUCKETS}
    bucketed["other"] = {"seconds": 0.0, "calls": 0}
    rev = {sc: b for b, scs in BUCKETS.items() for sc in scs}
    for name, r in rows.items():
        b = rev.get(name, "other")
        bucketed[b]["seconds"] += r["seconds"]
        bucketed[b]["calls"] += r["calls"]

    io_secs = sum(bucketed[b]["seconds"] for b in IO_BUCKETS)
    tot = total["seconds"] or sum(r["seconds"] for r in rows.values()) or 1e-12
    return {
        "total_syscall_s": round(tot, 4),
        "total_calls": total["calls"],
        "io_s": round(io_secs, 4),
        "io_pct_of_syscalls": round(100 * io_secs / tot, 2),
        "gpu_driver_s": round(bucketed["gpu_driver"]["seconds"], 4),
        "buckets": {b: {"seconds": round(v["seconds"], 4), "calls": v["calls"]}
                    for b, v in bucketed.items()},
        "top_syscalls": sorted(
            ({"syscall": k, **v} for k, v in rows.items()),
            key=lambda d: d["seconds"], reverse=True)[:12],
    }


def main():
    cpu_p, gpu_p, results = sys.argv[1], sys.argv[2], sys.argv[3]
    cpu, gpu = parse_summary(cpu_p), parse_summary(gpu_p)
    out = {"cpu": cpu, "gpu": gpu}
    with open(os.path.join(results, "io_profile.json"), "w") as fh:
        json.dump(out, fh, indent=2)

    def cell(d, *keys):
        for k in keys:
            d = (d or {}).get(k) if isinstance(d, dict) else None
        return d

    md = ["# I/O Profile (strace syscall accounting)", "",
          "Short instrumented run, plotting ENABLED. `seconds` = wall time spent",
          "inside the syscall (summed across threads).", "",
          "| Metric | CPU | GPU |", "|---|---|---|"]

    def r(label, ckeys, gkeys, fmt="{}"):
        cv = cell(cpu, *ckeys) if cpu else None
        gv = cell(gpu, *gkeys) if gpu else None
        cs = fmt.format(cv) if cv is not None else "-"
        gs = fmt.format(gv) if gv is not None else "-"
        return f"| {label} | {cs} | {gs} |"

    md += [
        r("Total in-syscall time (s)", ["total_syscall_s"], ["total_syscall_s"], "{:.3f}"),
        r("I/O time (write+read+sync+meta) (s)", ["io_s"], ["io_s"], "{:.3f}"),
        r("I/O %% of syscall time", ["io_pct_of_syscalls"], ["io_pct_of_syscalls"], "{:.1f}"),
        r("GPU driver ioctl time (s)", ["gpu_driver_s"], ["gpu_driver_s"], "{:.3f}"),
        r("Total syscalls", ["total_calls"], ["total_calls"], "{}"),
        "",
        "## Time by bucket (s)", "",
        "| Bucket | CPU | GPU |", "|---|---|---|",
    ]
    for b in ("file_write", "file_read", "sync", "metadata", "gpu_driver", "memory", "other"):
        md.append(r(b, ["buckets", b, "seconds"], ["buckets", b, "seconds"], "{:.4f}"))
    md.append("")

    for tag, d in (("CPU", cpu), ("GPU", gpu)):
        if not d:
            continue
        md += [f"## Top syscalls by time -- {tag}", "",
               "| syscall | seconds | calls | errors |", "|---|---|---|---|"]
        for s in d["top_syscalls"]:
            md.append(f"| {s['syscall']} | {s['seconds']:.4f} | {s['calls']} | {s['errors']} |")
        md.append("")

    with open(os.path.join(results, "io_profile.md"), "w") as fh:
        fh.write("\n".join(md))

    # chart: I/O bucket time CPU vs GPU
    cats = ["write", "read", "sync", "metadata", "gpu_ioctl"]
    keys = ["file_write", "file_read", "sync", "metadata", "gpu_driver"]
    def vals(d):
        return [cell(d, "buckets", k, "seconds") for k in keys] if d else [None]*5
    if cpu or gpu:
        svgchart.grouped_bars(
            os.path.join(results, "chart_io.svg"),
            "Time inside syscalls by bucket",
            cats, {"CPU": vals(cpu), "GPU": vals(gpu)}, ylabel="seconds")
    print("\n".join(md))


if __name__ == "__main__":
    main()
