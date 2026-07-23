#!/usr/bin/env python3
"""Side-by-side CPU vs GPU comparison of AMReX TinyProfiler region timings.

Parses two TinyProfiler tables (inclusive average wall-clock per BL_PROFILE
region) and prints a table sorted by the slower of the two, with the GPU
speedup -- the holdup breakdown across compute / IO / solver / regrid.

Usage:  compare_tinyprofiler.py <cpu.tinyprof.txt> <gpu.tinyprof.txt>
"""
import re
import sys
import tempfile

ROW = re.compile(
    r"^(?P<name>.+?)\s+(?P<ncalls>\d+)\s+"
    r"(?P<mn>[0-9.eE+-]+)\s+(?P<avg>[0-9.eE+-]+)\s+(?P<mx>[0-9.eE+-]+)\s+"
    r"(?P<pct>[0-9.]+)%\s*$"
)


def parse(path):
    out = {}
    try:
        fh = open(path, encoding="utf-8", errors="replace")
    except OSError:
        return out
    with fh:
        for line in fh:
            m = ROW.match(line.rstrip("\n"))
            if not m:
                continue
            name = m.group("name").strip()
            if not name or name.lower().startswith("name"):
                continue
            try:
                row = {"ncalls": int(m.group("ncalls")), "inclusive_wall_seconds": float(m.group("avg"))}
                if name not in out or row["inclusive_wall_seconds"] > out[name]["inclusive_wall_seconds"]:
                    out[name] = row
            except ValueError:
                pass
    return out


def main() -> int:
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write("Operator::Elastic::Fapply() 7 1.0 2.5 3.0 10.0%\n")
            path = f.name
        got = parse(path)
        return 0 if got["Operator::Elastic::Fapply()"]["ncalls"] == 7 and got["Operator::Elastic::Fapply()"]["inclusive_wall_seconds"] == 2.5 else 1
    if len(sys.argv) < 3:
        sys.stderr.write(__doc__)
        return 2
    cpu = parse(sys.argv[1])
    gpu = parse(sys.argv[2])
    names = sorted(set(cpu) | set(gpu),
                   key=lambda n: -max(cpu.get(n, {"inclusive_wall_seconds": 0.0})["inclusive_wall_seconds"],
                                      gpu.get(n, {"inclusive_wall_seconds": 0.0})["inclusive_wall_seconds"]))
    if not names:
        print("(no regions parsed -- were both runs built with --profile?)")
        return 1
    print(f"{'region':<42} {'CPU calls':>10} {'CPU inclusive wall s':>20} {'GPU calls':>10} {'GPU inclusive wall s':>20} {'speedup':>9}")
    print("-" * 118)
    for n in names[:25]:
        c = cpu.get(n, {"ncalls": 0, "inclusive_wall_seconds": 0.0})
        g = gpu.get(n, {"ncalls": 0, "inclusive_wall_seconds": 0.0})
        sp = f"{c['inclusive_wall_seconds'] / g['inclusive_wall_seconds']:6.2f}x" if g["inclusive_wall_seconds"] > 0 and c["inclusive_wall_seconds"] > 0 else "    -"
        short = n if len(n) <= 42 else n[:39] + "..."
        print(f"{short:<42} {c['ncalls']:>10} {c['inclusive_wall_seconds']:>20.4f} {g['ncalls']:>10} {g['inclusive_wall_seconds']:>20.4f} {sp:>9}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
