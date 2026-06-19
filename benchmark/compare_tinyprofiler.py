#!/usr/bin/env python3
"""Side-by-side CPU vs GPU comparison of AMReX TinyProfiler region timings.

Parses two TinyProfiler tables (inclusive average wall-clock per BL_PROFILE
region) and prints a table sorted by the slower of the two, with the GPU
speedup -- the holdup breakdown across compute / IO / solver / regrid.

Usage:  compare_tinyprofiler.py <cpu.tinyprof.txt> <gpu.tinyprof.txt>
"""
import re
import sys

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
                out[name] = max(out.get(name, 0.0), float(m.group("avg")))
            except ValueError:
                pass
    return out


def main() -> int:
    if len(sys.argv) < 3:
        sys.stderr.write(__doc__)
        return 2
    cpu = parse(sys.argv[1])
    gpu = parse(sys.argv[2])
    names = sorted(set(cpu) | set(gpu),
                   key=lambda n: -max(cpu.get(n, 0.0), gpu.get(n, 0.0)))
    if not names:
        print("(no regions parsed -- were both runs built with --profile?)")
        return 1
    print(f"{'region':<42} {'CPU s':>10} {'GPU s':>10} {'speedup':>9}")
    print("-" * 74)
    for n in names[:25]:
        c = cpu.get(n, 0.0)
        g = gpu.get(n, 0.0)
        sp = f"{c / g:6.2f}x" if g > 0 and c > 0 else "    -"
        short = n if len(n) <= 42 else n[:39] + "..."
        print(f"{short:<42} {c:>10.4f} {g:>10.4f} {sp:>9}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
