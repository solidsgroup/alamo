#!/usr/bin/env python3
"""Convert an AMReX TinyProfiler table into Brendan-Gregg folded-stack format.

TinyProfiler emits flat (non-hierarchical) per-region timings, so the produced
flame graph is a flat set of frames -- one bar per BL_PROFILE region, width =
inclusive average time. That is the honest representation of what TinyProfiler
knows; for a true nested CPU/GPU call-stack flame graph use `nsys` (GPU) or
`perf record` + FlameGraph (CPU).

Usage:  tinyprofiler_to_folded.py <tinyprof.txt>  > out.folded
"""
import re
import sys

# Matches:  <name>  <NCalls> <min> <avg> <max> <pct>%
ROW = re.compile(
    r"^(?P<name>.+?)\s+(?P<ncalls>\d+)\s+"
    r"(?P<mn>[0-9.eE+-]+)\s+(?P<avg>[0-9.eE+-]+)\s+(?P<mx>[0-9.eE+-]+)\s+"
    r"(?P<pct>[0-9.]+)%\s*$"
)


def main() -> int:
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        return 2
    seen = {}
    with open(sys.argv[1], encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = ROW.match(line.rstrip("\n"))
            if not m:
                continue
            name = m.group("name").strip()
            if not name or name.lower().startswith("name"):
                continue
            try:
                avg_s = float(m.group("avg"))
            except ValueError:
                continue
            # Keep the larger reading if a region appears in both tables.
            micros = int(round(avg_s * 1e6))
            if micros <= 0:
                micros = 1
            seen[name] = max(seen.get(name, 0), micros)
    for name, micros in sorted(seen.items(), key=lambda kv: -kv[1]):
        # Sanitize: folded format uses ';' as the stack separator.
        frame = name.replace(";", ":")
        print(f"{frame} {micros}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
