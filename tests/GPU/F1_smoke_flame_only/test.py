#!/usr/bin/env python3
import os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
import testlib_gpu

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"

binary = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-2d-cuda86-g++")
if binary is None:
    env_bin = os.environ.get("ALAMO_GPU_BIN")
    binary = Path(env_bin) if env_bin else None
if binary is None:
    print("SKIP: no GPU binary found (bin/alamo_gpu-2d-cuda86-g++)")
    sys.exit(77)

rc, log, elapsed = testlib_gpu.run_alamo(
    binary, input_file, outdir,
    extra_args=["max_step=20"],
    timeout=120)

print(log[-2000:])

if rc != 0:
    print(f"FAIL: alamo exited {rc}")
    sys.exit(1)

thermo = testlib_gpu.parse_thermo(outdir / "plot" / "thermo.dat")
ok, issues = testlib_gpu.check_sanity(thermo)
if not ok:
    for i in issues:
        print("FAIL:", i)
    sys.exit(1)

pressure = thermo.get("chamber_pressure", [])
if pressure and pressure[-1] <= 0:
    print(f"FAIL: chamber_pressure not positive ({pressure[-1]})")
    sys.exit(1)

steps = len(thermo.get("time", [])) - 1
print(f"PASS: F1 smoke flame-only ({elapsed:.1f}s, {steps} steps)")
sys.exit(0)
