#!/usr/bin/env python3
import os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
import testlib_gpu

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"

binary = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-2d-nofast-cuda*-g++")
if binary is None:
    env_bin = os.environ.get("ALAMO_GPU_STRICT_BIN")
    binary = Path(env_bin) if env_bin else None
if binary is None:
    print("SKIP: no GPU strict binary found (bin/alamo_gpu-2d-nofast-cuda*-g++)")
    sys.exit(77)

rc, log, elapsed = testlib_gpu.run_alamo(
    binary, input_file, outdir,
    extra_args=["max_step=50"],
    timeout=300)

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

steps = len(thermo.get("time", [])) - 1
if steps < 10:
    print(f"FAIL: only {steps} steps completed (expected >= 10)")
    sys.exit(1)

print(f"PASS: F2 smoke elastic ({elapsed:.1f}s, {steps} steps)")
sys.exit(0)
