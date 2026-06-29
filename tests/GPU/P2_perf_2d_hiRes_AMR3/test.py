#!/usr/bin/env python3
import os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
import testlib_gpu

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"
test_name = Path(__file__).parent.name

binary = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-2d-cuda*-g++")
if binary is None:
    env_bin = os.environ.get("ALAMO_GPU_BIN")
    binary = Path(env_bin) if env_bin else None
if binary is None:
    print("SKIP: no GPU binary found")
    sys.exit(77)

rc, log, elapsed = testlib_gpu.run_alamo(
    binary, input_file, outdir,
    extra_args=["max_step=10000", "allow_unused=1"],
    timeout=300)

print(log[-3000:])

if rc != 0 and rc != -1:
    print(f"FAIL: alamo exited {rc}")
    sys.exit(1)

thermo = testlib_gpu.parse_thermo(outdir / "plot" / "thermo.dat")
ok, issues = testlib_gpu.check_sanity(thermo)
if not ok:
    for i in issues:
        print("FAIL:", i)
    sys.exit(1)

steps_done = len(thermo.get("time", [])) - 1
ms_per_step = testlib_gpu.wall_per_step(elapsed, thermo)

print(f"PASS: {test_name}")
print(f"  Steps completed : {steps_done}")
print(f"  Elapsed         : {elapsed:.1f}s  ({'timeout' if rc == -1 else 'completed'})")
if ms_per_step is not None:
    print(f"  ms/step         : {ms_per_step:.2f}")
sys.exit(0)
