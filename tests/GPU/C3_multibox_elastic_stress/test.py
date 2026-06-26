#!/usr/bin/env python3
import os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
import testlib_gpu

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"

gpu_bin = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-2d-nofast-cuda86-g++")
if gpu_bin is None:
    env = os.environ.get("ALAMO_GPU_STRICT_BIN")
    gpu_bin = Path(env) if env else None
if gpu_bin is None:
    print("SKIP: no GPU strict binary found")
    sys.exit(77)

NUM_STEPS = 200  # -> ~40 elastic solves at interval=5

print(f"Running C3: {NUM_STEPS} steps, elastic every 5, max_grid_size=32 (multi-box)...")
rc, log, elapsed = testlib_gpu.run_alamo(
    gpu_bin, input_file, outdir,
    extra_args=["max_step={}".format(NUM_STEPS)],
    timeout=600)

print(log[-4000:])

if rc != 0:
    print(f"FAIL: alamo exited {rc}")
    sys.exit(1)

# Scan log for divergence signals
lower_log = log.lower()
diverge_signals = [
    "failed to converge", "diverged", "nan detected", "inf detected",
    "mlmg failed", "elastic solver failed", "abort",
]
for sig in diverge_signals:
    if sig in lower_log:
        print(f"FAIL: divergence signal in log: '{sig}'")
        sys.exit(1)

thermo = testlib_gpu.parse_thermo(outdir / "plot" / "thermo.dat")
ok, issues = testlib_gpu.check_sanity(thermo)
if not ok:
    for i in issues:
        print("FAIL:", i)
    sys.exit(1)

steps_done = len(thermo.get("time", [])) - 1
if steps_done < int(NUM_STEPS * 0.95):
    print(f"FAIL: only {steps_done}/{NUM_STEPS} steps completed")
    sys.exit(1)

# Physical sanity: chamber pressure must stay in 0.05-50 MPa
pressure = thermo.get("chamber_pressure", [])
if pressure:
    min_p, max_p = min(pressure), max(pressure)
    if min_p < 5e4 or max_p > 5e7:
        print(f"FAIL: chamber_pressure out of range [{min_p:.3e}, {max_p:.3e}] Pa")
        sys.exit(1)

elastic_solves = steps_done // 5
print(f"PASS: C3 multi-box elastic stress ({elapsed:.1f}s)")
print(f"  Steps: {steps_done}  Elastic solves: ~{elastic_solves}")
sys.exit(0)
