#!/usr/bin/env python3
import os, sys, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
sys.path.insert(0, str(ROOT / "scripts"))
import testlib_gpu

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"

cpu_bin = testlib_gpu.find_binary(ROOT, "bin/alamo-2d-g++")
if cpu_bin is None:
    env = os.environ.get("ALAMO_CPU_BIN")
    cpu_bin = Path(env) if env else None

gpu_bin = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-2d-nofast-cuda*-g++")
if gpu_bin is None:
    env = os.environ.get("ALAMO_GPU_STRICT_BIN")
    gpu_bin = Path(env) if env else None

if cpu_bin is None:
    print("SKIP: CPU binary not found")
    sys.exit(77)
if gpu_bin is None:
    print("SKIP: GPU strict binary not found")
    sys.exit(77)

try:
    import testlib
    import yt
    yt.set_log_level(50)
except ImportError as e:
    print(f"SKIP: scientific Python stack not available ({e})")
    sys.exit(77)

NUM_STEPS = 50
cpu_out = outdir / "cpu"
gpu_out = outdir / "gpu"

print("Running CPU...")
rc_cpu, log_cpu, _ = testlib_gpu.run_alamo(
    cpu_bin, input_file, cpu_out,
    extra_args=["max_step={}".format(NUM_STEPS)],
    timeout=300)
if rc_cpu != 0:
    print("FAIL: CPU run exited", rc_cpu)
    print(log_cpu[-2000:])
    sys.exit(1)

print("Running GPU strict...")
rc_gpu, log_gpu, _ = testlib_gpu.run_alamo(
    gpu_bin, input_file, gpu_out,
    extra_args=["max_step={}".format(NUM_STEPS)],
    timeout=300)
if rc_gpu != 0:
    print("FAIL: GPU run exited", rc_gpu)
    print(log_gpu[-2000:])
    sys.exit(1)

def find_latest_plot(base: Path):
    # AMReX writes plotfiles as <plot_file>/NNNNcell/
    candidates = sorted(base.glob("plot/*cell/"), key=lambda p: p.name)
    if candidates:
        return candidates[-1]
    # Fall back: plot/ itself might be the plotfile directory
    if (base / "plot" / "Header").exists():
        return base / "plot"
    return None

cpu_plot = find_latest_plot(cpu_out)
gpu_plot = find_latest_plot(gpu_out)

if cpu_plot is None or gpu_plot is None:
    print("FAIL: could not find plot directories")
    print("  CPU plot:", cpu_plot, "  GPU plot:", gpu_plot)
    sys.exit(1)

try:
    # Domain: x in [0, 0.004] m; cross-section at y=0
    cpu_df = testlib.readContours(str(cpu_plot),
                                   start=[0.0, 0.0, 0.0],
                                   end=[0.004, 0.0, 0.0],
                                   vars=["eta", "temp"])
    gpu_df = testlib.readContours(str(gpu_plot),
                                   start=[0.0, 0.0, 0.0],
                                   end=[0.004, 0.0, 0.0],
                                   vars=["eta", "temp"])
except Exception as e:
    print("FAIL: could not read plot data:", e)
    sys.exit(1)

import numpy as np
tolerance = 0.05
all_ok = True

for var in ["eta", "temp"]:
    try:
        cpu_vals = np.array(cpu_df[var])
        gpu_vals = np.array(gpu_df[var])
        cpu_x    = np.array(cpu_df["x"])
        gpu_x    = np.array(gpu_df["x"])
    except KeyError:
        print(f"WARN: variable '{var}' not found in output, skipping")
        continue

    # Sort by x, interpolate GPU onto CPU x-points
    cpu_ord = np.argsort(cpu_x)
    gpu_ord = np.argsort(gpu_x)
    gpu_at_cpu = np.interp(cpu_x[cpu_ord], gpu_x[gpu_ord], gpu_vals[gpu_ord])
    cpu_s = cpu_vals[cpu_ord]

    denom = np.maximum(np.abs(cpu_s) + np.abs(gpu_at_cpu), 1e-30)
    rel_err = float(np.max(np.abs(cpu_s - gpu_at_cpu) / denom))
    status = "OK" if rel_err <= tolerance else "FAIL"
    print(f"  {var}: max rel error = {rel_err:.4f} ({status})")
    if rel_err > tolerance:
        all_ok = False

if not all_ok:
    print("FAIL: C4 AMR correctness — GPU/CPU field mismatch exceeds 5%")
    sys.exit(1)

print("PASS: C4 AMR correctness")
sys.exit(0)
