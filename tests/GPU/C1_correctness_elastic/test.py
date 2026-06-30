#!/usr/bin/env python3
import os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
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
    print("SKIP: CPU binary not found (bin/alamo-2d-g++)")
    sys.exit(77)
if gpu_bin is None:
    print("SKIP: GPU strict binary not found (bin/alamo_gpu-2d-nofast-cuda*-g++)")
    sys.exit(77)

NUM_STEPS = 30  # 30 steps -> 6 elastic solves at interval=5

cpu_out = outdir / "cpu"
gpu_out = outdir / "gpu"

print("Running CPU reference...")
rc_cpu, log_cpu, elapsed_cpu = testlib_gpu.run_alamo(
    cpu_bin, input_file, cpu_out,
    extra_args=["max_step={}".format(NUM_STEPS)],
    timeout=300)
if rc_cpu != 0:
    print("FAIL: CPU run exited", rc_cpu)
    print(log_cpu[-2000:])
    sys.exit(1)

print("Running GPU strict...")
rc_gpu, log_gpu, elapsed_gpu = testlib_gpu.run_alamo(
    gpu_bin, input_file, gpu_out,
    extra_args=["max_step={}".format(NUM_STEPS)],
    timeout=300)
if rc_gpu != 0:
    print("FAIL: GPU run exited", rc_gpu)
    print(log_gpu[-2000:])
    sys.exit(1)

cpu_thermo = testlib_gpu.parse_thermo(cpu_out / "plot" / "thermo.dat")
gpu_thermo = testlib_gpu.parse_thermo(gpu_out / "plot" / "thermo.dat")

ok_cpu, cpu_issues = testlib_gpu.check_sanity(cpu_thermo)
ok_gpu, gpu_issues = testlib_gpu.check_sanity(gpu_thermo)

if not ok_cpu:
    print("FAIL: CPU thermo insane:", cpu_issues)
    sys.exit(1)
if not ok_gpu:
    print("FAIL: GPU thermo insane:", gpu_issues)
    sys.exit(1)

# Compare deterministic physics and elastic boundary diagnostic columns. The
# traction/displacement diagnostics are included here to guard the
# Base::Mechanics::Integrate() device reduction path that previously triggered
# CUDA error 700 and then had multi-box face-sampling semantics repaired.
PHYSICS_COLS = [
    "time", "volume", "area", "mass_flux", "max_temp",
    "mdot_max", "heatflux_max", "L_max", "eta_min",
    "disp_xhi_x", "disp_xhi_y", "disp_yhi_x", "disp_yhi_y",
    "trac_xhi_x", "trac_xhi_y", "trac_yhi_x", "trac_yhi_y",
]
cpu_phys = {c: cpu_thermo[c] for c in PHYSICS_COLS if c in cpu_thermo}
gpu_phys = {c: gpu_thermo[c] for c in PHYSICS_COLS if c in gpu_thermo}
missing_cpu = [c for c in PHYSICS_COLS if c not in cpu_thermo]
missing_gpu = [c for c in PHYSICS_COLS if c not in gpu_thermo]
if missing_cpu or missing_gpu:
    print("FAIL: missing required thermo columns")
    if missing_cpu:
        print("  CPU missing:", missing_cpu)
    if missing_gpu:
        print("  GPU missing:", missing_gpu)
    sys.exit(1)

ok, diff_issues = testlib_gpu.compare_thermo(cpu_phys, gpu_phys, rel_tol=0.05, abs_tol=1e-10)
if not ok:
    print("FAIL: CPU/GPU thermo mismatch (physics columns):")
    for i in diff_issues[:20]:
        print(" ", i)
    sys.exit(1)

rows = len(cpu_thermo.get("time", []))
print(f"PASS: C1 correctness elastic parity")
print(f"  CPU elapsed: {elapsed_cpu:.1f}s  GPU elapsed: {elapsed_gpu:.1f}s")
print(f"  Physics columns compared: {sorted(cpu_phys)}")
print(f"  Rows compared: {rows}  Elastic solves: ~{NUM_STEPS // 5}")
sys.exit(0)
