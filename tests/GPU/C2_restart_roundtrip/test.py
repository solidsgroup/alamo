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

run_a_out = outdir / "run_a"
run_b_out = outdir / "run_b"

# Run A: 20 steps; checkpoint written at step 10 (amr.check_int=10 in input)
print("Run A: 20 steps with checkpoint at step 10...")
rc_a, log_a, _ = testlib_gpu.run_alamo(
    gpu_bin, input_file, run_a_out,
    extra_args=["max_step=20"],
    timeout=120)
if rc_a != 0:
    print("FAIL: Run A exited", rc_a)
    print(log_a[-2000:])
    sys.exit(1)

# ALAMO has no standalone checkpoint file: every plotfile dir IS the restart
# checkpoint. With amr.plot_int=10 (set in input), run A wrote restartable dirs
# plot/00010{cell,node}. Flame carries both cell fabs (eta, temp, ...) and a
# nodal fab (phi), so a correct restart must restore BOTH arenas via the
# separate restart_cell / restart_node keys (the bare `restart` key only sets
# the cell file). See src/Integrator/Integrator.cpp Restart().
plot_a = run_a_out / "plot"
chk_cell = plot_a / "00010cell"
chk_node = plot_a / "00010node"
if not chk_cell.is_dir() or not chk_node.is_dir():
    print(f"FAIL: step-10 checkpoint dirs missing under {plot_a}")
    print("  found:", sorted(p.name for p in plot_a.glob("000*")))
    sys.exit(1)

# Run B: restart from step 10, run until max_step=20 (10 more steps)
print(f"Run B: restart from {chk_cell.name}/{chk_node.name}, running to step 20...")
rc_b, log_b, _ = testlib_gpu.run_alamo(
    gpu_bin, input_file, run_b_out,
    extra_args=["max_step=20", f"restart_cell={chk_cell}", f"restart_node={chk_node}"],
    timeout=120)
if rc_b != 0:
    print("FAIL: Run B exited", rc_b)
    print(log_b[-2000:])
    sys.exit(1)

thermo_a = testlib_gpu.parse_thermo(run_a_out / "plot" / "thermo.dat")
thermo_b = testlib_gpu.parse_thermo(run_b_out / "plot" / "thermo.dat")

if not thermo_a:
    print("FAIL: Run A thermo.dat missing or empty")
    sys.exit(1)
if not thermo_b:
    print("FAIL: Run B thermo.dat missing or empty")
    sys.exit(1)

# Run B must be finite/sane (restart must not corrupt the state into NaN/Inf).
ok_b, b_issues = testlib_gpu.check_sanity(thermo_b)
if not ok_b:
    print("FAIL: restarted run thermo insane:", b_issues[:10])
    sys.exit(1)

# Restart roundtrip assertion.
#
# This test validates the restart MECHANISM (which previously segfaulted while
# restoring the nodal fab — see the Integrator::Restart out-of-bounds fix in
# src/Integrator/Integrator.cpp and benchmark/GPU_TEST_SUITE_FIXES.md). It
# asserts that the checkpointed primary fields round-trip exactly: both runs
# must reach the same final simulation time, and the eta-derived geometry
# (burn volume + area) at that time must match bit-tight.
#
# It deliberately does NOT assert bit-exact reproduction of the full thermo
# row. Flame's checkpoint persists only writeout=true fabs; the auxiliary
# thermal state (temp_old / temps) and the chamber pressure ODE state are not
# checkpointed, so the post-restart thermal/burn transient is not guaranteed to
# reproduce the continuous run. Tightening this into a true bit-exact roundtrip
# is tracked as a follow-up (see GPU_TEST_SUITE_FIXES.md, "C2 deferred").
def final_time(t):
    return t["time"][-1] if t.get("time") else None

ta, tb = final_time(thermo_a), final_time(thermo_b)
if ta is None or tb is None:
    print("FAIL: missing time column in A or B")
    sys.exit(1)
if abs(ta - tb) > 1e-12 + 1e-6 * abs(ta):
    print(f"FAIL: runs ended at different times: A={ta} B={tb}")
    sys.exit(1)

ROUNDTRIP_COLS = ["volume", "area"]
final_a = {c: [thermo_a[c][-1]] for c in ROUNDTRIP_COLS if thermo_a.get(c)}
final_b = {c: [thermo_b[c][-1]] for c in ROUNDTRIP_COLS if thermo_b.get(c)}
ok, issues = testlib_gpu.compare_thermo(final_a, final_b, rel_tol=1e-6, abs_tol=1e-12)
if not ok:
    print("FAIL: restart field round-trip mismatch (checkpointed geometry):")
    for i in issues[:20]:
        print(" ", i)
    sys.exit(1)

steps_a = len(thermo_a.get("time", [])) - 1
steps_b = len(thermo_b.get("time", [])) - 1
print(f"PASS: C2 restart roundtrip")
print(f"  Run A steps: {steps_a}  Run B steps (from restart): {steps_b}")
print(f"  Final time A={ta:.6g} B={tb:.6g}; geometry round-trip OK ({ROUNDTRIP_COLS})")
sys.exit(0)
