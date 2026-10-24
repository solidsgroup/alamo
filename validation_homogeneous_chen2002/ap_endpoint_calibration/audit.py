"""Execution and fixed-property checks for the pure-AP follow-up."""
import json
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
sys.path.insert(0, str(STUDY / "analysis"))
from compare_pure_htpb import closure, digest
import numpy as np

THERMAL = ("density_kg_m3", "cp_J_kg_K", "conductivity_W_m_K", "heat_release_J_kg")


def fixed_properties(parameters):
    baseline = json.loads((HERE / "parameters_baseline.json").read_text())
    assert parameters["binder"] == baseline["binder"], "Binder calibration must remain fixed"
    assert parameters["T0_K"] == baseline["T0_K"]
    assert parameters["arrhenius_calibration_history"] == baseline["arrhenius_calibration_history"]
    for key in THERMAL:
        assert parameters["AP"][key] == baseline["AP"][key], key
    for stage in parameters.get("ap_arrhenius_calibration_history", []):
        assert stage["q_cal_cm2_s"] == [200., 1000.], "Held-out flux entered fit"
        assert all(m["q_cal_cm2_s"] in (200., 1000.) and m["ap_volume_fraction"] == 1.
                   for m in stage["measurements"]), "Mixture data entered fit"


def audit(rows, parameters, relaxations=6):
    fixed_properties(parameters)
    binary_hash = digest(STUDY.parent / "bin/lowmach-2d-clang++")
    previous_manifest = json.loads((STUDY / "binder_q66_calibration/density_correction/launch_manifest.json").read_text())
    assert binary_hash == previous_manifest["binary_sha256"], "Solver changed since baseline"
    evidence = []
    for row in rows:
        folder = STUDY / "runs" / row["name"]
        case = json.loads((folder / "case.json").read_text())
        receipt = json.loads((folder / "run.json").read_text())
        assert case["parameters"] == parameters, row["name"]
        for key, expected in dict(density_mixing_rule="volume_additive", width_ratio=.25,
                cells_per_width=8, dt_scale=.2, relaxations=relaxations, nx=32,
                gas_conductivity_W_m_K=100., gas_cp_J_kg_K=1000., pressure_Pa=1.e8,
                advect_temperature=False).items():
            assert case[key] == expected, (row["name"], key)
        assert receipt["returncode"] == 0 and receipt["binary_sha256"] == binary_hash
        assert receipt["input_sha256"] == case["input_sha256"] == digest(folder / "input")
        deck = (folder / "input").read_text()
        assert re.findall(r"^\s*chemistry\.model\.type\s*=\s*(\S+)\s*$", deck, re.M) == ["frozen"]
        assert re.findall(r"^\s*mechanisms\.names\s*=\s*(.*?)\s*$", deck, re.M) == ["binder_regression"]
        assert not re.search(r"^\s*[^#\n]*temperature_cutoff\s*=", deck, re.M)
        log = (folder / "output/out.log").read_text()
        steps = re.findall(r"STEP\s+(\d+) ends\. TIME = ([\d.eE+\-]+) DT = ([\d.eE+\-]+)", log)
        assert steps and "AMReX (26.06) finalized" in log
        step, end, dt = map(float, steps[-1])
        assert end >= case["duration_s"] and row["end_time_s"] >= .999 * case["duration_s"]
        assert abs(end-row["end_time_s"]) <= 1.01*dt
        r, T, _ = closure(parameters, row["ap_volume_fraction"], row["q_cal_cm2_s"], volume_density=True)
        assert np.isclose(r, row["expected_r_eq14_19_cm_s"], rtol=1.e-10)
        assert np.isclose(T, case["reference"]["surface_temperature_K"], rtol=1.e-10)
        evidence.append(dict(name=row["name"], returncode=0, binary_sha256=binary_hash,
            input_sha256=receipt["input_sha256"], parameters_match=True, binder_unchanged=True,
            all_thermal_properties_unchanged=True, chemistry_frozen=True, phase_change_only=True,
            density_mixing_rule="volume_additive", requested_duration_s=case["duration_s"],
            final_plot_time_s=row["end_time_s"], final_log_time_s=end, final_step=int(step),
            last_dt_s=dt, closure_recomputed=True, command=receipt["command"]))
    return evidence
