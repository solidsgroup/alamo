"""Verify cp-only restoration, binder-only fitting, and completed simulations."""
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
sys.path.insert(0, str(STUDY / "analysis"))
from compare_pure_htpb import closure, digest
from report_calibration import audit as base_audit

FIXED = HERE / "fixed_properties.json"
PREVIOUS = STUDY / "ap_endpoint_calibration/parameters_frozen.json"


def audit(rows, parameters, relaxations=6, previous=False):
    assert digest(PREVIOUS) == digest(HERE / "parameters_baseline.json")
    baseline = json.loads(PREVIOUS.read_text())
    assert digest(STUDY.parent / "bin/lowmach-2d-clang++") == json.loads(
        (STUDY / "ap_endpoint_calibration/launch_manifest.json").read_text())["binary_sha256"]
    if previous:
        assert parameters == baseline
        fixed = PREVIOUS
    else:
        fixed = FIXED
        original = json.loads((STUDY / "reference/pure_htpb.json").read_text())
        assert parameters["binder"]["cp_J_kg_K"] == original["binder"]["cp_J_kg_K"] == 2418.29
        for key in ("density_kg_m3", "conductivity_W_m_K", "heat_release_J_kg"):
            assert parameters["binder"][key] == baseline["binder"][key]
        assert parameters["binder"]["heat_release_J_kg"] == -66 * 4184
        assert parameters["AP"] == baseline["AP"]
        assert parameters["ap_arrhenius_calibration_history"] == baseline["ap_arrhenius_calibration_history"]
        assert parameters["prior_binder_arrhenius_calibration_history"] == baseline["arrhenius_calibration_history"]
        for stage in parameters["arrhenius_calibration_history"]:
            assert min(stage["balance_temperature_K"]) > parameters["T0_K"]
            for measurement in stage["measurements"]:
                case = json.loads((STUDY / "runs" / measurement["name"] / "case.json").read_text())
                assert case["reference"]["ap_volume_fraction"] == 0., "Mixture entered binder fit"
    for row in rows:
        assert row["surface_temperature_K"] > parameters["T0_K"]
        assert closure(parameters, row["ap_volume_fraction"], row["q_cal_cm2_s"],
                       volume_density=True)[1] > parameters["T0_K"]
    return base_audit(rows, parameters, arrhenius=True, fixed_parameters=fixed,
                      relaxations=relaxations, volume_density=True)
