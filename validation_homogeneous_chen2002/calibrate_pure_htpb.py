#!/usr/bin/env python3
"""Fit only binder Arrhenius A and E/R to Chen t=0 endpoints; q=500 held out.

All thermal properties and all AP parameters remain fixed. A simulation-based
update compensates endpoint solver/reference ratios without fitting mixtures.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", type=Path, default=STUDY/"reference/pure_htpb.json")
    parser.add_argument("--fixed-parameters", type=Path, default=STUDY/"reference/pure_htpb.json",
                        help="Immutable thermal/AP baseline for this calibration")
    parser.add_argument("--measurements", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Preserve existing calibration stage: {args.output}")
    parameters = json.loads(args.parameters.read_text())
    original = json.loads(args.fixed_parameters.read_text())
    assert parameters["AP"] == original["AP"]
    assert parameters["T0_K"] == original["T0_K"], "Cold temperature must remain fixed"
    b = parameters["binder"]
    for key in ("density_kg_m3", "cp_J_kg_K", "conductivity_W_m_K", "heat_release_J_kg"):
        assert b[key] == original["binder"][key], "Thermal properties must remain fixed"
    curves = list(csv.DictReader((STUDY/"analysis/figure4_curves.csv").open()))
    targets = sorted([r for r in curves if r["kind"] == "eq15_solid"
                      and float(r["t"]) == 0 and float(r["q_cal_cm2_s"]) in (200, 1000)],
                     key=lambda r: float(r["q_cal_cm2_s"]))
    assert len(targets) == 2
    q = np.array([float(r["q_cal_cm2_s"]) for r in targets])
    article_rate = np.array([float(r["r_cm_s"])/100 for r in targets])
    ratios = np.ones(2)
    measurements = []
    if args.measurements:
        rows = json.loads(args.measurements.read_text())
        for i, flux in enumerate(q):
            selected = [r for r in rows if r["ap_volume_fraction"] == 0
                        and r["q_cal_cm2_s"] == flux]
            assert len(selected) == 1, "One completed pure-binder case is required per fit flux"
            row = selected[0]
            case = json.loads((STUDY/"runs"/row["name"]/"case.json").read_text())
            assert case["parameters"] == parameters, "Measurements must use the previous stage"
            assert abs(row["late_speed_drift_percent"]) < .2, "Calibration run is not sufficiently steady"
            ratios[i] = row["r_cm_s"]/row["expected_r_eq14_19_cm_s"]
            measurements.append({k: row[k] for k in ("name", "q_cal_cm2_s", "r_cm_s",
                                 "expected_r_eq14_19_cm_s", "late_speed_drift_percent")})
    desired_rate = article_rate/ratios
    temperature = parameters["T0_K"]+(q*41840/(b["density_kg_m3"]*desired_rate)
                                      + b["heat_release_J_kg"])/b["cp_J_kg_K"]
    assert np.all(temperature > 0), "Fixed thermal data require a nonphysical absolute temperature"
    log_A, activation = np.linalg.solve(np.column_stack((np.ones(2), -1/temperature)), np.log(desired_rate))
    A = float(np.exp(log_A))
    assert A > 0 and activation > 0, "Fit requires positive Arrhenius parameters"
    before = {k: b[k] for k in ("A_m_s", "activation_temperature_K")}
    b.update(A_m_s=A, activation_temperature_K=float(activation))
    parameters["sources"].update(
        status="Arrhenius-only pure-binder calibration to Chen Figure 4 model curves; not experimental material measurements. q=500 and all AP-containing rates are excluded from fitting.",
        kinetics="Only binder A and E/R fitted to q=200 and 1000 pure-binder endpoints; AP kinetics unchanged. Generator converts physical A to the local phase-field rate multiplier.",
        density_cp_conductivity=f"All constituent density, heat capacity and conductivity values remain exactly at {args.fixed_parameters}; none are fitted.",
        heat_release=f"Fixed baseline heats: binder {b['heat_release_J_kg']/4184:g} cal/g and AP {parameters['AP']['heat_release_J_kg']/4184:g} cal/g. Heat applied once during phase change with frozen gas chemistry; neither heat is fitted.",
        calibration="Two endpoint Arrhenius fit using temperatures from the fixed thermal balance; simulation/analytic-rate correction. q=500 held out. This conditional effective fit is not a physical-property validation.")
    parameters.setdefault("arrhenius_calibration_history", []).append(dict(
        previous_parameters=str(args.parameters), fixed_parameters=str(args.fixed_parameters),
        previous_binder_kinetics=before,
        measured_source=str(args.measurements) if args.measurements else None,
        measurements=measurements, q_cal_cm2_s=q.tolist(), target_r_cm_s=(100*article_rate).tolist(),
        solver_to_analytic_ratios=ratios.tolist(), balance_temperature_K=temperature.tolist(),
        fitted_A_m_s=A, fitted_activation_temperature_K=float(activation)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(parameters, indent=2)+"\n")
    print(f"{args.output}: binder A={A:.10g} m/s, E/R={activation:.10g} K")
    print(f"Fixed-thermal balance temperatures: {temperature.tolist()} K; q=500 held out.")


if __name__ == "__main__":
    main()
