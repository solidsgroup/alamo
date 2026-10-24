#!/usr/bin/env python3
"""Fit pure-AP A and E/R using the established fixed-thermal two-flux method."""
import argparse
import csv
import json
from pathlib import Path

from audit import HERE, STUDY, audit, closure, digest, fixed_properties
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parameters", type=Path, required=True)
    p.add_argument("--measurements", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    assert not args.output.exists(), "Preserve previous parameter stages"
    parameters = json.loads(args.parameters.read_text())
    fixed_properties(parameters)
    ap = parameters["AP"]
    curves = list(csv.DictReader((STUDY / "analysis/figure4_curves.csv").open()))
    targets = sorted([r for r in curves if r["kind"] == "eq15_solid"
        and float(r["t"]) == 1. and float(r["q_cal_cm2_s"]) in (200., 1000.)],
        key=lambda r: float(r["q_cal_cm2_s"]))
    assert len(targets) == 2
    fluxes = np.array([float(r["q_cal_cm2_s"]) for r in targets])
    target_rates = np.array([float(r["r_cm_s"])/100 for r in targets])
    ratios = np.ones(2)
    measurements = []
    if args.measurements:
        rows = json.loads(args.measurements.read_text())
        for i, q in enumerate(fluxes):
            selected = [r for r in rows if r["ap_volume_fraction"] == 1.
                        and r["q_cal_cm2_s"] == q and r["relaxations"] == 6]
            assert len(selected) == 1
            row = selected[0]
            audit(selected, parameters)
            assert abs(row["late_speed_drift_percent"]) < .2
            ratios[i] = row["r_cm_s"]/row["expected_r_eq14_19_cm_s"]
            measurements.append({k: row[k] for k in ("name", "q_cal_cm2_s", "ap_volume_fraction",
                "r_cm_s", "expected_r_eq14_19_cm_s", "late_speed_drift_percent")})
    desired = target_rates/ratios
    T = parameters["T0_K"] + (41840*fluxes/(ap["density_kg_m3"]*desired)
                                 + ap["heat_release_J_kg"])/ap["cp_J_kg_K"]
    assert np.all(T > 0)
    lnA, activation = np.linalg.solve(np.column_stack((np.ones(2), -1/T)), np.log(desired))
    A = float(np.exp(lnA))
    assert np.isfinite(A) and A > 0 and activation > 0
    before = {k: ap[k] for k in ("A_m_s", "activation_temperature_K")}
    ap.update(A_m_s=A, activation_temperature_K=float(activation))
    parameters["sources"].update(
        scientific_reference="Chen et al. (2002), Figure 4. Binder calibration retained; new fit uses pure-AP t=1 endpoints.",
        status="AP A/E-only calibration to q=200 and 1000 pure-AP model-curve endpoints; q=500 and all mixtures excluded from the new fit.",
        density_cp_conductivity="All thermal properties unchanged from ap_endpoint_calibration/parameters_baseline.json.",
        mixing="Volume-additive density; mass-weighted cp and Q; volume-weighted ln(A) and E/R; existing Chen conductivity rule.",
        kinetics="Binder kinetics remain frozen. Only AP A and E/R are fitted to pure-AP endpoints.",
        calibration="Fixed thermal-balance temperatures, two-point Arrhenius inversion, and measured solver/analytic-rate correction; q=500 held out.")
    parameters.setdefault("ap_arrhenius_calibration_history", []).append(dict(
        previous_parameters=str(args.parameters), previous_parameters_sha256=digest(args.parameters),
        fixed_parameters=str(HERE / "parameters_baseline.json"), previous_AP_kinetics=before,
        measured_source=str(args.measurements) if args.measurements else None,
        measured_source_sha256=digest(args.measurements) if args.measurements else None,
        measurements=measurements, q_cal_cm2_s=fluxes.tolist(), target_r_cm_s=(100*target_rates).tolist(),
        solver_to_analytic_ratios=ratios.tolist(), balance_temperature_K=T.tolist(),
        fitted_A_m_s=A, fitted_activation_temperature_K=float(activation)))
    fixed_properties(parameters)
    for q, r in zip(fluxes, desired):
        assert np.isclose(closure(parameters, 1., q, volume_density=True)[0], 100*r, rtol=1.e-10)
    args.output.write_text(json.dumps(parameters, indent=2)+"\n")
    print(f"AP A={A:.10g} m/s; E/R={activation:.10g} K; E={activation*8.31446261815324/1000:.8g} kJ/mol")
    print("Thermal-balance temperatures [K]:", T.tolist())
    print("Measured solver/analytic ratios:", ratios.tolist())


if __name__ == "__main__":
    main()
