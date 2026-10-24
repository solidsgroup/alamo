#!/usr/bin/env python3
"""REJECTED METHOD: historical cp/Q fit, retained only for provenance.

Seed from the sharp-interface balance; an update compensates the measured
solver/reference rate ratios. The fitted values are conditional effective
parameters, not experimental material measurements. q=500 is held out.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

STUDY = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameters", type=Path, default=STUDY/"reference/pure_htpb.json")
    parser.add_argument("--measurements", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Preserve existing calibration stage: {args.output}")
    parameters = json.loads(args.parameters.read_text())
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
            assert len(selected) == 1, "Exactly one completed pure-binder run per fit flux is required"
            row = selected[0]
            assert abs(row["late_speed_drift_percent"]) < .2, "Calibration run is not sufficiently steady"
            ratios[i] = row["r_cm_s"]/row["expected_r_eq14_19_cm_s"]
            measurements.append({k: row[k] for k in ("name", "q_cal_cm2_s", "r_cm_s",
                                 "expected_r_eq14_19_cm_s", "late_speed_drift_percent")})
    desired_rate = article_rate/ratios
    b = parameters["binder"]
    temperature = b["activation_temperature_K"]/np.log(b["A_m_s"]/desired_rate)
    cp, Q = np.linalg.solve(np.column_stack((temperature-parameters["T0_K"], -np.ones(2))),
                            q*41840/(b["density_kg_m3"]*desired_rate))
    assert cp > 0 and Q < 0, "Fit did not retain positive heat capacity and endothermic decomposition"
    before = {k: b[k] for k in ("cp_J_kg_K", "heat_release_J_kg")}
    b.update(cp_J_kg_K=float(cp), heat_release_J_kg=float(Q))
    parameters["sources"].update(
        status="Pure-binder-only calibration to Chen Figure 4 model curves; not experimental measurements. AP and mixed-composition data are excluded from fitting.",
        density_cp_conductivity="Density and conductivity fixed at original provisional constituent values. Only binder cp and decomposition heat are fitted; all AP properties remain fixed.",
        heat_release="Binder heat fitted to Chen t=0 endpoints; Gross -300 cal/g is no longer imposed. AP stays -100 cal/g. Heat is applied once during phase change with frozen chemistry.",
        calibration="q=200 and 1000 cal/(cm2 s) fit binder cp/Q with rho=920, k=0.13, A=10.36, E/R=7500, T0=300 fixed. q=500 is held out. Simulation-based adjustments are resolution/surrogate dependent.")
    parameters.setdefault("calibration_history", []).append(dict(
        previous_parameters=str(args.parameters), previous_binder=before,
        measured_source=str(args.measurements) if args.measurements else None,
        measurements=measurements, q_cal_cm2_s=q.tolist(), target_r_cm_s=(100*article_rate).tolist(),
        solver_to_analytic_ratios=ratios.tolist(), fitted_cp_J_kg_K=float(cp), fitted_Q_J_kg=float(Q)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(parameters, indent=2)+"\n")
    print(f"{args.output}: binder cp={cp:.9g} J/(kg K), Q={Q/4184:.9g} cal/g")
    print("Fit fluxes: 200, 1000; held-out flux: 500. AP parameters unchanged.")


if __name__ == "__main__":
    main()
