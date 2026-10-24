#!/usr/bin/env python3
"""Audit pure-binder stages; with --final, report two fit points and one holdout."""
import argparse
import csv
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
sys.path.insert(0, str(STUDY / "analysis"))
from compare_pure_htpb import closure, digest, table
from report_calibration import audit
from extract_runs import write_csv
import matplotlib.pyplot as plt
import numpy as np

R = 8.31446261815324  # J/(mol K)
FIXED = HERE / "fixed_properties.json"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parameters", type=Path, required=True)
    p.add_argument("--summary", type=Path, nargs="+", required=True)
    p.add_argument("--final", action="store_true")
    args = p.parse_args()
    parameters = json.loads(args.parameters.read_text())
    rows = sorted([r for path in args.summary for r in json.loads(path.read_text())],
                  key=lambda r: r["q_cal_cm2_s"])
    expected_fluxes = {200., 500., 1000.} if args.final else {200., 1000.}
    assert len(rows) == len(expected_fluxes)
    assert {(r["q_cal_cm2_s"], r["ap_volume_fraction"]) for r in rows} == {
        (q, 0.) for q in expected_fluxes}
    fixed = json.loads(FIXED.read_text())
    assert fixed["AP"] == json.loads((STUDY / "reference/pure_htpb.json").read_text())["AP"]
    assert {k: fixed["binder"][k] for k in (
        "density_kg_m3", "cp_J_kg_K", "conductivity_W_m_K", "heat_release_J_kg")} == {
        "density_kg_m3": 920, "cp_J_kg_K": 2130,
        "conductivity_W_m_K": .213, "heat_release_J_kg": -66 * 4184}
    evidence = audit(rows, parameters, arrhenius=True, fixed_parameters=FIXED)
    out = HERE if args.final else args.summary[0].parent
    (out / "execution_audit.json").write_text(json.dumps(evidence, indent=2) + "\n")
    for r in rows:
        print(f"q={r['q_cal_cm2_s']:g}: {r['r_cm_s']:.8f} cm/s; "
              f"target error {r['error_figure4_solid_percent']:+.5f}%; "
              f"drift {r['late_speed_drift_percent']:+.4f}%; "
              f"Ts={r['surface_temperature_K']:.2f} K; "
              f"incident flux error {r['incident_flux_error_percent']:+.3f}%")
    if not args.final:
        return
    assert max(abs(r["late_speed_drift_percent"]) for r in rows) < .2
    assert max(abs(r["error_figure4_solid_percent"]) for r in rows
               if r["q_cal_cm2_s"] != 500) < 1
    initial = json.loads((HERE / "parameters_00.json").read_text())
    b, b0 = parameters["binder"], initial["binder"]
    comparison = []
    for r in rows:
        q, target = r["q_cal_cm2_s"], r["figure4_solid_r_cm_s"]
        initial_r, initial_T, _ = closure(initial, 0, q)
        final_r, final_T, _ = closure(parameters, 0, q)
        comparison.append(dict(q_cal_cm2_s=q, role="held out" if q == 500 else "fit",
            figure4_r_cm_s=target, initial_analytic_r_cm_s=initial_r,
            initial_analytic_error_percent=100*(initial_r/target-1),
            initial_analytic_Ts_K=initial_T, final_analytic_r_cm_s=final_r,
            final_analytic_Ts_K=final_T, lowmach_r_cm_s=r["r_cm_s"],
            lowmach_error_percent=r["error_figure4_solid_percent"],
            lowmach_Ts_K=r["surface_temperature_K"]))
    write_csv(HERE / "comparison.csv", comparison)
    write_csv(HERE / "simulation_summary.csv", rows)
    (HERE / "simulation_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    history = []
    for path in sorted((STUDY / "analysis/binder_q66").glob("stage*/simulation_summary.json")):
        stage_rows = json.loads(path.read_text())
        stage_parameters = json.loads((STUDY / "runs" / stage_rows[0]["name"] / "case.json").read_text())["parameters"]
        audit(stage_rows, stage_parameters, arrhenius=True, fixed_parameters=FIXED)
        for r in stage_rows:
            history.append(dict(stage=path.parent.name, A_m_s=stage_parameters["binder"]["A_m_s"],
                activation_temperature_K=stage_parameters["binder"]["activation_temperature_K"],
                q_cal_cm2_s=r["q_cal_cm2_s"], r_cm_s=r["r_cm_s"],
                target_error_percent=r["error_figure4_solid_percent"]))
    write_csv(HERE / "stage_history.csv", history)
    provenance_path = STUDY / "analysis/figure4_provenance.json"
    provenance = json.loads(provenance_path.read_text())
    source = Path(provenance["source"])
    assert digest(source) == provenance["sha256"]
    paths = [args.parameters, *args.summary, FIXED, source, provenance_path,
             STUDY / "analysis/figure4_curves.csv", HERE / "parameters_00.json",
             HERE / "analyze.py", STUDY / "calibrate_pure_htpb.py", STUDY / "run_study.py",
             STUDY / "analysis/extract_runs.py", STUDY / "analysis/report_calibration.py"]
    (HERE / "provenance.json").write_text(json.dumps({str(f): digest(f) for f in paths}, indent=2) + "\n")
    # Sensitivity of the analytic fit to a full printed stroke width at each
    # endpoint. This is a graphical perturbation envelope, not a confidence interval.
    resolution = provenance["conservative_graphical_resolution"]["r_cm_s"]
    target = np.array([r["figure4_r_cm_s"] for r in comparison if r["role"] == "fit"])
    sensitivity = []
    for signs in itertools.product((-1, 1), repeat=2):
        rate = (target + resolution*np.array(signs))/100
        T = fixed["T0_K"] + (41840*np.array([200, 1000])/(b["density_kg_m3"]*rate)
                               + b["heat_release_J_kg"])/b["cp_J_kg_K"]
        lnA, activation = np.linalg.solve(np.column_stack((np.ones(2), -1/T)), np.log(rate))
        sensitivity.append(dict(q200_sign=signs[0], q1000_sign=signs[1],
                                A_m_s=float(np.exp(lnA)), activation_temperature_K=float(activation)))
    write_csv(HERE / "graphical_sensitivity.csv", sensitivity)
    grid = np.linspace(180, 1050, 180)
    fig, ax = plt.subplots(figsize=(6.6, 4.3), layout="constrained")
    ax.plot(grid, [closure(initial, 0, q)[0] for q in grid], color="#0072B2",
            label="Analytic endpoint fit")
    ax.plot(grid, [closure(parameters, 0, q)[0] for q in grid], "--", color="#D55E00",
            label="Solver-corrected parameters: analytic prediction")
    ax.errorbar([r["q_cal_cm2_s"] for r in comparison],
                [r["figure4_r_cm_s"] for r in comparison], yerr=resolution,
                fmt="o", mfc="white", color="black", capsize=3, label="Chen Fig. 4, pure binder")
    for role, marker in (("fit", "x"), ("held out", "s")):
        selected = [r for r in comparison if r["role"] == role]
        ax.scatter([r["q_cal_cm2_s"] for r in selected], [r["lowmach_r_cm_s"] for r in selected],
                   marker=marker, s=45, color="#D55E00", label=f"LowMach: {role}", zorder=5)
    ax.set(xlabel="Prescribed heat flux (cal cm⁻² s⁻¹)", ylabel="Regression speed (cm s⁻¹)",
           title="Pure binder: fixed thermal properties, Q = −66 cal/g")
    ax.legend(frameon=False, fontsize=8)
    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"comparison.{ext}", dpi=220)
    held = next(r for r in comparison if r["role"] == "held out")
    comparison_table = table(["q [cal/(cm² s)]", "Use", "Chen [cm/s]", "Initial analytic [cm/s]", "Final LowMach [cm/s]", "LowMach error [%]", "LowMach Ts [K]"],
        [[f"{r['q_cal_cm2_s']:g}", r["role"], f"{r['figure4_r_cm_s']:.6f}",
          f"{r['initial_analytic_r_cm_s']:.6f}", f"{r['lowmach_r_cm_s']:.6f}",
          f"{r['lowmach_error_percent']:+.4f}", f"{r['lowmach_Ts_K']:.2f}"] for r in comparison])
    history_table = table(["Stage", "A [m/s]", "E/R [K]", "q", "LowMach [cm/s]", "Target error [%]"],
        [[r["stage"], f"{r['A_m_s']:.7f}", f"{r['activation_temperature_K']:.3f}",
          f"{r['q_cal_cm2_s']:g}", f"{r['r_cm_s']:.6f}", f"{r['target_error_percent']:+.4f}"] for r in history])
    report = f'''# Pure-binder Arrhenius calibration with fixed user properties

For the existing LowMach calibration setup, the final fitted values are **A_binder = {b['A_m_s']:.8g} m/s** and **E_binder/R = {b['activation_temperature_K']:.8g} K**, equivalent to **E_binder = {b['activation_temperature_K']*R/1000:.7g} kJ/mol** ({b['activation_temperature_K']*R/4184:.7g} kcal/mol). The physical speed law is `r = A_binder exp[-(E_binder/R)/Ts]`.

The direct analytic fit, before correcting for the finite interface and heating surrogate, gives **A_binder = {b0['A_m_s']:.9g} m/s**, **E_binder/R = {b0['activation_temperature_K']:.9g} K**, and **E_binder = {b0['activation_temperature_K']*R/1000:.8g} kJ/mol**. Use this pair for Chen's ideal steady surface balance; use the final pair to reproduce the executed LowMach setup. Numerical correction is conditional on that setup, and is not a separate material measurement.

## Fixed properties and source

Density is 920 kg/m³ (0.92 g/cm³), conductivity 0.213 W/(m K), specific heat 2130 J/(kg K), and Q = −276144 J/kg (−66 cal/g). These are exactly the user's supplied values; only A and E/R were fitted. Initial/deep-solid temperature remains 300 K as in the existing study. All AP values remain those in `reference/pure_htpb.json`. Pure-binder calibration alone does not validate predictions for mixtures.

Targets are the t=0 endpoints of Figure 4 in Chen et al. (2002), printed p. 2927 / PDF page 5, extracted from the original vector paths in `{source}`. Its SHA-256 matches the existing digitization provenance. [Publisher record](https://doi.org/10.1016/S1540-7489(02)80357-1). These are model/DNS curves, not experimental pure-binder measurements. Fit fluxes are 200 and 1000 cal/(cm² s); 500 is held out throughout fitting.

## Method and results

Chen equations (12)–(14) give `q = rho*r*[cp*(Ts-T0) - Q]`. With q converted using 1 cal/(cm² s) = 41840 W/m², each fitting rate fixes `Ts = T0 + [q/(rho*r) + Q]/cp`. Solve the two equations `ln(r) = ln(A) - (E/R)/Ts`. Initial fitting temperatures are {comparison[0]['initial_analytic_Ts_K']:.4f} and {comparison[-1]['initial_analytic_Ts_K']:.4f} K. Conductivity fixes the thermal profile length `delta = k/(rho*cp*r)`; it does not enter this steady endpoint inversion.

The existing `calibrate_pure_htpb.py` now accepts `--fixed-parameters` so the immutable baseline can be supplied explicitly. Subsequent stages divide each target rate by its measured LowMach/analytic rate ratio and repeat the same two-equation inversion. Both endpoint runs must have late speed drift below 0.2%; final endpoint errors must be below 1%. The held-out case is generated only after the parameters are frozen.

{comparison_table}

The initial analytic held-out error is {held['initial_analytic_error_percent']:+.4f}%; final LowMach held-out error is {held['lowmach_error_percent']:+.4f}%. All rates measured from LowMach are obtained from raw eta=0.5 interface-position slopes over the final 40% of each run, with a separate volume-loss check.

{history_table}

## Numerical interpretation and uncertainty

The inherited setup uses ell/delta=0.25, eight cells per ell, 32 transverse cells, dt scale=0.2, and six thermal relaxation times. Heating is supplied by a fixed gas slab, gas conductivity 100 W/(m K), gas cp=1000 J/(kg K), pressure 100 MPa, frozen chemistry and temperature advection disabled. It is a prescribed-flux condensed-phase surrogate. Final incident-flux error ranges from {min(r['incident_flux_error_percent'] for r in rows):+.3f}% to {max(r['incident_flux_error_percent'] for r in rows):+.3f}%; maximum late speed drift is {max(abs(r['late_speed_drift_percent']) for r in rows):.4f}%. Gas thermal storage ranges from {min(r['gas_storage_percent_input'] for r in rows):.3f}% to {max(r['gas_storage_percent_input'] for r in rows):.3f}% of input. Volume-loss and interface-position rates differ by at most {max(abs(100*(r['volume_r_cm_s']/r['r_cm_s']-1)) for r in rows):.4f}%. These numerical effects are partly absorbed by the final fitted pair; this is not a grid-convergence study.

One printed stroke corresponds to about {resolution:.4f} cm/s. Perturbing both fitting endpoints independently by ±one stroke gives an analytic-fit envelope A={min(s['A_m_s'] for s in sensitivity):.3f}–{max(s['A_m_s'] for s in sensitivity):.3f} m/s and E/R={min(s['activation_temperature_K'] for s in sensitivity):.1f}–{max(s['activation_temperature_K'] for s in sensitivity):.1f} K. This is graphical sensitivity, not a statistical confidence interval; extra stored decimals support reproducibility only.

## Using the parameters

`parameters_frozen.json` contains the final constituent parameters for `run_study.py`; `parameters_00.json` contains the direct analytic fit. `activation_temperature` in LowMach is E/R in kelvin, not E in J/mol. Physical A in m/s must be normalized for the chosen phase field: for the study's lambda=mobility=w1=1, w12=2, kappa=3*ell² configuration, `rate_multiplier = A/(1.5*ell)`. The generator performs that conversion for each input. A is not directly interchangeable with `rate_multiplier`, and changes to the phase-field parameters require their corresponding conversion.

For Chen's existing volume-fraction mixing rule, `ln(A_blend)=(1-t)*ln(A_binder)+t*ln(A_AP)` and `(E/R)_blend=(1-t)*(E/R)_binder+t*(E/R)_AP`. Supply pure constituent values once. No AP/blend data enter this fit.

`comparison.csv`, `comparison.pdf`, `simulation_summary.json`, `stage_history.csv`, `execution_audit.json` and `provenance.json` record the results. Previous rejected studies and production input files are preserved. See `REPRODUCE.md` for commands.
'''
    (HERE / "REPORT.md").write_text(report)
    print(comparison_table)


if __name__ == "__main__":
    main()
