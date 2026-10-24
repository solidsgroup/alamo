#!/usr/bin/env python3
"""Audit and report pure-binder A/E fitting with restored original cp and Q."""
import argparse
import csv
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
sys.path.insert(0, str(STUDY / "analysis"))
from compare_pure_htpb import closure, digest, table
from report_calibration import audit as base_audit
from extract_runs import summarize, write_csv
import matplotlib.pyplot as plt
import numpy as np

FIXED = HERE / "fixed_properties.json"
PREVIOUS = STUDY / "ap_endpoint_calibration/parameters_frozen.json"
R = 8.31446261815324


def audit(rows, parameters, relaxations=6):
    baseline = json.loads(FIXED.read_text())
    previous = json.loads(PREVIOUS.read_text())
    original = json.loads((STUDY / "reference/pure_htpb.json").read_text())
    assert parameters["AP"] == previous["AP"]
    assert parameters["ap_arrhenius_calibration_history"] == previous["ap_arrhenius_calibration_history"]
    assert parameters["prior_binder_arrhenius_calibration_history"] == previous["arrhenius_calibration_history"]
    assert parameters["baseline_change"]["source_sha256"] == digest(PREVIOUS)
    for k in ("cp_J_kg_K", "heat_release_J_kg"):
        assert parameters["binder"][k] == original["binder"][k]
    for k in ("density_kg_m3", "conductivity_W_m_K"):
        assert parameters["binder"][k] == previous["binder"][k]
    assert digest(STUDY.parent / "bin/lowmach-2d-clang++") == json.loads(
        (STUDY / "ap_endpoint_calibration/launch_manifest.json").read_text())["binary_sha256"]
    assert all(r["ap_volume_fraction"] == 0. for r in rows)
    return base_audit(rows, parameters, arrhenius=True, fixed_parameters=FIXED,
                      relaxations=relaxations, volume_density=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parameters", type=Path, required=True)
    p.add_argument("--summaries", type=Path, nargs="+", required=True)
    p.add_argument("--freeze", action="store_true")
    p.add_argument("--final", action="store_true")
    p.add_argument("--long-summary", type=Path)
    args = p.parse_args()
    assert not (args.freeze and args.final)
    parameters = json.loads(args.parameters.read_text())
    rows = sorted([r for path in args.summaries for r in json.loads(path.read_text())],
                  key=lambda r: r["q_cal_cm2_s"])
    wanted = {200., 500., 1000.} if args.final else {200., 1000.}
    assert len(rows) == len(wanted)
    assert {(r["q_cal_cm2_s"], r["ap_volume_fraction"]) for r in rows} == {(q, 0.) for q in wanted}
    evidence = audit(rows, parameters)
    out = HERE if args.final else args.summaries[0].parent
    (out / "execution_audit.json").write_text(json.dumps(evidence, indent=2)+"\n")
    for r in rows:
        print(f"q={r['q_cal_cm2_s']:g}: {r['r_cm_s']:.8f} cm/s; error {r['error_figure4_solid_percent']:+.5f}%; "
              f"drift {r['late_speed_drift_percent']:+.5f}%; Ts={r['surface_temperature_K']:.3f} K")
    assert max(abs(r["late_speed_drift_percent"]) for r in rows) < .2
    if args.freeze or args.final:
        assert max(abs(r["error_figure4_solid_percent"]) for r in rows if r["q_cal_cm2_s"] != 500.) < 1.
    if args.freeze:
        frozen = HERE / "parameters_frozen.json"
        assert not frozen.exists(), "Preserve previous frozen parameters"
        frozen.write_bytes(args.parameters.read_bytes())
        (HERE / "freeze_record.json").write_text(json.dumps(dict(
            source_parameters=str(args.parameters), source_summary=str(args.summaries[0]),
            parameters_sha256=digest(frozen), summary_sha256=digest(args.summaries[0]),
            heldout_flux_cal_cm2_s=500., execution_audit=evidence), indent=2)+"\n")
        print("Frozen:", frozen)
    if not args.final:
        return
    freeze = json.loads((HERE / "freeze_record.json").read_text())
    assert digest(args.parameters) == freeze["parameters_sha256"]
    assert digest(Path(freeze["source_summary"])) == freeze["summary_sha256"]
    long_rows = json.loads(args.long_summary.read_text())
    assert len(long_rows) == 1 and long_rows[0]["q_cal_cm2_s"] == 500.
    long_evidence = audit(long_rows, parameters, relaxations=12)
    (HERE / "duration_execution_audit.json").write_text(json.dumps(long_evidence, indent=2)+"\n")
    original = json.loads((STUDY / "reference/pure_htpb.json").read_text())
    previous = json.loads(PREVIOUS.read_text())
    direct = json.loads((HERE / "parameters_00.json").read_text())
    units = []
    for name, m in (("original binder", original["binder"]), ("previous binder Q=-66", previous["binder"]),
                    ("restored cp/Q: direct analytic binder fit", direct["binder"]),
                    ("restored cp/Q: final binder fit", parameters["binder"]),
                    ("original AP", original["AP"]), ("current AP: unchanged", parameters["AP"])):
        units.append(dict(stage=name, **m, heat_release_cal_g=m["heat_release_J_kg"]/4184,
            E_kJ_mol=m["activation_temperature_K"]*R/1000,
            E_kcal_mol=m["activation_temperature_K"]*R/4184))
    write_csv(HERE / "parameter_comparison.csv", units)
    comparisons = []
    for r in rows:
        analytic, T, _ = closure(parameters, 0., r["q_cal_cm2_s"], volume_density=True)
        initial, Ti, _ = closure(direct, 0., r["q_cal_cm2_s"], volume_density=True)
        comparisons.append(dict(q_cal_cm2_s=r["q_cal_cm2_s"], role="held out" if r["q_cal_cm2_s"] == 500. else "fit",
            chen_r_cm_s=r["figure4_solid_r_cm_s"], lowmach_r_cm_s=r["r_cm_s"],
            error_percent=r["error_figure4_solid_percent"], analytic_r_cm_s=analytic,
            analytic_error_percent=100*(analytic/r["figure4_solid_r_cm_s"]-1),
            direct_analytic_r_cm_s=initial, direct_analytic_Ts_K=Ti,
            analytic_Ts_K=T, lowmach_Ts_K=r["surface_temperature_K"]))
    write_csv(HERE / "comparison.csv", comparisons)
    write_csv(HERE / "simulation_summary.csv", rows)
    (HERE / "simulation_summary.json").write_text(json.dumps(rows, indent=2)+"\n")
    windows = []
    for r in rows+long_rows:
        path = next(path.parent/f"timeseries_{r['name']}.csv" for path in [*args.summaries, args.long_summary]
                    if (path.parent/f"timeseries_{r['name']}.csv").exists())
        raw = [{k: v if k == "plotfile" else float(v) if v else float("nan")
                for k, v in row.items()} for row in csv.DictReader(path.open())]
        for fraction in (.2, .3, .4, .5):
            result = summarize(raw, fraction)
            windows.append(dict(name=r["name"], relaxations=r["relaxations"], late_fraction=fraction,
                r_cm_s=result["r_cm_s"], change_percent=100*(result["r_cm_s"]/r["r_cm_s"]-1)))
    write_csv(HERE / "fit_window_sensitivity.csv", windows)
    held = next(r for r in rows if r["q_cal_cm2_s"] == 500.)
    long = long_rows[0]
    duration = dict(short_duration_s=held["end_time_s"], long_duration_s=long["end_time_s"],
        short_r_cm_s=held["r_cm_s"], long_r_cm_s=long["r_cm_s"],
        rate_change_percent=100*(long["r_cm_s"]/held["r_cm_s"]-1),
        long_error_percent=long["error_figure4_solid_percent"], long_late_drift_percent=long["late_speed_drift_percent"])
    (HERE / "duration_comparison.json").write_text(json.dumps(duration, indent=2)+"\n")
    fig, ax = plt.subplots(figsize=(7, 4.7), layout="constrained")
    grid = np.linspace(180, 1050, 150)
    ax.plot(grid, [closure(direct, 0., q, volume_density=True)[0] for q in grid],
            color=".4", lw=1.2, label="Direct analytic fit")
    ax.plot(grid, [closure(parameters, 0., q, volume_density=True)[0] for q in grid],
            "--", color="#0072B2", label="Solver-calibrated analytic closure")
    ax.scatter([r["q_cal_cm2_s"] for r in comparisons], [r["chen_r_cm_s"] for r in comparisons],
               facecolors="none", edgecolors="black", marker="o", s=65, label="Chen Fig. 4: pure binder")
    for role, marker in (("fit", "s"), ("held out", "D")):
        selected = [r for r in comparisons if r["role"] == role]
        ax.scatter([r["q_cal_cm2_s"] for r in selected], [r["lowmach_r_cm_s"] for r in selected],
                   marker=marker, s=30, color="#0072B2", label=f"LowMach: {role}")
        for r in selected:
            ax.annotate(f"{r['error_percent']:+.2f}%".replace("-", "−"),
                (r["q_cal_cm2_s"], r["lowmach_r_cm_s"]), xytext=(-8, 12), textcoords="offset points",
                ha="right", color="#0072B2", fontsize=9)
    ax.set(xlabel="Prescribed heat flux (cal cm⁻² s⁻¹)", ylabel="Regression speed (cm s⁻¹)",
           title="Pure binder: original specific heat and Q restored")
    ax.legend(frameon=False, fontsize=8)
    for ext in ("png", "pdf"):
        fig.savefig(HERE/f"comparison.{ext}", dpi=220)
    b, ap = parameters["binder"], parameters["AP"]
    E = b["activation_temperature_K"]*R
    low = comparisons[0]
    warm_limit = 100*200*41840/(b["density_kg_m3"]*(-b["heat_release_J_kg"]))
    unit_table = table(["Parameter set", "A [m/s]", "E/R [K]", "E [kJ/mol]", "E [kcal/mol]"],
        [[r["stage"], f"{r['A_m_s']:.8g}", f"{r['activation_temperature_K']:.7f}",
          f"{r['E_kJ_mol']:.7f}", f"{r['E_kcal_mol']:.7f}"] for r in units])
    result_table = table(["q", "Use", "Chen [cm/s]", "LowMach [cm/s]", "Error [%]", "LowMach Ts [K]"],
        [[f"{r['q_cal_cm2_s']:g}", r["role"], f"{r['chen_r_cm_s']:.7f}", f"{r['lowmach_r_cm_s']:.7f}",
          f"{r['error_percent']:+.5f}", f"{r['lowmach_Ts_K']:.3f}"] for r in comparisons])
    report = f'''# Pure-binder recalibration with original cp and Q

The solver-calibrated binder values are **A={b['A_m_s']:.9g} m/s**, **E/R={b['activation_temperature_K']:.9g} K**, **E={E/1000:.8g} kJ/mol={E/4184:.8g} kcal/mol**. Only binder A and E/R were fitted at AP volume fraction zero. AP remains at A={ap['A_m_s']:.9g} m/s and E/R={ap['activation_temperature_K']:.9g} K, or E={ap['activation_temperature_K']*R/1000:.8g} kJ/mol={ap['activation_temperature_K']*R/4184:.8g} kcal/mol.

The original binder specific heat **2418.29 J/(kg K)** and Q=**−300 cal/g=−1255200 J/kg** were restored. Density=920 kg/m³, current conductivity=0.213 W/(m K), T0=300 K and all calibrated AP parameters were retained. Conductivity was not restored to its older 0.13 value because this request restored only cp and Q. Gas chemistry remains frozen and negative Q is applied only as an endothermic phase-change heat. The current volume-additive density and mass-weighted cp/Q rules remain configured.

## Fitted rates and unit reference

{result_table}

![Pure-binder comparison](comparison.png)

{unit_table}

E is converted using R={R:g} J/(mol K) and 1 kcal=4.184 kJ. `activation_temperature` is E/R in kelvin; it is not E in energy units. A is the physical speed prefactor in `r=A*exp[-(E/R)/Ts]`; `run_study.py` converts it to the phase-field rate multiplier for each interface width.

## Interpretation of the restored heat

The direct steady heat balance `q=rho*r*[cp*(Ts-T0)-Q]` at the Chen q=200 target requires **Ts={low['direct_analytic_Ts_K']:.3f} K**, below the 300 K deep-solid temperature. With Ts constrained to at least T0, the supplied heat could support at most **{warm_limit:.6f} cm/s**, below Chen's **{low['chen_r_cm_s']:.6f} cm/s** at that flux. The mathematical fit therefore uses heat from cooling the initially warmer solid in addition to the prescribed gas-side input. The model has its existing zero temperature cutoff, so subambient regression is permitted. The measured q=200 surface temperature is {low['lowmach_Ts_K']:.3f} K.

This explains the much lower fitted activation energy. These are effective parameters for this restored-property surrogate, not independent physical measurements of binder pyrolysis kinetics. No separate surface-temperature constraint was imposed. Comparisons are to the supplied Chen Figure 4 model curves, with graphical resolution approximately 0.005 cm/s. [Chen et al. (2002)](https://doi.org/10.1016/S1540-7489(02)80357-1).

## Calibration and runtime evidence

The established two-flux method fits q=200 and 1000, using temperatures inferred from the fixed thermal balance. A direct analytic fit is first executed; subsequent stages divide each target rate by the previous measured solver/analytic ratio and repeat the two-point Arrhenius inversion. New stages retain their measured-source history. The parameters are frozen only after fresh endpoint errors are below 1% and late-window drift below 0.2%. The q=500 case is held out of every fit and generated after freezing. No mixture or AP run is used in this recalibration, and the earlier full-mixture sweep remains a result for the preceding binder properties.

All three final rates are raw eta=0.5 position slopes over the last 40% of six thermal relaxation times. The numerical setup remains ell/delta=0.25, eight cells per ell, 32 transverse cells, dt scale=0.2, pressure=100 MPa, gas conductivity=100 W/(m K), gas cp=1000 J/(kg K), temperature advection disabled and frozen gas chemistry. Maximum late speed drift is {max(abs(r['late_speed_drift_percent']) for r in rows):.5f}%. Fitting the final 20%, 30%, 40% or 50% changes six-time rates by at most {max(abs(r['change_percent']) for r in windows if r['relaxations']==6):.5f}%. Position and volume-loss rates agree within {max(abs(100*(r['volume_r_cm_s']/r['r_cm_s']-1)) for r in rows):.5f}%.

The held-out q=500 run was also repeated for twelve thermal relaxation times with a deeper cold boundary. Its rate changes from {held['r_cm_s']:.8f} to {long['r_cm_s']:.8f} cm/s ({duration['rate_change_percent']:+.5f}%), with late drift {duration['long_late_drift_percent']:+.5f}%. This is a combined duration/depth check at one flux, and does not replace the standard held-out result. Incident-flux errors for the three standard runs span {min(r['incident_flux_error_percent'] for r in rows):+.4f}% to {max(r['incident_flux_error_percent'] for r in rows):+.4f}%. Rate-fit scatter does not establish total physical or grid-convergence uncertainty.

Cheaper `gpt-5.6-luna` agents launched the simulations; the current/root model performed calibration, raw-output extraction, audits and reporting. All accepted runs have successful observed exits, matching input/executable hashes, finalized logs and full requested duration. The executable and production source were not changed. `parameters_frozen.json`, `parameter_comparison.csv`, `comparison.csv`, execution audits, `duration_comparison.json` and `provenance.json` retain the evidence. See `REPRODUCE.md` for commands.
'''
    (HERE / "REPORT.md").write_text(report)
    provenance_file = STUDY / "analysis/figure4_provenance.json"
    source = json.loads(provenance_file.read_text())
    assert digest(Path(source["source"])) == source["sha256"]
    paths = [FIXED, PREVIOUS, args.parameters, *args.summaries, args.long_summary,
        HERE / "parameters_00.json", HERE / "freeze_record.json", HERE / "REPRODUCE.md", Path(__file__),
        STUDY / "calibrate_pure_htpb.py", STUDY / "run_study.py", STUDY / "analysis/report_calibration.py",
        STUDY / "analysis/compare_pure_htpb.py", STUDY / "analysis/extract_runs.py",
        STUDY / "analysis/figure4_curves.csv", provenance_file, Path(source["source"])]
    (HERE / "provenance.json").write_text(json.dumps({str(path): digest(path) for path in paths}, indent=2)+"\n")
    print(unit_table)
    print(result_table)
    print("Duration change [%]:", duration["rate_change_percent"])


if __name__ == "__main__":
    main()
