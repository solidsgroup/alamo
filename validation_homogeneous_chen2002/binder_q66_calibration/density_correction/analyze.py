#!/usr/bin/env python3
"""Root-model postprocessing of the corrected-density sweep and duration check."""
import argparse
import csv
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BINDER = HERE.parent
STUDY = BINDER.parent
sys.path.insert(0, str(STUDY / "analysis"))
from compare_pure_htpb import closure, digest, table
from report_calibration import audit
from extract_runs import summarize, write_csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    parameters_path = BINDER / "parameters_frozen.json"
    params = json.loads(parameters_path.read_text())
    manifest = json.loads((HERE / "launch_manifest.json").read_text())
    assert digest(parameters_path) == manifest["parameters_sha256"]
    assert digest(STUDY.parent / "bin/lowmach-2d-clang++") == manifest["binary_sha256"]
    all_rows = json.loads(args.summary.read_text())
    rows = sorted([r for r in all_rows if r["relaxations"] == 6],
                  key=lambda r: (r["q_cal_cm2_s"], r["ap_volume_fraction"]))
    long_rows = [r for r in all_rows if r["relaxations"] == 12]
    fluxes, fractions = (200, 500, 1000), (0., .2, .4, .6, .8, 1.)
    assert len(all_rows) == 19 and len(rows) == 18 and len(long_rows) == 1
    assert {(r["q_cal_cm2_s"], r["ap_volume_fraction"]) for r in rows} == {
        (q, t) for q in fluxes for t in fractions}
    assert {(r["q_cal_cm2_s"], r["ap_volume_fraction"]) for r in long_rows} == {(500, .4)}
    evidence = audit(rows, params, arrhenius=True,
                     fixed_parameters=BINDER / "fixed_properties.json", volume_density=True)
    long_evidence = audit(long_rows, params, arrhenius=True,
        fixed_parameters=BINDER / "fixed_properties.json", volume_density=True, relaxations=12)
    (HERE / "execution_audit.json").write_text(json.dumps(evidence, indent=2) + "\n")
    (HERE / "duration_execution_audit.json").write_text(json.dumps(long_evidence, indent=2) + "\n")
    (HERE / "simulation_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    write_csv(HERE / "simulation_summary.csv", rows)
    prior_path = BINDER / "sweep/simulation_summary.json"
    prior = {(r["q_cal_cm2_s"], r["ap_volume_fraction"]): r
             for r in json.loads(prior_path.read_text())}
    assert len(prior) == 18
    old_hash = digest(HERE / "baseline_snapshot/bin/lowmach-2d-clang++")
    for old in prior.values():
        folder = STUDY / "runs" / old["name"]
        receipt = json.loads((folder / "run.json").read_text())
        case = json.loads((folder / "case.json").read_text())
        assert receipt["returncode"] == 0 and receipt["binary_sha256"] == old_hash
        assert receipt["input_sha256"] == case["input_sha256"] == digest(folder / "input")
        assert case["parameters"] == params
        assert np.isclose(closure(params, old["ap_volume_fraction"], old["q_cal_cm2_s"])[0],
                          old["expected_r_eq14_19_cm_s"], rtol=1.e-10)
    comparisons, properties = [], []
    for r in rows:
        t, q = r["ap_volume_fraction"], r["q_cal_cm2_s"]
        case = json.loads((STUDY / "runs" / r["name"] / "case.json").read_text())
        material = case["reference"]
        b, ap = params["binder"], params["AP"]
        rho = (1-t)*b["density_kg_m3"] + t*ap["density_kg_m3"]
        assert np.isclose(material["density_kg_m3"], rho, rtol=1.e-12)
        for key in ("cp_J_kg_K", "heat_release_J_kg"):
            assert np.isclose(rho*material[key],
                (1-t)*b["density_kg_m3"]*b[key]+t*ap["density_kg_m3"]*ap[key], rtol=1.e-12)
        analytic, Ts, _ = closure(params, t, q, volume_density=True)
        old = prior[q, t]
        if t in (0., 1.):
            assert np.isclose(r["r_cm_s"], old["r_cm_s"], rtol=1.e-9), "Pure endpoint changed"
        comparisons.append(dict(q_cal_cm2_s=q, ap_volume_fraction=t,
            chen_r_cm_s=r["figure4_solid_r_cm_s"], old_r_cm_s=old["r_cm_s"],
            corrected_r_cm_s=r["r_cm_s"], old_error_percent=old["error_figure4_solid_percent"],
            corrected_error_percent=r["error_figure4_solid_percent"],
            analytic_r_cm_s=analytic, analytic_error_percent=100*(analytic/r["figure4_solid_r_cm_s"]-1),
            solver_vs_analytic_percent=100*(r["r_cm_s"]/analytic-1),
            corrected_surface_temperature_K=r["surface_temperature_K"], analytic_surface_temperature_K=Ts,
            late_speed_drift_percent=r["late_speed_drift_percent"],
            incident_flux_error_percent=r["incident_flux_error_percent"], name=r["name"]))
        if q == 200:
            properties.append({k: material[k] for k in ("ap_volume_fraction", "ap_mass_fraction",
                "density_kg_m3", "cp_J_kg_K", "conductivity_W_m_K", "heat_release_J_kg",
                "A_m_s", "activation_temperature_K")})
    write_csv(HERE / "comparison.csv", comparisons)
    write_csv(HERE / "mixed_properties.csv", properties)
    stats = []
    for q in (*fluxes, "all"):
        selected = [r for r in comparisons if q == "all" or r["q_cal_cm2_s"] == q]
        current = np.array([r["corrected_error_percent"] for r in selected])
        old = np.array([r["old_error_percent"] for r in selected])
        stats.append(dict(q_cal_cm2_s=q, n=len(selected), old_mean_absolute_error_percent=float(np.abs(old).mean()),
            corrected_mean_absolute_error_percent=float(np.abs(current).mean()),
            corrected_max_absolute_error_percent=float(np.abs(current).max()),
            corrected_mean_signed_error_percent=float(current.mean()),
            corrected_rms_relative_error_percent=float(np.sqrt(np.mean(current**2)))))
    write_csv(HERE / "error_summary.csv", stats)
    windows = []
    for r in all_rows:
        path = args.summary.parent / f"timeseries_{r['name']}.csv"
        raw = [{k: v if k == "plotfile" else float(v) if v else float("nan")
                for k, v in item.items()} for item in csv.DictReader(path.open())]
        for fraction in (.2, .3, .4, .5):
            s = summarize(raw, fraction)
            windows.append(dict(name=r["name"], relaxations=r["relaxations"], late_fraction=fraction,
                rate_cm_s=s["r_cm_s"], change_from_default_percent=100*(s["r_cm_s"]/r["r_cm_s"]-1),
                n_late=s["n_late"], fit_standard_error_cm_s=s["fit_standard_error_cm_s"]))
    write_csv(HERE / "fit_window_sensitivity.csv", windows)
    short = next(r for r in rows if r["q_cal_cm2_s"] == 500 and r["ap_volume_fraction"] == .4)
    long = long_rows[0]
    duration = dict(q_cal_cm2_s=500, ap_volume_fraction=.4,
        short_duration_s=short["end_time_s"], long_duration_s=long["end_time_s"],
        short_r_cm_s=short["r_cm_s"], long_r_cm_s=long["r_cm_s"],
        rate_change_percent=100*(long["r_cm_s"]/short["r_cm_s"]-1),
        short_late_drift_percent=short["late_speed_drift_percent"],
        long_late_drift_percent=long["late_speed_drift_percent"])
    (HERE / "duration_comparison.json").write_text(json.dumps(duration, indent=2) + "\n")
    curve_path = STUDY / "analysis/figure4_curves.csv"
    marker_path = STUDY / "analysis/figure4_markers.csv"
    curves, markers = [list(csv.DictReader(path.open())) for path in (curve_path, marker_path)]
    grid = np.linspace(0, 1, 201)
    fig, ax = plt.subplots(figsize=(9.5, 6.2), layout="constrained")
    analytic_rows = []
    for q, color in zip(fluxes, ("#0072B2", "#D55E00", "#009E73")):
        pts = [r for r in curves if r["kind"] == "eq15_solid" and float(r["q_cal_cm2_s"]) == q]
        ax.plot([float(r["t"]) for r in pts], [float(r["r_cm_s"]) for r in pts],
                    color=color, lw=1.3, label=f"q = {q} cal cm⁻² s⁻¹")
        analytic = [closure(params, t, q, volume_density=True)[0] for t in grid]
        ax.plot(grid, analytic, "--", color=color, lw=1)
        analytic_rows += [dict(q_cal_cm2_s=q, ap_volume_fraction=float(t), r_cm_s=r)
                          for t, r in zip(grid, analytic)]
        for kind, marker in (("dns_2d_circle", "o"), ("dns_3d_asterisk", "*")):
            pts = [r for r in markers if r["kind"] == kind and float(r["q_cal_cm2_s"]) == q]
            ax.scatter([float(r["t"]) for r in pts], [float(r["r_cm_s"]) for r in pts],
                           marker=marker, s=15, facecolors="none", edgecolors=color, alpha=.45)
        selected = [r for r in comparisons if r["q_cal_cm2_s"] == q]
        x = [r["ap_volume_fraction"] for r in selected]
        ax.plot(x, [r["corrected_r_cm_s"] for r in selected], "s", color=color, ms=5)
        ax.plot(x, [r["old_r_cm_s"] for r in selected], "x", color=color, ms=4, alpha=.5)
        for r in selected:
            t = r["ap_volume_fraction"]
            ax.annotate(f"{r['corrected_error_percent']:+.2f}%".replace("-", "−"),
                (t, r["corrected_r_cm_s"]),
                xytext=(-7 if t == 1. else 7, 8 if t == 0. else -10),
                textcoords="offset points", ha="right" if t == 1. else "left",
                va="bottom" if t == 0. else "top", fontsize=8.5, color=color,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=1.2))
    ax.set(title="Chen Figure 4: corrected mixture density", ylabel="Regression speed (cm s⁻¹)",
           xlabel="AP volume fraction", ylim=(0, 3.5), xlim=(-.025, 1.025), xticks=fractions)
    ax.text(.015, .97, "Labels: corrected LowMach error relative to Chen (%)",
            transform=ax.transAxes, ha="left", va="top", fontsize=9, color=".3")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.legend(handles=[Line2D([], [], color=".25", ls=ls, marker=m, label=label)
        for ls, m, label in (("-", None, "Chen Fig. 4 solid curves (Eq. 15)"),
            ("--", None, "Corrected analytic closure"), ("none", "s", "Corrected LowMach"),
            ("none", "x", "Previous LowMach"),
            ("none", "o", "Chen 2D DNS"), ("none", "*", "Chen 3D DNS"))],
        loc="outside lower center", ncol=3, frameon=False, fontsize=8)
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"figure4_comparison.{ext}", dpi=220)
    write_csv(HERE / "analytic_curves.csv", analytic_rows)
    example = next(r for r in comparisons if r["q_cal_cm2_s"] == 500 and r["ap_volume_fraction"] == .4)
    worst = max(comparisons, key=lambda r: abs(r["corrected_error_percent"]))
    point_table = table(["q", "AP volume fraction", "Chen [cm/s]", "Previous [cm/s]", "Corrected [cm/s]", "Previous error [%]", "Corrected error [%]"],
        [[f"{r['q_cal_cm2_s']:g}", f"{r['ap_volume_fraction']:.1f}", f"{r['chen_r_cm_s']:.5f}",
          f"{r['old_r_cm_s']:.5f}", f"{r['corrected_r_cm_s']:.5f}", f"{r['old_error_percent']:+.3f}",
          f"{r['corrected_error_percent']:+.3f}"] for r in comparisons])
    stats_table = table(["q", "Previous mean absolute error [%]", "Corrected mean absolute error [%]", "Corrected maximum error [%]"],
        [[str(s["q_cal_cm2_s"]), f"{s['old_mean_absolute_error_percent']:.3f}",
          f"{s['corrected_mean_absolute_error_percent']:.3f}", f"{s['corrected_max_absolute_error_percent']:.3f}"] for s in stats])
    report = f'''# Figure 4 sweep with volume-additive mixture density

The density correction reduces the 18-case mean absolute error against Chen Figure 4 from **{stats[-1]['old_mean_absolute_error_percent']:.3f}% to {stats[-1]['corrected_mean_absolute_error_percent']:.3f}%**. At 40% AP and q=500 cal/(cm² s), the measured LowMach error changes from **{example['old_error_percent']:+.3f}% to {example['corrected_error_percent']:+.3f}%**. The corrected analytic prediction at that point has error **{example['analytic_error_percent']:+.3f}%**. The earlier approximately −1% claim referred to that analytic result; the actual corrected solver result is reported separately here.

The largest corrected absolute error is {abs(worst['corrected_error_percent']):.3f}% at q={worst['q_cal_cm2_s']:g}, AP volume fraction {worst['ap_volume_fraction']:g}. Parameters were frozen throughout; no A/E, cp, Q, conductivity or AP endpoint was refitted. All six pure-endpoint reruns reproduce their previous rates within relative tolerance 1e-9.

![Figure 4 comparison](figure4_comparison.png)

The single regression-speed chart labels every corrected LowMach square with its signed percentage error relative to Chen's solid curve, `100*(r_LowMach/r_Chen-1)`. The dashed curves show the corrected analytic closure. Previous LowMach results remain faint crosses.

{stats_table}

{point_table}

## Change and fixed assumptions

`PhaseChange.H` now uses `rho=(1-t)*rho_binder+t*rho_AP`, equivalently `1/rho=(1-w)/rho_binder+w/rho_AP`. Specific heat and Q remain mass-weighted. Consequently rho*cp and rho*Q equal the sums of constituent heat capacities and phase-change energies per initial volume. The study generator uses this density consistently for its coupled reference solution, initial solid density, interface scaling and thermal relaxation time. Tests verify additive volume, heat capacity, phase-change heat, product mass conservation and thermal coupling.

Binder rho=920 kg/m³, cp=2130 J/(kg K), k=0.213 W/(m K), Q=−66 cal/g, A=24.60833296 m/s and E/R=5568.850642 K remain fixed. AP retains rho=1950 kg/m³, cp=1297.90 J/(kg K), k=0.4186 W/(m K), Q=−100 cal/g, A=948 m/s and E/R=11000 K. T0=300 K. E/R and ln(A) are volume-weighted; conductivity retains the existing Chen two-dimensional rule. Q is parsed as energy per mass and multiplied by transferred mass.

Gas chemistry remains frozen and the only phase-change mechanism is binder regression. The gas-slab source supplies the prescribed heat flux; negative Q absorbs heat. Surface temperature is computed from the coupled energy/kinetic problem, rather than prescribed independently. The article comparator is the vector-extracted solid curve (Eq. 15), interpolated at each selected volume fraction. DNS markers are overlaid but are not treated as independent experimental observations. [Chen et al. (2002)](https://doi.org/10.1016/S1540-7489(02)80357-1).

## Runtime and numerical checks

All 18 sweep points were rerun with the corrected executable, including pure endpoints. The six-relaxation-time setup retains ell/delta=0.25, eight cells per ell, 32 transverse cells and dt scale=0.2. Maximum late-window rate drift is {max(abs(r['late_speed_drift_percent']) for r in rows):.4f}%. Changing the fitting window to the final 20%, 30%, 40% or 50% changes baseline rates by at most {max(abs(r['change_from_default_percent']) for r in windows if r['relaxations']==6):.4f}%. Rates are deterministic eta=0.5 position slopes, cross-checked against integrated solid-volume loss; maximum disagreement is {max(abs(100*(r['volume_r_cm_s']/r['r_cm_s']-1)) for r in rows):.4f}%.

At 40% AP and q=500, a separate twelve-relaxation-time run changes the rate from {duration['short_r_cm_s']:.8f} to {duration['long_r_cm_s']:.8f} cm/s ({duration['rate_change_percent']:+.4f}%). The deeper lower boundary maintains the final cold-boundary separation, so this is a combined duration/deep-boundary check. Its late-window drift is {duration['long_late_drift_percent']:+.4f}%. Only this representative point has doubled-duration verification. The longer result does not replace a sweep point.

Corrected solver/analytic errors span {min(r['solver_vs_analytic_percent'] for r in comparisons):+.3f}% to {max(r['solver_vs_analytic_percent'] for r in comparisons):+.3f}%; incident-flux errors span {min(r['incident_flux_error_percent'] for r in rows):+.3f}% to {max(r['incident_flux_error_percent'] for r in rows):+.3f}%. Matching the analytic density rule does not remove the gas heating/interface surrogate's numerical bias or the unfitted AP endpoint error. This is not a grid-convergence or experimental validation study. Figure stroke resolution remains about 0.005 cm/s.

## Execution provenance

Three `gpt-5.6-luna` agents launched the simulations, each owning one flux group. The q=500 agent also launched the duration check. The root/current model performed all raw-output extraction, audits, fit-window analysis, error calculations and plotting. A premature completion receipt for the duration check was removed while that run continued; incomplete data were excluded. Final acceptance requires a successful observed exit, matching input/executable hashes, a finalized solver log and completed requested duration. The corrected executable SHA-256 is `{manifest['binary_sha256']}`. Fixed parameter SHA-256 is `{manifest['parameters_sha256']}`.

The old executable and changed source/input files are retained under `baseline_snapshot/`; old sweep outputs and reports are preserved. `density_change.patch` isolates this correction from earlier workspace changes. Build and solver/unit-test logs are retained under `checks/`. `comparison.csv`, `error_summary.csv`, `mixed_properties.csv`, `fit_window_sensitivity.csv`, `duration_comparison.json`, execution audits and `provenance.json` hold the machine-readable evidence. See `REPRODUCE.md` for commands.
'''
    (HERE / "REPORT.md").write_text(report)
    paths = [args.summary, prior_path, parameters_path, BINDER / "fixed_properties.json",
        HERE / "launch_manifest.json", curve_path, marker_path,
        STUDY / "analysis/figure4_provenance.json", Path(__file__),
        STUDY / "run_study.py", STUDY / "analysis/compare_pure_htpb.py",
        STUDY / "analysis/report_calibration.py", STUDY / "analysis/extract_runs.py",
        STUDY.parent / "src/Model/Mechanism/PhaseChange.H",
        HERE / "density_change.patch", HERE / "REPRODUCE.md",
        *sorted((HERE / "checks").glob("*.log"))]
    (HERE / "provenance.json").write_text(json.dumps({str(p): digest(p) for p in paths}, indent=2) + "\n")
    print(stats_table)
    print(point_table)
    print("Doubled-duration change [%]:", duration["rate_change_percent"])


if __name__ == "__main__":
    main()
