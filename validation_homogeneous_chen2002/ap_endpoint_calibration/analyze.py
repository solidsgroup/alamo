#!/usr/bin/env python3
"""Audit and report the frozen AP fit and the new Figure 4 mixture sweep."""
import argparse
import csv
import json
from pathlib import Path

from audit import HERE, STUDY, audit, closure, digest
from compare_pure_htpb import table
from extract_runs import summarize, write_csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

R = 8.31446261815324
PRIOR = STUDY / "binder_q66_calibration/density_correction"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summaries", type=Path, nargs="+", required=True)
    args = p.parse_args()
    params_path = HERE / "parameters_frozen.json"
    params = json.loads(params_path.read_text())
    baseline = json.loads((HERE / "parameters_baseline.json").read_text())
    freeze = json.loads((HERE / "freeze_record.json").read_text())
    assert digest(params_path) == freeze["parameters_sha256"]
    assert digest(Path(freeze["source_summary"])) == freeze["summary_sha256"]
    manifest = json.loads((HERE / "launch_manifest.json").read_text())
    assert manifest["parameters_sha256"] == digest(params_path)
    assert manifest["binary_sha256"] == digest(STUDY.parent / "bin/lowmach-2d-clang++")
    all_rows = [r for path in args.summaries for r in json.loads(path.read_text())]
    rows = sorted([r for r in all_rows if r["relaxations"] == 6],
                  key=lambda r: (r["q_cal_cm2_s"], r["ap_volume_fraction"]))
    long_rows = [r for r in all_rows if r["relaxations"] == 12]
    fluxes, fractions = (200., 500., 1000.), (0., .2, .4, .6, .8, 1.)
    assert len(all_rows) == 19 and len(rows) == 18 and len(long_rows) == 1
    assert {(r["q_cal_cm2_s"], r["ap_volume_fraction"]) for r in rows} == {
        (q, t) for q in fluxes for t in fractions}
    assert {(r["q_cal_cm2_s"], r["ap_volume_fraction"]) for r in long_rows} == {(500., .4)}
    fitted_names = {r["name"] for r in freeze["execution_audit"]}
    assert {r["name"] for r in all_rows} == fitted_names | {r["name"] for r in manifest["cases"]}
    for case in manifest["cases"]:
        assert case["input_sha256"] == digest(STUDY / "runs" / case["name"] / "input")
    evidence = audit(rows, params)
    long_evidence = audit(long_rows, params, relaxations=12)
    assert max(abs(r["late_speed_drift_percent"]) for r in rows) < .2, "Sweep not sufficiently steady"
    for name, data in (("execution_audit.json", evidence), ("duration_execution_audit.json", long_evidence),
                       ("simulation_summary.json", rows)):
        (HERE / name).write_text(json.dumps(data, indent=2)+"\n")
    write_csv(HERE / "simulation_summary.csv", rows)
    prior_rows = json.loads((PRIOR / "simulation_summary.json").read_text())
    audit(prior_rows, baseline)
    prior = {(r["q_cal_cm2_s"], r["ap_volume_fraction"]): r for r in prior_rows}
    comparisons = []
    for r in rows:
        q, t = r["q_cal_cm2_s"], r["ap_volume_fraction"]
        old = prior[q, t]
        if t == 0.:
            assert np.isclose(old["r_cm_s"], r["r_cm_s"], rtol=1.e-9), "Pure binder changed"
        case = json.loads((STUDY / "runs" / r["name"] / "case.json").read_text())
        m = case["reference"]
        b, ap = params["binder"], params["AP"]
        rho = (1-t)*b["density_kg_m3"]+t*ap["density_kg_m3"]
        assert np.isclose(m["density_kg_m3"], rho, rtol=1.e-12)
        for key in ("cp_J_kg_K", "heat_release_J_kg"):
            assert np.isclose(rho*m[key], (1-t)*b["density_kg_m3"]*b[key]
                              +t*ap["density_kg_m3"]*ap[key], rtol=1.e-12)
        analytic, Ts, _ = closure(params, t, q, volume_density=True)
        comparisons.append(dict(q_cal_cm2_s=q, ap_volume_fraction=t,
            role="AP fit" if t == 1. and q != 500. else "AP held out" if t == 1. else
                 "unchanged binder control" if t == 0. else "unfitted mixture prediction",
            chen_r_cm_s=r["figure4_solid_r_cm_s"], previous_r_cm_s=old["r_cm_s"],
            calibrated_r_cm_s=r["r_cm_s"], previous_error_percent=old["error_figure4_solid_percent"],
            calibrated_error_percent=r["error_figure4_solid_percent"], analytic_r_cm_s=analytic,
            analytic_error_percent=100*(analytic/r["figure4_solid_r_cm_s"]-1),
            solver_vs_analytic_percent=100*(r["r_cm_s"]/analytic-1),
            surface_temperature_K=r["surface_temperature_K"], analytic_surface_temperature_K=Ts,
            late_speed_drift_percent=r["late_speed_drift_percent"], name=r["name"]))
    write_csv(HERE / "comparison.csv", comparisons)
    pure_ap = [r for r in comparisons if r["ap_volume_fraction"] == 1.]
    assert max(abs(r["calibrated_error_percent"]) for r in pure_ap if r["role"] == "AP fit") < 1.
    stats = []
    for group in (*fluxes, "all", "mixtures_only"):
        selected = [r for r in comparisons if group == "all" or
                    group == "mixtures_only" and 0 < r["ap_volume_fraction"] < 1 or r["q_cal_cm2_s"] == group]
        old = np.array([r["previous_error_percent"] for r in selected])
        new = np.array([r["calibrated_error_percent"] for r in selected])
        stats.append(dict(group=group, n=len(selected), previous_mean_absolute_error_percent=float(np.abs(old).mean()),
            calibrated_mean_absolute_error_percent=float(np.abs(new).mean()),
            previous_max_absolute_error_percent=float(np.abs(old).max()),
            calibrated_max_absolute_error_percent=float(np.abs(new).max())))
    write_csv(HERE / "error_summary.csv", stats)
    windows = []
    for r in all_rows:
        caches = [path.parent/f"timeseries_{r['name']}.csv" for path in args.summaries]
        cache = next(path for path in caches if path.exists())
        raw = [{k: v if k == "plotfile" else float(v) if v else float("nan")
                for k, v in row.items()} for row in csv.DictReader(cache.open())]
        for fraction in (.2, .3, .4, .5):
            s = summarize(raw, fraction)
            windows.append(dict(name=r["name"], relaxations=r["relaxations"], late_fraction=fraction,
                r_cm_s=s["r_cm_s"], change_from_default_percent=100*(s["r_cm_s"]/r["r_cm_s"]-1),
                n_late=s["n_late"]))
    write_csv(HERE / "fit_window_sensitivity.csv", windows)
    assert max(abs(r["change_from_default_percent"]) for r in windows if r["relaxations"] == 6) < .1, "Rate depends on fitting window"
    short = next(r for r in rows if r["q_cal_cm2_s"] == 500. and r["ap_volume_fraction"] == .4)
    long = long_rows[0]
    duration = dict(short_duration_s=short["end_time_s"], long_duration_s=long["end_time_s"],
        short_r_cm_s=short["r_cm_s"], long_r_cm_s=long["r_cm_s"],
        rate_change_percent=100*(long["r_cm_s"]/short["r_cm_s"]-1),
        long_late_drift_percent=long["late_speed_drift_percent"])
    (HERE / "duration_comparison.json").write_text(json.dumps(duration, indent=2)+"\n")
    curves = list(csv.DictReader((STUDY / "analysis/figure4_curves.csv").open()))
    markers = list(csv.DictReader((STUDY / "analysis/figure4_markers.csv").open()))
    fig, ax = plt.subplots(figsize=(9.5, 6.2), layout="constrained")
    grid = np.linspace(0, 1, 201)
    analytic_rows = []
    for q, color in zip(fluxes, ("#0072B2", "#D55E00", "#009E73")):
        pts = [r for r in curves if r["kind"] == "eq15_solid" and float(r["q_cal_cm2_s"]) == q]
        ax.plot([float(r["t"]) for r in pts], [float(r["r_cm_s"]) for r in pts], color=color,
                lw=1.3, label=f"q = {q:g} cal cm⁻² s⁻¹")
        analytic = [closure(params, t, q, volume_density=True)[0] for t in grid]
        analytic_rows += [dict(q_cal_cm2_s=q, ap_volume_fraction=float(t), r_cm_s=r) for t, r in zip(grid, analytic)]
        ax.plot(grid, analytic, "--", color=color, lw=1)
        for kind, marker in (("dns_2d_circle", "o"), ("dns_3d_asterisk", "*")):
            pts = [r for r in markers if r["kind"] == kind and float(r["q_cal_cm2_s"]) == q]
            ax.scatter([float(r["t"]) for r in pts], [float(r["r_cm_s"]) for r in pts],
                       marker=marker, s=15, facecolors="none", edgecolors=color, alpha=.45)
        selected = [r for r in comparisons if r["q_cal_cm2_s"] == q]
        ax.plot([r["ap_volume_fraction"] for r in selected], [r["previous_r_cm_s"] for r in selected],
                "x", color=color, ms=4, alpha=.5)
        for r in selected:
            t = r["ap_volume_fraction"]
            ax.plot(t, r["calibrated_r_cm_s"], "s", color=color, ms=5)
            ax.annotate(f"{r['calibrated_error_percent']:+.2f}%".replace("-", "−"),
                (t, r["calibrated_r_cm_s"]), xytext=(-7 if t == 1. else 7, 8 if t == 0. else -10),
                textcoords="offset points", ha="right" if t == 1. else "left",
                va="bottom" if t == 0. else "top", fontsize=8.5, color=color,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=1.2))
    ax.set(title="Chen Figure 4: pure-AP kinetics recalibrated", ylabel="Regression speed (cm s⁻¹)",
           xlabel="AP volume fraction", ylim=(0, 3.5), xlim=(-.025, 1.025), xticks=fractions)
    ax.text(.015, .97, "Labels: LowMach error relative to Chen (%)", transform=ax.transAxes,
            ha="left", va="top", fontsize=9, color=".3")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    handles = [Line2D([], [], color=".25", ls=ls, marker=m, label=label,
               markerfacecolor="none" if m in ("o", "*") else ".25") for ls, m, label in (
        ("-", None, "Chen Fig. 4 solid curves (Eq. 15)"), ("--", None, "Recalibrated analytic closure"),
        ("none", "s", "LowMach: AP recalibrated"), ("none", "x", "LowMach: before AP fit"),
        ("none", "o", "Chen 2D DNS"), ("none", "*", "Chen 3D DNS"))]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=8)
    for ext in ("png", "pdf"):
        fig.savefig(HERE/f"figure4_comparison.{ext}", dpi=220)
    write_csv(HERE / "analytic_curves.csv", analytic_rows)
    ap, ap0 = params["AP"], baseline["AP"]
    analytic_params = json.loads((HERE / "parameters_analytic.json").read_text())["AP"]
    kinetics = [dict(stage=name, A_m_s=m["A_m_s"], activation_temperature_K=m["activation_temperature_K"],
                E_kJ_mol=R*m["activation_temperature_K"]/1000) for name, m in
                (("previous", ap0), ("direct analytic fit", analytic_params), ("final solver-calibrated", ap))]
    write_csv(HERE / "kinetics.csv", kinetics)
    stats_table = table(["Flux/group", "N", "Previous mean absolute error [%]", "New mean absolute error [%]", "New maximum absolute error [%]"],
        [[r["group"], r["n"], f"{r['previous_mean_absolute_error_percent']:.4f}",
          f"{r['calibrated_mean_absolute_error_percent']:.4f}", f"{r['calibrated_max_absolute_error_percent']:.4f}"] for r in stats])
    point_table = table(["q", "AP fraction", "Use", "Chen [cm/s]", "New LowMach [cm/s]", "Previous error [%]", "New error [%]"],
        [[f"{r['q_cal_cm2_s']:g}", f"{r['ap_volume_fraction']:.1f}", r["role"], f"{r['chen_r_cm_s']:.6f}",
          f"{r['calibrated_r_cm_s']:.6f}", f"{r['previous_error_percent']:+.4f}", f"{r['calibrated_error_percent']:+.4f}"]
         for r in comparisons])
    total = next(r for r in stats if r["group"] == "all")
    blend = next(r for r in stats if r["group"] == "mixtures_only")
    held = next(r for r in pure_ap if r["role"] == "AP held out")
    report = f'''# Pure-AP Arrhenius calibration and Figure 4 sweep

The new pure-AP fit gives **A_AP={ap['A_m_s']:.9g} m/s**, **E_AP/R={ap['activation_temperature_K']:.9g} K**, or **E_AP={R*ap['activation_temperature_K']/1000:.8g} kJ/mol**. Previously A_AP={ap0['A_m_s']:g} m/s and E_AP/R={ap0['activation_temperature_K']:g} K. The physical law is `r=A*exp[-(E/R)/Ts]`; `run_study.py` converts A to the phase-field multiplier for each interface width.

Across all 18 points, mean absolute error versus Chen changes from **{total['previous_mean_absolute_error_percent']:.4f}% to {total['calibrated_mean_absolute_error_percent']:.4f}%**. For the 12 unfitted interior mixtures alone it changes from **{blend['previous_mean_absolute_error_percent']:.4f}% to {blend['calibrated_mean_absolute_error_percent']:.4f}%**. The held-out pure-AP q=500 case has error **{held['calibrated_error_percent']:+.4f}%**.

![Figure 4 comparison](figure4_comparison.png)

The single chart labels each new LowMach square with signed error `100*(r_LowMach/r_Chen-1)`. Solid curves are vector-extracted Chen Eq. 15 curves; dashed curves are our coupled analytic closure with the new AP parameters. Faint crosses are the previous volume-additive-density sweep. DNS circles/asterisks retain the paper's values. The reference curves and markers are unchanged.

{stats_table}

{point_table}

## Calibration and fixed assumptions

Only AP A and E/R were adjusted using the pure-AP endpoints at q=200 and 1000 cal/(cm² s). The q=500 pure-AP point and every interior mixture were excluded from the fit. Existing completed pure-AP runs supplied the initial measured solver/analytic ratios. For each target rate, the same established inversion computes `Ts=T0+[q/(rho*r)+Q]/cp`, then solves `ln(r)=ln(A)-(E/R)/Ts`; the solver-corrected target is the article rate divided by its previously measured solver/analytic ratio. Fresh endpoint runs independently verify the new fit. `freeze_record.json` records acceptance before the held-out case and mixture sweep were generated. The two fresh accepted fit runs are reused as the sweep's q=200/1000 pure-AP endpoints; all other 16 sweep points were run with the frozen file.

The direct analytic fit is A_AP={analytic_params['A_m_s']:.9g} m/s and E_AP/R={analytic_params['activation_temperature_K']:.9g} K (`parameters_analytic.json`). The final solver-calibrated pair is conditional on the executed grid, interface and gas heating surrogate; it includes compensation for numerical bias.

The complete binder property dictionary and its calibration history are unchanged: rho=920 kg/m³, cp=2130 J/(kg K), k=0.213 W/(m K), Q=−66 cal/g, A=24.6083329584 m/s, E/R=5568.85064240 K. AP thermal properties remain rho=1950 kg/m³, cp=1297.9 J/(kg K), k=0.4186 W/(m K), Q=−100 cal/g. T0=300 K. Density is volume-additive; cp and Q remain mass-weighted; ln(A) and E/R remain volume-weighted; conductivity retains the existing Chen two-dimensional rule. Three fresh pure-binder controls reproduce the previous rates within relative tolerance 1e-9.

Gas chemistry is frozen. The prescribed gas-slab heating and negative phase-change Q are the configured thermal sources; there is no gas-reaction heat. Surface temperature is solved through the thermal/kinetic coupling. These are comparisons to the supplied article's model curves and DNS, not experimental material measurements. [Chen et al. (2002)](https://doi.org/10.1016/S1540-7489(02)80357-1).

## Runtime and numerical checks

Every sweep run reaches six thermal relaxation times, with ell/delta=0.25, eight cells per ell, 32 transverse cells and dt scale=0.2. Maximum late-window speed drift is {max(abs(r['late_speed_drift_percent']) for r in rows):.4f}%. Changing the rate-fit window from the final 40% to the final 20%, 30% or 50% changes rates by at most {max(abs(r['change_from_default_percent']) for r in windows if r['relaxations']==6):.4f}%. Rates come from raw eta=0.5 positions; independent volume-loss rates agree within {max(abs(100*(r['volume_r_cm_s']/r['r_cm_s']-1)) for r in rows):.4f}%.

At q=500 and 40% AP, a new twelve-relaxation-time run changes the measured rate from {duration['short_r_cm_s']:.8f} to {duration['long_r_cm_s']:.8f} cm/s ({duration['rate_change_percent']:+.4f}%). Its late drift is {duration['long_late_drift_percent']:+.4f}%. The deeper lower boundary preserves final cold-boundary separation, making this a combined duration/depth check; only this representative point was doubled. It does not replace a sweep point.

Solver/analytic rate differences span {min(r['solver_vs_analytic_percent'] for r in comparisons):+.4f}% to {max(r['solver_vs_analytic_percent'] for r in comparisons):+.4f}%; incident-flux errors span {min(r['incident_flux_error_percent'] for r in rows):+.4f}% to {max(r['incident_flux_error_percent'] for r in rows):+.4f}%. Fit scatter is not a statistical uncertainty estimate. This does not establish grid convergence. Figure stroke resolution is approximately 0.005 cm/s; extra reported digits preserve reproducibility.

## Reproduction and provenance

Cheaper `gpt-5.6-luna` agents launched the simulations and recorded observed terminal exits; the root/current model performed the calibration, raw-output extraction, audits, statistics and plotting. All accepted runs have matching input/executable hashes, successful receipts, finalized logs, full duration, fixed thermal/binder properties and independently recomputed closures. No solver source or executable changed for this follow-up.

`parameters_baseline.json` preserves the previous constituent values; `parameters_frozen.json`, `freeze_record.json`, `kinetics.csv`, `comparison.csv`, `error_summary.csv`, the execution audits, `fit_window_sensitivity.csv`, `duration_comparison.json` and `provenance.json` record this study. Earlier reports and figures remain intact. See `REPRODUCE.md` for commands.
'''
    (HERE / "REPORT.md").write_text(report)
    provenance_file = STUDY / "analysis/figure4_provenance.json"
    source = json.loads(provenance_file.read_text())
    assert digest(Path(source["source"])) == source["sha256"]
    paths = [params_path, HERE / "parameters_baseline.json", HERE / "parameters_analytic.json",
        HERE / "freeze_record.json", HERE / "launch_manifest.json", HERE / "REPRODUCE.md", Path(__file__), HERE / "audit.py",
        HERE / "calibrate.py", HERE / "freeze.py", *args.summaries, PRIOR / "simulation_summary.json",
        STUDY / "run_study.py", STUDY / "analysis/compare_study.py", STUDY / "analysis/extract_runs.py",
        STUDY / "analysis/compare_pure_htpb.py", provenance_file, Path(source["source"]),
        STUDY / "analysis/figure4_curves.csv", STUDY / "analysis/figure4_markers.csv"]
    (HERE / "provenance.json").write_text(json.dumps({str(path): digest(path) for path in paths}, indent=2)+"\n")
    print(stats_table)
    print(point_table)
    print("Duration change [%]:", duration["rate_change_percent"])


if __name__ == "__main__":
    main()
