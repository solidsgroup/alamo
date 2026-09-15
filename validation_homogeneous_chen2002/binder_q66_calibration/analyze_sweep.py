#!/usr/bin/env python3
"""Compare the frozen pure-constituent volume-fraction sweep with Chen Fig. 4."""
import argparse
import csv
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
from matplotlib.lines import Line2D
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, nargs="+", required=True)
    args = parser.parse_args()
    parameters_path = HERE / "parameters_frozen.json"
    parameters = json.loads(parameters_path.read_text())
    rows = sorted([r for path in args.summary for r in json.loads(path.read_text())],
                  key=lambda r: (r["q_cal_cm2_s"], r["ap_volume_fraction"]))
    fluxes, fractions = (200, 500, 1000), (0., .2, .4, .6, .8, 1.)
    assert len(rows) == 18
    assert {(r["q_cal_cm2_s"], r["ap_volume_fraction"]) for r in rows} == {
        (q, t) for q in fluxes for t in fractions}
    evidence = audit(rows, parameters, arrhenius=True,
                     fixed_parameters=HERE / "fixed_properties.json")
    out = HERE / "sweep"
    out.mkdir(exist_ok=True)
    (out / "execution_audit.json").write_text(json.dumps(evidence, indent=2) + "\n")
    (out / "simulation_summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    write_csv(out / "simulation_summary.csv", rows)
    comparison = []
    for r in rows:
        q, t = r["q_cal_cm2_s"], r["ap_volume_fraction"]
        analytic, T, _ = closure(parameters, t, q)
        comparison.append(dict(q_cal_cm2_s=q, ap_volume_fraction=t,
            ap_mass_fraction=(t*parameters["AP"]["density_kg_m3"] /
                (t*parameters["AP"]["density_kg_m3"]+(1-t)*parameters["binder"]["density_kg_m3"])),
            role="binder fit" if t == 0 and q != 500 else "binder held out" if t == 0
                else "fixed AP prediction" if t == 1 else "blend prediction",
            lowmach_r_cm_s=r["r_cm_s"], chen_r_cm_s=r["figure4_solid_r_cm_s"],
            error_cm_s=r["r_cm_s"]-r["figure4_solid_r_cm_s"],
            error_percent=r["error_figure4_solid_percent"],
            analytic_r_cm_s=analytic, solver_error_vs_analytic_percent=100*(r["r_cm_s"]/analytic-1),
            surface_temperature_K=r["surface_temperature_K"],
            late_speed_drift_percent=r["late_speed_drift_percent"],
            incident_flux_error_percent=r["incident_flux_error_percent"], name=r["name"]))
    write_csv(out / "comparison.csv", comparison)
    mixed_properties = []
    for t in fractions:
        row = next(r for r in rows if r["q_cal_cm2_s"] == 200 and r["ap_volume_fraction"] == t)
        material = json.loads((STUDY / "runs" / row["name"] / "case.json").read_text())["reference"]
        mixed_properties.append({k: material[k] for k in (
            "ap_volume_fraction", "ap_mass_fraction", "density_kg_m3", "cp_J_kg_K",
            "conductivity_W_m_K", "heat_release_J_kg", "A_m_s", "activation_temperature_K")})
    write_csv(out / "mixed_properties.csv", mixed_properties)
    stats = []
    for q in (*fluxes, "all"):
        selected = [r for r in comparison if q == "all" or r["q_cal_cm2_s"] == q]
        errors = np.array([r["error_percent"] for r in selected])
        stats.append(dict(q_cal_cm2_s=q, n=len(selected), mean_signed_error_percent=float(errors.mean()),
            mean_absolute_error_percent=float(np.abs(errors).mean()),
            rms_relative_error_percent=float(np.sqrt(np.mean(errors**2))),
            max_absolute_error_percent=float(np.abs(errors).max())))
    write_csv(out / "error_summary.csv", stats)
    curves_path = STUDY / "analysis/figure4_curves.csv"
    markers_path = STUDY / "analysis/figure4_markers.csv"
    curves = list(csv.DictReader(curves_path.open()))
    markers = list(csv.DictReader(markers_path.open()))
    fig, axs = plt.subplots(2, 3, figsize=(11.5, 6.4), layout="constrained",
                            sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    grid = np.linspace(0, 1, 201)
    analytic_rows = []
    for i, q in enumerate(fluxes):
        ax, residual = axs[:, i]
        pts = [r for r in curves if r["kind"] == "eq15_solid" and float(r["q_cal_cm2_s"]) == q]
        ax.plot([float(r["t"]) for r in pts], [float(r["r_cm_s"]) for r in pts],
                color="black", lw=1.4, label="Chen Fig. 4 solid curve")
        for kind, marker, label in (("dns_2d_circle", "o", "Chen 2D DNS"),
                                    ("dns_3d_asterisk", "*", "Chen 3D DNS")):
            subset = [r for r in markers if r["kind"] == kind and float(r["q_cal_cm2_s"]) == q]
            ax.scatter([float(r["t"]) for r in subset], [float(r["r_cm_s"]) for r in subset],
                       marker=marker, s=19, facecolors="none", edgecolors=".55", label=label)
        analytic = [closure(parameters, t, q)[0] for t in grid]
        ax.plot(grid, analytic, "--", color="#0072B2", lw=1.1, label="Frozen analytic closure")
        analytic_rows += [dict(q_cal_cm2_s=q, ap_volume_fraction=float(t), r_cm_s=r)
                          for t, r in zip(grid, analytic)]
        selected = [r for r in comparison if r["q_cal_cm2_s"] == q]
        x = [r["ap_volume_fraction"] for r in selected]
        ax.plot(x, [r["lowmach_r_cm_s"] for r in selected], "s", color="#D55E00",
                ms=5, label="LowMach, frozen constituents")
        residual.plot(x, [r["error_percent"] for r in selected], "s-", color="#D55E00", ms=4)
        residual.axhline(0, color=".4", lw=.7)
        ax.set(title=f"q = {q} cal cm⁻² s⁻¹")
        residual.set(xlabel="AP volume fraction", xlim=(-.025, 1.025), xticks=fractions)
    axs[0, 0].set_ylabel("Regression speed (cm s⁻¹)")
    axs[1, 0].set_ylabel("Error vs Chen (%)")
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False, fontsize=8)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure4_comparison.{ext}", dpi=220)
    fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
    for q, color in zip(fluxes, ("#0072B2", "#D55E00", "#009E73")):
        pts = [r for r in curves if r["kind"] == "eq15_solid" and float(r["q_cal_cm2_s"]) == q]
        ax.plot([float(r["t"]) for r in pts], [float(r["r_cm_s"]) for r in pts],
                color=color, lw=1.4, label=f"q = {q} cal cm⁻² s⁻¹")
        ax.plot(grid, [closure(parameters, t, q)[0] for t in grid], "--", color=color, lw=1)
        for kind, marker in (("dns_2d_circle", "o"), ("dns_3d_asterisk", "*")):
            subset = [r for r in markers if r["kind"] == kind and float(r["q_cal_cm2_s"]) == q]
            ax.scatter([float(r["t"]) for r in subset], [float(r["r_cm_s"]) for r in subset],
                       marker=marker, s=16, facecolors="none", edgecolors=color, alpha=.55)
        selected = [r for r in comparison if r["q_cal_cm2_s"] == q]
        ax.plot([r["ap_volume_fraction"] for r in selected], [r["lowmach_r_cm_s"] for r in selected],
                "s", color=color, ms=5)
    ax.set(xlabel="AP volume fraction", ylabel="Regression speed (cm s⁻¹)",
           title="Chen Figure 4 and frozen-constituent LowMach sweep", xlim=(-.02, 1.02), ylim=(0, 3.5))
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    fig.legend(handles=[Line2D([], [], color=".25", ls=ls, marker=marker, label=label)
        for ls, marker, label in (("-", None, "Chen solid curves"), ("--", None, "Frozen analytic closure"),
                                  ("none", "s", "LowMach"), ("none", "o", "Chen 2D DNS"),
                                  ("none", "*", "Chen 3D DNS"))],
        loc="outside lower center", frameon=False, fontsize=8, ncol=3)
    for ext in ("pdf", "png"):
        fig.savefig(out / f"figure4_overlay.{ext}", dpi=220)
    write_csv(out / "analytic_curves.csv", analytic_rows)
    paths = [parameters_path, HERE / "fixed_properties.json", *args.summary, curves_path,
             markers_path, STUDY / "run_study.py", Path(__file__)]
    (out / "provenance.json").write_text(json.dumps({str(p): digest(p) for p in paths}, indent=2) + "\n")
    comparison_table = table(["q [cal/(cm² s)]", "AP volume fraction", "LowMach [cm/s]", "Chen [cm/s]", "Error [%]"],
        [[f"{r['q_cal_cm2_s']:g}", f"{r['ap_volume_fraction']:.1f}",
          f"{r['lowmach_r_cm_s']:.5f}", f"{r['chen_r_cm_s']:.5f}",
          f"{r['error_percent']:+.3f}"] for r in comparison])
    stats_table = table(["q", "Mean signed error [%]", "Mean absolute error [%]", "RMS relative error [%]", "Max absolute error [%]"],
        [[str(r["q_cal_cm2_s"]), f"{r['mean_signed_error_percent']:+.3f}",
          f"{r['mean_absolute_error_percent']:.3f}", f"{r['rms_relative_error_percent']:.3f}",
          f"{r['max_absolute_error_percent']:.3f}"] for r in stats])
    worst = max(comparison, key=lambda r: abs(r["error_percent"]))
    drift = max(abs(r["late_speed_drift_percent"]) for r in rows)
    report = f'''# Frozen-constituent sweep compared with Chen Figure 4

All 18 cases completed: AP volume fractions 0, 0.2, 0.4, 0.6, 0.8 and 1 at q=200, 500 and 1000 cal/(cm² s). The three pure-binder cases are reused from the calibration/held-out assessment; the remaining 15 cases are new predictions with frozen constituent properties. No mixture data were used to change any parameter.

Across all 18 cases, mean absolute relative error against Chen's Figure 4 solid curves is **{stats[-1]['mean_absolute_error_percent']:.3f}%**; the largest absolute error is **{abs(worst['error_percent']):.3f}%**, at q={worst['q_cal_cm2_s']:g}, AP volume fraction {worst['ap_volume_fraction']:g}, where LowMach gives {worst['lowmach_r_cm_s']:.5f} cm/s versus Chen {worst['chen_r_cm_s']:.5f} cm/s. Signed error is `100*(r_LowMach/r_Chen - 1)`; negative values mean underprediction.

![Figure 4 comparison](figure4_comparison.png)

`figure4_overlay.pdf` / `.png` reproduce the original single-axis arrangement
of all three flux curves, with the new simulation points overlaid.

## Pointwise errors

{comparison_table}

{stats_table}

The comparator is linear interpolation of the vector-extracted **solid curves** in Figure 4 (equation 15), including exact pure-material endpoints. It is not an error calculation against individual scattered DNS markers. Both 2D and 3D DNS markers from the paper are shown in the overlay. These are published model/DNS results, not experimental burn-rate data. One printed stroke is approximately 0.005 cm/s; the displayed numerical precision does not imply equivalent source accuracy. See [Chen et al. (2002)](https://doi.org/10.1016/S1540-7489(02)80357-1).

## Frozen properties and mixing

Binder: rho=920 kg/m³, cp=2130 J/(kg K), k=0.213 W/(m K), Q=−276144 J/kg (−66 cal/g), A={parameters['binder']['A_m_s']:.8g} m/s, E/R={parameters['binder']['activation_temperature_K']:.8g} K. AP retains `reference/pure_htpb.json`: rho=1950 kg/m³, cp=1297.90 J/(kg K), k=0.4186 W/(m K), Q=−418400 J/kg, A=948 m/s, E/R=11000 K. T0=300 K.

The existing homogeneous implementation converts AP volume fraction t to mass fraction `w=rho_AP*t/[rho_AP*t+rho_binder*(1-t)]`, then uses mass-weighted density, cp and Q. Activation temperature is volume-weighted and ln(A) is volume-weighted. Conductivity follows the existing two-dimensional Chen mixing formula. Pure constituent values are supplied once. The density mixing convention is inherited from the existing study and implementation; no mixing rule is refitted here. The arithmetic mass-weighted density differs from ordinary volume-additive mixture density, so these errors characterize this implemented closure.

The blue dashed curves are independent analytic evaluations of equations (14)/(19) using the frozen mixed properties. Their separation from the black Chen curves measures closure/constituent differences. The additional difference between the LowMach squares and blue curves measures the numerical heating/interface surrogate. `comparison.csv` records both errors separately. Matching pure binder does not force agreement for mixtures or the unfitted pure AP endpoint.

## Numerical checks and reproducibility

All cases use the inherited settings: ell/delta=0.25, eight cells per ell, 32 transverse cells, dt scale=0.2, six relaxation times, gas k=100 W/(m K), gas cp=1000 J/(kg K), pressure 100 MPa, prescribed gas-slab heating, frozen chemistry, and disabled temperature advection. Rates are measured from raw eta=0.5 interface-position fits over the final 40% of each simulation, cross-checked against volume loss.

Maximum late-window speed drift is {drift:.4f}%; the 0.2% steadiness criterion is {'met for all cases' if drift < .2 else 'not met for every case; inspect the tabulated drift before interpreting steady rates'}. Maximum disagreement between volume-loss and interface-position rates is {max(abs(100*(r['volume_r_cm_s']/r['r_cm_s']-1)) for r in rows):.4f}%. Incident-flux errors span {min(r['incident_flux_error_percent'] for r in rows):+.3f}% to {max(r['incident_flux_error_percent'] for r in rows):+.3f}%. Solver/analytic rate errors span {min(r['solver_error_vs_analytic_percent'] for r in comparison):+.3f}% to {max(r['solver_error_vs_analytic_percent'] for r in comparison):+.3f}%. These checks do not establish grid convergence or validate an experimental pressure-dependent gas flame.

`execution_audit.json` verifies completed durations, successful solver exits, input/binary hashes, unchanged constituent parameters and independently recomputed analytic targets for every case. The frozen parameter file is `../parameters_frozen.json`, SHA-256 `{digest(parameters_path)}`. `comparison.csv` contains pointwise results; `error_summary.csv` contains error aggregates; `mixed_properties.csv` lists the thermal and Arrhenius properties at each composition; `figure4_comparison.pdf` is the exportable overlay.

See `../REPRODUCE.md` for the calibration and sweep commands.
'''
    (out / "SWEEP_REPORT.md").write_text(report)
    print(stats_table)
    print(comparison_table)


if __name__ == "__main__":
    main()
