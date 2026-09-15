#!/usr/bin/env python3
"""Snapshot nine completed baselines of the user-stopped pure-constituent study.

First run compare_study.py --tag _pure_htpb --output analysis/pure_htpb.
This independently recomputes the analytic closure, audits execution receipts,
and compares raw-output-derived rates with the original study and Chen Fig. 4.
"""
import csv
import hashlib
import json
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/chen-pure-htpb-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.optimize import brentq
from extract_runs import write_csv

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
OUT = HERE / "pure_htpb"


def closure(parameters, t, q, volume_density=False):
    """Independent SI evaluation of Chen Eqs. 14/19 and the Eq. 15 rate."""
    b, ap = parameters["binder"], parameters["AP"]
    w = t*ap["density_kg_m3"] / (t*ap["density_kg_m3"]+(1-t)*b["density_kg_m3"])
    material = {k: (1-w)*b[k]+w*ap[k] for k in
                ("density_kg_m3", "cp_J_kg_K", "heat_release_J_kg")}
    if volume_density:
        material["density_kg_m3"] = (1-t)*b["density_kg_m3"]+t*ap["density_kg_m3"]
    material["A_m_s"] = b["A_m_s"]**(1-t)*ap["A_m_s"]**t
    material["activation_temperature_K"] = (1-t)*b["activation_temperature_K"]+t*ap["activation_temperature_K"]

    def solve(m):
        rho, cp, Q, A, E = [m[k] for k in ("density_kg_m3", "cp_J_kg_K",
                              "heat_release_J_kg", "A_m_s", "activation_temperature_K")]
        T = brentq(lambda T: rho*A*math.exp(-E/T)*(cp*(T-parameters["T0_K"])-Q)-41840*q,
                   max(1., parameters["T0_K"]+Q/cp), 10000., xtol=1.e-10)
        return 100*A*math.exp(-E/T), T

    r, T = solve(material)
    rb, _ = solve(b)
    rap, _ = solve(ap)
    return r, T, 1/((1-t)/rb+t/rap)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "|" + "---|"*len(headers)] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def main():
    new = json.loads((OUT/"simulation_summary.json").read_text())
    old = json.loads((HERE/"simulation_summary.json").read_text())
    parameters = json.loads((STUDY/"reference/pure_htpb.json").read_text())
    old_parameters = json.loads((STUDY/"reference/parameters.json").read_text())
    cp_only = {**old_parameters, "binder": {**old_parameters["binder"],
                "cp_J_kg_K": parameters["binder"]["cp_J_kg_K"]}}
    heat_only = {**old_parameters, "binder": {**old_parameters["binder"],
                  "heat_release_J_kg": parameters["binder"]["heat_release_J_kg"]}}
    expected_names = {p.parent.name for p in (STUDY/"runs").glob("*_pure_htpb*/case.json")}
    assert len(expected_names) == 11, "Expected the original eleven planned cases"
    assert len(new) == 9 and all(s["name"].endswith("_pure_htpb_nx32") for s in new), "Stopped study report includes only the nine completed baselines"
    base = sorted([s for s in new if s["name"].endswith("_pure_htpb_nx32")], key=lambda s: (s["q_cal_cm2_s"], s["ap_volume_fraction"]))
    prior = {(s["q_cal_cm2_s"], s["ap_volume_fraction"]): s for s in old if s["name"].endswith("_nx32")}
    assert len(base) == 9 and len(prior) == 9
    binary_hash = digest(STUDY.parent/"bin/lowmach-2d-clang++")
    audit = []
    for s in new:
        folder = STUDY/"runs"/s["name"]
        case = json.loads((folder/"case.json").read_text())
        run = json.loads((folder/"run.json").read_text())
        input_hash = digest(folder/"input")
        assert run["returncode"] == 0 and input_hash == run["input_sha256"] == case["input_sha256"], s["name"]
        assert run["binary_sha256"] == binary_hash, s["name"]
        assert case["parameters"] == parameters, s["name"]
        assert s["end_time_s"] >= .999*case["duration_s"], s["name"]
        log = (folder/"output/out.log").read_text()
        steps = re.findall(r"STEP\s+(\d+) ends\. TIME = ([\d.eE+\-]+) DT = ([\d.eE+\-]+)", log)
        assert steps and "finalized" in log, s["name"]
        final_step, final_time, final_dt = map(float, steps[-1])
        assert final_time >= case["duration_s"] and abs(final_time-s["end_time_s"]) <= 1.01*final_dt, s["name"]
        r, T, r15 = closure(parameters, s["ap_volume_fraction"], s["q_cal_cm2_s"])
        assert np.isclose(r, s["expected_r_eq14_19_cm_s"], rtol=1.e-10)
        assert np.isclose(T, s["expected_surface_temperature_K"], rtol=1.e-10)
        audit.append(dict(name=s["name"], returncode=run["returncode"], input_sha256=input_hash,
                          binary_sha256=run["binary_sha256"], final_plot_time_s=s["end_time_s"],
                          binary_matches_current_executable=True,
                          requested_duration_s=case["duration_s"], input_hash_matches=True,
                          final_log_time_s=final_time, final_step=int(final_step), last_dt_s=final_dt,
                          material_matches=True, independently_recomputed_closure_matches=True,
                          elapsed_s=run.get("elapsed_s"), command=run["command"]))
    (OUT/"execution_audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    stopped = []
    for name in sorted(expected_names-{s["name"] for s in new}):
        folder = STUDY/"runs"/name
        receipt = json.loads((folder/"run.json").read_text()) if (folder/"run.json").exists() else None
        stopped.append(dict(name=name, status="stopped at user request; partial outputs excluded", receipt=receipt))
    (OUT/"stopped_cases.json").write_text(json.dumps(stopped,indent=2)+"\n")
    source_paths = [STUDY/"reference/pure_htpb.json", STUDY/"reference/parameters.json",
                    HERE/"simulation_summary.json", HERE/"figure4_curves.csv",
                    Path("/home/esandall/Downloads/homogeneous_model_paper.pdf"),
                    Path("/home/esandall/Downloads/gross_fourflame_model.pdf")]
    (OUT/"provenance.json").write_text(json.dumps({str(p):digest(p) for p in source_paths},indent=2)+"\n")
    comparison = []
    for s in base:
        q, t = s["q_cal_cm2_s"], s["ap_volume_fraction"]
        p = prior[q,t]
        r, T, r15 = closure(parameters, t, q)
        comparison.append(dict(q_cal_cm2_s=q, ap_volume_fraction=t, pure_htpb_r_cm_s=s["r_cm_s"],
            pure_htpb_closure_r_cm_s=r, pure_htpb_eq15_r_cm_s=r15,
            prior_r_cm_s=p["r_cm_s"], prior_closure_r_cm_s=p["expected_r_eq14_19_cm_s"],
            chen_figure4_r_cm_s=s["figure4_solid_r_cm_s"],
            numerical_error_percent=100*(s["r_cm_s"]/r-1),
            change_from_prior_percent=100*(s["r_cm_s"]/p["r_cm_s"]-1),
            closure_error_chen_percent=100*(r/s["figure4_solid_r_cm_s"]-1),
            simulation_error_chen_percent=s["error_figure4_solid_percent"]))
    write_csv(OUT/"parameter_comparison.csv", comparison)
    (OUT/"parameter_comparison.json").write_text(json.dumps(comparison, indent=2)+"\n")
    density_sensitivity = [dict(q_cal_cm2_s=q, ap_volume_fraction=.5,
        executed_density_rule="arithmetic mass-weighted", arithmetic_density_closure_cm_s=closure(parameters,.5,q)[0],
        volume_density_closure_cm_s=closure(parameters,.5,q,volume_density=True)[0],
        alternative_is_analytic_only=True) for q in (200,500,1000)]
    write_csv(OUT/"density_sensitivity.csv",density_sensitivity)
    curves = list(csv.DictReader((HERE/"figure4_curves.csv").open()))
    plt.rcParams.update({"font.size":9, "pdf.fonttype":42, "axes.spines.top":False, "axes.spines.right":False})
    fig, axs = plt.subplots(2,3,figsize=(12,6.8),layout="constrained",sharex=True,
                            gridspec_kw={"height_ratios":[2,1]})
    grid = np.linspace(0,1,201)
    analytic = []
    for i,q in enumerate((200,500,1000)):
        pts = [r for r in curves if r["kind"] == "eq15_solid" and float(r["q_cal_cm2_s"]) == q]
        axs[0,i].plot([float(p["t"]) for p in pts], [float(p["r_cm_s"]) for p in pts], color="black",lw=1.3)
        new_curve = [closure(parameters,t,q) for t in grid]
        old_curve = [closure(old_parameters,t,q) for t in grid]
        axs[0,i].plot(grid,[r[0] for r in new_curve],color="#0072B2",lw=1.4)
        axs[0,i].plot(grid,[r[2] for r in new_curve],color="#0072B2",ls=":",lw=1.2)
        axs[0,i].plot(grid,[r[0] for r in old_curve],color="#D55E00",ls="--",lw=1)
        for t,n,o in zip(grid,new_curve,old_curve):
            analytic.append(dict(q_cal_cm2_s=q, ap_volume_fraction=t, pure_htpb_eq14_19_r_cm_s=n[0],
                                 pure_htpb_eq15_r_cm_s=n[2], prior_eq14_19_r_cm_s=o[0]))
        for s in [s for s in base if s["q_cal_cm2_s"] == q]:
            t=s["ap_volume_fraction"]
            axs[0,i].plot(t,s["r_cm_s"],"o",color="#0072B2",ms=5)
            axs[0,i].plot(t,prior[q,t]["r_cm_s"],"x",color="#D55E00",ms=6)
        for s in [s for s in new if s["q_cal_cm2_s"] == q]:
            mark = "D" if s["width_ratio"] < .249 else "^" if s["dt_scale"] < .99 else "o"
            axs[1,i].plot(s["ap_volume_fraction"],s["error_eq14_19_percent"],marker=mark,color="#0072B2",ls="none",ms=5)
        axs[0,i].set(title=f"q = {q} cal cm⁻² s⁻¹")
        axs[1,i].axhline(0,color=".4",lw=.7)
        axs[1,i].set(xlabel="AP volume fraction, t", xlim=(-.03,1.03))
    error_limits = (min(0., min(s["error_eq14_19_percent"] for s in new))-.25,
                    max(0., max(s["error_eq14_19_percent"] for s in new))+.25)
    for ax in axs[1]: ax.set_ylim(*error_limits)
    axs[0,0].set(ylabel="Regression speed (cm s⁻¹)")
    axs[1,0].set(ylabel="New run error vs\nnew Eq. 14+19 (%)")
    handles = [Line2D([],[],color=c,ls=ls,marker=m,label=label) for c,ls,m,label in [
        ("black","-",None,"Chen Fig. 4 solid curve"), ("#0072B2","-",None,"Pure-HTPB Eq. 14+19"),
        ("#0072B2",":",None,"Pure-HTPB Eq. 15"), ("#D55E00","--",None,"Prior fitted Eq. 14+19"),
        ("#0072B2","none","o","New LowMach baseline"), ("#D55E00","none","x","Prior LowMach baseline"),
        ]]
    fig.legend(handles=handles,loc="outside lower center",ncol=4,frameon=False,fontsize=8)
    fig.suptitle("Stopped pure-constituent study: nine completed baselines")
    for suffix in ("pdf","png"): fig.savefig(OUT/f"parameter_comparison.{suffix}",dpi=220)
    write_csv(OUT/"analytic_curves.csv",analytic)
    key = sorted([s for s in new if s["q_cal_cm2_s"] == 500 and s["ap_volume_fraction"] == .5],
                 key=lambda s:(s["width_ratio"] < .249,s["dt_scale"] < .99))
    labels = ["Baseline","Smaller step","Thinner interface + smaller step"]
    errors = [s["error_eq14_19_percent"] for s in base]
    material_table = table(["Constituent","ρ [kg/m³]","cp [J/(kg K)]","k [W/(m K)]","Q [J/kg]","A [m/s]","E/R [K]"],
        [[label]+[f"{m[k]:g}" for k in ("density_kg_m3","cp_J_kg_K","conductivity_W_m_K","heat_release_J_kg","A_m_s","activation_temperature_K")]
         for label,m in [("HTPB",parameters["binder"]),("AP",parameters["AP"])]])
    baseline_table = table(["q [cal/(cm² s)]","t","New LowMach","New Eq. 14+19","Error [%]","Prior LowMach","Chen Fig. 4","New vs Chen [%]"],
        [[f"{s['q_cal_cm2_s']:g}",f"{s['ap_volume_fraction']:.1f}"]+[f"{s[k]:.5f}" for k in
          ("pure_htpb_r_cm_s","pure_htpb_closure_r_cm_s")]+[f"{s['numerical_error_percent']:+.2f}"]+
          [f"{s[k]:.5f}" for k in ("prior_r_cm_s","chen_figure4_r_cm_s")]+[f"{s['simulation_error_chen_percent']:+.2f}"] for s in comparison])
    refinement_table = table(["Setup","ell/delta","dy [µm]","Rate [cm/s]","Closure error [%]","Late drift [%]"],
        [[label,f"{s['width_ratio']:.3f}",f"{1e6*s['dx_m']:.4f}",f"{s['r_cm_s']:.5f}",f"{s['error_eq14_19_percent']:+.3f}",f"{s['late_speed_drift_percent']:+.3f}"] for label,s in zip(labels,key)])
    thermal_table = table(["Setup","Incident flux error [% q]","Gas storage [% q]","Snapshot residual [% q]","Steady-solid residual [% q]"],
        [[label]+[f"{s[k]:+.3f}" for k in ("incident_flux_error_percent","gas_storage_percent_input","thermal_balance_residual_percent_input","steady_profile_balance_residual_percent_input")] for label,s in zip(labels,key)])
    report = f'''# STOPPED: pure-HTPB constituent rerun of the Chen heat-flux sweep

The user stopped this parameter study because its agreement with Chen Figure 4 became worse. Nine baseline flux/composition combinations completed; the two q=500, t=0.5 refinements were stopped and their partial outputs are excluded. No new simulations are requested by this report. The completed baselines differ from their own analytic Eq. 14+19 closure by {min(errors):+.2f}% to {max(errors):+.2f}%. This documents numerical consistency for the selected parameters, without a completed refinement assessment.

The user's new requested direction is calibration to experimental pure-HTPB data. Chen Figure 4 supplies model/DNS endpoints, while the supplied Gross paper's propellant experimental comparisons concern AP/HTPB mixtures. Neither is a directly interchangeable experimental pure-HTPB calibration target. The specific pure-HTPB dataset and its conditions must be identified before fitting; no experimental calibration is performed here.

The new material set is an AP-free HTPB constituent engineering baseline, with AP mixed in once through the homogeneous model. It is not fitted to Chen Figure 4 or independently validated for a particular cured HTPB formulation. Consequently, changed agreement with that figure measures changed material assumptions as well as numerical error. The prior fitted study remains intact in the parent analysis directory.

## Parameters and provenance

{material_table}

T0=300 K. `../../reference/pure_htpb.json` is the exact parameter record. Densities, heat capacities, and conductivities retain the original input's provisional pure-constituent values. Chen (2002), Figs. 5–6 supplies the Arrhenius endpoints. Gross et al. (2013), section 4.2 supplies the selected endothermic heats, −300 cal/g for HTPB and −100 cal/g for AP (1 cal/g=4184 J/kg). The sources are the user-supplied PDFs; no additional scientific references were used.

In the prior Figure-4-fitted realization, HTPB cp=1255.24 J/(kg K), Q=−196590 J/kg and k=0.271382 W/(m K). The new HTPB heat capacity is larger and its decomposition heat sink much stronger. At the same incident heat flux, those changes lower the rate; conductivity also changes the thermal length and numerical setup. AP parameters change much less. Neither dataset identifies uniquely the physical properties used by Chen.

At q=500 for pure HTPB, an analytic-only change of cp alone lowers the prior {closure(old_parameters,0,500)[0]:.5f} cm/s rate to {closure(cp_only,0,500)[0]:.5f}; changing Q alone gives {closure(heat_only,0,500)[0]:.5f}; changing both gives {closure(parameters,0,500)[0]:.5f} cm/s. Both thermal changes matter. Conductivity does not enter this steady surface energy/kinetic balance directly, although it changes the resolved thermal profile and its numerical errors.

The independent analytic calculation uses q=ρr[cp(Ts−T0)−Q], r=A exp[−(E/R)/Ts]. For a blend, AP mass fraction is w=ρ_AP t/[ρ_AP t+ρ_B(1−t)]; ρ, cp and Q use the existing arithmetic mass-weighted rule, while log A and E/R use volume fraction t. Eq. 15 instead harmonically averages the two pure-component rates at equal q. The density-convention ambiguity documented in `../REFERENCE_NOTES.md` remains; it is not resolved by changing the binder properties.

The change also reverses the pure-endpoint ordering. At q=500, the new analytic HTPB rate is {closure(parameters,0,500)[0]:.5f} cm/s and AP is {closure(parameters,1,500)[0]:.5f} cm/s, while the prior HTPB rate was {closure(old_parameters,0,500)[0]:.5f} cm/s. New HTPB therefore regresses more slowly than AP under the same prescribed heat flux. The new t=0.5 Eq. 14+19 rate, {closure(parameters,.5,500)[0]:.5f} cm/s, is below both pure endpoints. This interior depression comes from the chosen combined thermal/kinetic mixing closure, conditional on its arithmetic mass-weighted density rule; it is not a feature forced by Eq. 15, whose harmonic mean remains between the pure rates. With unequal constituent densities, arithmetic mass weighting of density raises blend volumetric heat capacity relative to the simple volume-weighted constituent capacity, contributing to slower mixed regression. These conditional model predictions should not be interpreted as independently validated physical trends or a correction to Chen's article.

An **analytic-only density sensitivity** holds cp, Q, A and E/R fixed at the same blend values and substitutes volume-weighted density. At q=500, t=0.5 this gives {closure(parameters,.5,500,volume_density=True)[0]:.5f} cm/s, between the pure endpoints, compared with {closure(parameters,.5,500)[0]:.5f} cm/s for the executed rule. The arithmetic density is 12.88% larger here. This isolates a substantial density contribution to the interior depression without attributing it solely to material kinetics. `density_sensitivity.csv` retains all three fluxes. No simulations use this alternative; the source and all inputs retain the user-specified arithmetic mass-weighted density. This comparison does not identify which convention Chen used.

## Measured comparison

All rates are cm/s. New closure values are recomputed independently in `../compare_pure_htpb.py` and checked against every case's recorded targets. New LowMach rates come from raw AMReX eta=0.5 surface-position fits, not initialized rates. Prior rates are the retained raw-output-derived original summary. Chen values come from the previously extracted PDF vector curves; graphical resolution is approximately 0.005 cm/s, not a confidence interval.

{baseline_table}

`parameter_comparison.pdf` and `.png` overlay the two material sets, analytic closures, simulations, and Chen's solid curves. `parameter_comparison.csv` includes both simulation-versus-paper and analytic-versus-paper errors, so numerical discrepancy is separated from the material/closure discrepancy. `analytic_curves.csv` retains both reconstructed closures and new Eq. 15 curves.

## Baseline numerical and heating checks

{refinement_table}

Baseline ell=delta/4, dy=ell/8 and delta=k/(ρ cp r_ref); the periodic transverse strip has 32 cells. Runs last six delta/r_ref with approximately forty output intervals. The fit uses the final 40% of physical time. Maximum absolute baseline late-window drift is {max(abs(s['late_speed_drift_percent']) for s in base):.3f}%. Surface-position and integrated-solid-volume speeds differ by at most {max(abs(100*(s['volume_r_cm_s']/s['r_cm_s']-1)) for s in new):.3f}% across the nine completed cases. The two refinements are incomplete; no convergence or refined-rate claim is made for this parameter set.

Baseline eta=0.5 surface-temperature offsets from the sharp-interface analytic targets span {min(s['error_surface_temperature_K'] for s in base):+.2f} to {max(s['error_surface_temperature_K'] for s in base):+.2f} K. Distributed kinetics across a finite thermal/interface transition are not identical to an Arrhenius law evaluated at a single interpolated temperature. Fit standard errors describe positional scatter and do not represent total numerical or material uncertainty.

{thermal_table}

The input-only fixed gas slab supplies integrated heat; gas-side conductive flux is measured at eta=10⁻⁶ using harmonic face conductivity. The artificial conducting layer retains k=100 W/(m K), cp=1000 J/(kg K), R=319.787 J/(kg K), P=100 MPa and no temperature advection. Gas storage and incident-flux errors quantify the actual delivered heating. The thermal residual uses source, latent heat from measured volume loss, lower-boundary loss and reconstructed thermal storage. It is a snapshot quadrature diagnostic; the separate steady-solid residual replaces solid storage by measured speed times the thermal-profile gradient. Sparse output can bias moving-interface storage even when solver steps are small. These diagnostics do not establish independence from surrogate conductivity or pressure.

Across the nine baselines, incident-flux errors span {min(s['incident_flux_error_percent'] for s in base):+.2f}% to {max(s['incident_flux_error_percent'] for s in base):+.2f}% of prescribed heat input, and gas storage spans {min(s['gas_storage_percent_input'] for s in base):.2f}% to {max(s['gas_storage_percent_input'] for s in base):.2f}%. These measured surrogate effects are part of the remaining speed bias; they should not be treated as material disagreement with the paper.

## Execution evidence and scope

Every selected baseline run has returncode=0, an input SHA-256 matching both execution receipt and case metadata, and the exact new material dictionary. Each solver log records AMReX finalization after reaching the requested duration. The last raw plot is within one final solver step of that time and at least 99.9% of the requested duration; some output writers retain the preceding step rather than the final step. `execution_audit.json` records commands, hashes, durations and independent-closure checks; `provenance.json` hashes the two source PDFs, parameter sets, original simulation summary and digitized curves. All nine completed receipts match the current executable hash: `{audit[0]['binary_sha256']}`. `stopped_cases.json` identifies the two excluded partial refinements. New extraction caches, summaries and figures are confined to this directory. The old run outputs, summaries and report are retained.

This is a planar condensed-phase heat-flux study using the existing executable, frozen gas chemistry, and one product gas. It does not reproduce heterogeneous DNS or validate reconciliation with Gross's premixed gas chemistry, composition-dependent gas reactions, flame heat feedback, or experimental pressure-dependent burning rates. Provisional cp/k and the changed constituent basis require separate physical validation before using this as a predictive combustion model.

## Reproduce postprocessing

From the repository root, with numpy/scipy/matplotlib/yt available:

```bash
python validation_homogeneous_chen2002/analysis/compare_study.py --tag _pure_htpb --output validation_homogeneous_chen2002/analysis/pure_htpb --refresh
python validation_homogeneous_chen2002/analysis/compare_pure_htpb.py
```

`simulation_summary.csv` / `.json` contain the nine completed baselines; `timeseries_*.csv` retain raw-field measurements. `numerical_diagnostics.pdf` / `.png` show the baseline key-case thermal checks. The original `write_report.py` must not be run on mixed material datasets. These commands only regenerate analysis of retained completed outputs; they do not resume the stopped simulations.
'''
    assert len({a["binary_sha256"] for a in audit}) == 1
    (OUT/"STUDY_REPORT.md").write_text(report)
    print(baseline_table)
    print(refinement_table)
    print(thermal_table)


if __name__ == "__main__":
    main()
