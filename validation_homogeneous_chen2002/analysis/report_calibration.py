#!/usr/bin/env python3
"""Audit stages of the rejected cp/Q calibration attempt.

Final reporting is disabled after the user clarified that A/E must be fitted.
--stage-only retains the ability to audit the completed rejected stages.
"""
import argparse
import csv
import json
import re
from pathlib import Path

from compare_pure_htpb import closure, digest, table
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from extract_runs import write_csv

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
OUT = HERE/"calibrated_htpb"


def audit(rows, parameters, arrhenius=False, fixed_parameters=None, relaxations=6,
          volume_density=False):
    original = json.loads((fixed_parameters or STUDY/"reference/pure_htpb.json").read_text())
    assert parameters["AP"] == original["AP"], "AP must remain untouched"
    fixed = ("density_kg_m3", "conductivity_W_m_K", "cp_J_kg_K", "heat_release_J_kg") if arrhenius else ("density_kg_m3", "conductivity_W_m_K", "A_m_s", "activation_temperature_K")
    for k in fixed:
        assert parameters["binder"][k] == original["binder"][k], k
    assert parameters["T0_K"] == original["T0_K"]
    for stage in parameters["arrhenius_calibration_history" if arrhenius else "calibration_history"]:
        assert stage["q_cal_cm2_s"] == [200.,1000.], "Held-out flux entered the fit"
        assert all(m["q_cal_cm2_s"] in (200.,1000.) for m in stage["measurements"])
    binary_hash = digest(STUDY.parent/"bin/lowmach-2d-clang++")
    result = []
    for s in rows:
        folder = STUDY/"runs"/s["name"]
        case = json.loads((folder/"case.json").read_text())
        run = json.loads((folder/"run.json").read_text())
        assert case["parameters"] == parameters, s["name"]
        if volume_density:
            assert case["density_mixing_rule"] == "volume_additive", s["name"]
        for k,v in {"width_ratio":.25,"cells_per_width":8,"dt_scale":.2,"relaxations":relaxations,"nx":32,
                    "gas_conductivity_W_m_K":100.,"gas_cp_J_kg_K":1000.,
                    "pressure_Pa":1.e8,"advect_temperature":False}.items():
            assert case[k] == v, (s["name"],k)
        assert run["returncode"] == 0 and run["binary_sha256"] == binary_hash, s["name"]
        assert run["input_sha256"] == case["input_sha256"] == digest(folder/"input"), s["name"]
        deck = (folder/"input").read_text()
        assert re.findall(r"^\s*chemistry\.model\.type\s*=\s*(\S+)\s*$", deck, re.M) == ["frozen"], s["name"]
        assert re.findall(r"^\s*mechanisms\.names\s*=\s*(.*?)\s*$", deck, re.M) == ["binder_regression"], s["name"]
        if arrhenius:
            assert not re.search(r"^\s*[^#\n]*temperature_cutoff\s*=",(folder/"input").read_text(),re.M), "Standalone study retains the default zero cutoff"
        log = (folder/"output/out.log").read_text()
        steps = re.findall(r"STEP\s+(\d+) ends\. TIME = ([\d.eE+\-]+) DT = ([\d.eE+\-]+)",log)
        assert steps and "finalized" in log, s["name"]
        step, end, dt = map(float,steps[-1])
        assert end >= case["duration_s"] and s["end_time_s"] >= .999*case["duration_s"]
        assert abs(end-s["end_time_s"]) <= 1.01*dt, s["name"]
        r,T,_ = closure(parameters,s["ap_volume_fraction"],s["q_cal_cm2_s"],
                        volume_density=volume_density)
        assert np.isclose(r,s["expected_r_eq14_19_cm_s"],rtol=1.e-10)
        result.append(dict(name=s["name"], returncode=0, input_sha256=run["input_sha256"],
            binary_sha256=binary_hash, requested_duration_s=case["duration_s"],
            final_plot_time_s=s["end_time_s"], final_log_time_s=end, final_step=int(step),
            last_dt_s=dt, parameters_match=True, untouched_AP_and_fixed_binder_properties=True,
            fitted_parameters="binder A and E/R only" if arrhenius else "binder cp and Q only (rejected)",
            temperature_cutoff_K=0.0 if arrhenius else None,
            gas_chemistry_frozen=True, phase_change_only_mechanism=True,
            density_mixing_rule="volume_additive" if volume_density else "arithmetic_mass_weighted",
            closure_recomputed=True, command=run["command"]))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",type=Path,required=True)
    parser.add_argument("--parameters",type=Path,required=True)
    parser.add_argument("--frozen-summary",type=Path,help="Seven frozen runs to combine with the two final-stage endpoints")
    parser.add_argument("--stage-only",action="store_true")
    args = parser.parse_args()
    if not args.stage_only:
        raise SystemExit("The user rejected the cp/Q fit; final reporting for that workflow is disabled.")
    rows = json.loads(args.summary.read_text())
    if args.frozen_summary:
        rows += json.loads(args.frozen_summary.read_text())
    parameters = json.loads(args.parameters.read_text())
    assert rows, "No completed measured cases"
    evidence = audit(rows,parameters)
    audit_dir = args.summary.parent if args.stage_only else OUT
    audit_dir.mkdir(parents=True,exist_ok=True)
    (audit_dir/"execution_audit.json").write_text(json.dumps(evidence,indent=2)+"\n")
    for s in rows:
        print(f"{s['name']}: {s['r_cm_s']:.7f} cm/s; target error {s['error_figure4_solid_percent']:+.4f}%; drift {s['late_speed_drift_percent']:+.4f}%; incident flux {s['incident_flux_error_percent']:+.3f}%")
    if args.stage_only:
        assert len(rows)==2 and {(s["q_cal_cm2_s"],s["ap_volume_fraction"]) for s in rows} == {(200.,0.),(1000.,0.)}
        return
    assert len(rows)==9 and {(s["q_cal_cm2_s"],s["ap_volume_fraction"]) for s in rows} == {(q,t) for q in (200.,500.,1000.) for t in (0.,.5,1.)}
    write_csv(OUT/"simulation_summary.csv",rows)
    (OUT/"simulation_summary.json").write_text(json.dumps(rows,indent=2)+"\n")
    provenance_paths=[args.parameters,args.summary,HERE/"figure4_curves.csv",HERE/"simulation_summary.json",
                      STUDY/"reference/pure_htpb.json"]
    if args.frozen_summary: provenance_paths.append(args.frozen_summary)
    (OUT/"provenance.json").write_text(json.dumps({str(p):digest(p) for p in provenance_paths},indent=2)+"\n")
    history = []
    for path in sorted(OUT.glob("stage*/simulation_summary.json")):
        stage_rows = json.loads(path.read_text())
        if not stage_rows: continue
        case = json.loads((STUDY/"runs"/stage_rows[0]["name"]/"case.json").read_text())
        stage_params = case["parameters"]
        audit(stage_rows,stage_params)
        for s in stage_rows:
            history.append(dict(stage=path.parent.name, name=s["name"], cp_J_kg_K=stage_params["binder"]["cp_J_kg_K"],
                Q_J_kg=stage_params["binder"]["heat_release_J_kg"], q_cal_cm2_s=s["q_cal_cm2_s"],
                r_cm_s=s["r_cm_s"], target_r_cm_s=s["figure4_solid_r_cm_s"],
                target_error_percent=s["error_figure4_solid_percent"], late_drift_percent=s["late_speed_drift_percent"]))
    write_csv(OUT/"stage_history.csv",history)
    old = {(s["q_cal_cm2_s"],s["ap_volume_fraction"]):s for s in json.loads((HERE/"simulation_summary.json").read_text()) if s["name"].endswith("_nx32")}
    stopped = {(s["q_cal_cm2_s"],s["ap_volume_fraction"]):s for s in json.loads((HERE/"pure_htpb/simulation_summary.json").read_text())}
    result = []
    for s in sorted(rows,key=lambda s:(s["q_cal_cm2_s"],s["ap_volume_fraction"])):
        q,t = s["q_cal_cm2_s"],s["ap_volume_fraction"]
        result.append(dict(q_cal_cm2_s=q,ap_volume_fraction=t,role="calibration endpoint" if t==0 and q!=500 else "held-out pure binder" if t==0 else "unfitted blend prediction" if t<1 else "fixed AP prediction",
            r_cm_s=s["r_cm_s"], chosen_closure_cm_s=s["expected_r_eq14_19_cm_s"], target_cm_s=s["figure4_solid_r_cm_s"],
            error_target_percent=s["error_figure4_solid_percent"], numerical_error_percent=s["error_eq14_19_percent"],
            prior_fitted_run_cm_s=old[q,t]["r_cm_s"], stopped_engineering_run_cm_s=stopped[q,t]["r_cm_s"]))
    write_csv(OUT/"comparison.csv",result)
    (OUT/"comparison.json").write_text(json.dumps(result,indent=2)+"\n")
    curves = list(csv.DictReader((HERE/"figure4_curves.csv").open()))
    fig,axs=plt.subplots(1,3,figsize=(11,3.8),layout="constrained")
    grid=np.linspace(0,1,201)
    for ax,q in zip(axs,(200,500,1000)):
        pts=[r for r in curves if r["kind"]=="eq15_solid" and float(r["q_cal_cm2_s"])==q]
        ax.plot([float(r["t"]) for r in pts],[float(r["r_cm_s"]) for r in pts],color="black")
        ax.plot(grid,[closure(parameters,t,q)[0] for t in grid],color="#0072B2",lw=1)
        for s in [s for s in result if s["q_cal_cm2_s"]==q]:
            t=s["ap_volume_fraction"]
            ax.plot(t,s["r_cm_s"],marker="s" if q==500 and t==0 else "o",color="#0072B2",ls="none",ms=5)
            ax.plot(t,s["prior_fitted_run_cm_s"],"x",color="#D55E00",ms=5)
        ax.set(title=f"q={q} cal cm⁻² s⁻¹",xlabel="AP volume fraction",xlim=(-.025,1.025))
    axs[0].set_ylabel("Regression speed (cm s⁻¹)")
    handles=[Line2D([],[],color=c,ls=ls,marker=m,label=l) for c,ls,m,l in [("black","-",None,"Chen Fig. 4"),("#0072B2","-",None,"Frozen analytic closure"),("#0072B2","none","o","New LowMach"),("#0072B2","none","s","Held-out pure HTPB"),("#D55E00","none","x","Original fitted study")]]
    fig.legend(handles=handles,loc="outside lower center",ncol=3,frameon=False,fontsize=8)
    for suffix in ("pdf","png"):fig.savefig(OUT/f"comparison.{suffix}",dpi=220)
    held=next(s for s in result if s["q_cal_cm2_s"]==500 and s["ap_volume_fraction"]==0)
    fit_max=max(abs(s["error_target_percent"]) for s in result if s["role"]=="calibration endpoint")
    b=parameters["binder"]
    summary_table=table(["q","AP volume fraction","Use","Measured [cm/s]","Chen [cm/s]","Error [%]"],[[f"{s['q_cal_cm2_s']:g}",f"{s['ap_volume_fraction']:g}",s["role"],f"{s['r_cm_s']:.5f}",f"{s['target_cm_s']:.5f}",f"{s['error_target_percent']:+.2f}"] for s in result])
    history_table=table(["Stage","cp [J/(kg K)]","Q [J/kg]","q","Measured [cm/s]","Target error [%]"],[[s["stage"],f"{s['cp_J_kg_K']:.3f}",f"{s['Q_J_kg']:.2f}",f"{s['q_cal_cm2_s']:g}",f"{s['r_cm_s']:.5f}",f"{s['target_error_percent']:+.3f}"] for s in history])
    report=f'''# Pure-HTPB-only calibration and frozen mixture predictions

The user selected Chen Figure 4's pure-binder model/DNS endpoint as the calibration target. This is a calibration to the supplied article, not to experimental pure-HTPB measurements. Only binder cp and decomposition heat Q were adjusted using q=200 and 1000 cal/(cm² s). The q=500 pure-binder point was held out, and no AP or mixed-composition rates entered the fit.

Final binder cp={b['cp_J_kg_K']:.6f} J/(kg K), Q={b['heat_release_J_kg']:.6f} J/kg ({b['heat_release_J_kg']/4184:.6f} cal/g). Fixed values are rho=920 kg/m³, k=0.13 W/(m K), A=10.36 m/s, E/R=7500 K and T0=300 K. All AP properties exactly equal `reference/pure_htpb.json`: rho=1950, cp=1297.90, k=0.4186, Q=−418400, A=948 and E/R=11000 in SI units. The original Gross −300 cal/g HTPB heat was released as a constraint; the AP −100 cal/g heat remains fixed.

The maximum absolute error at the two final fitting endpoints is {fit_max:.3f}%; the specified 1% endpoint acceptance criterion is {'met' if fit_max<1 else 'not met'}.

The held-out q=500 pure-binder run measures {held['r_cm_s']:.5f} cm/s against {held['target_cm_s']:.5f} cm/s, a {held['error_target_percent']:+.3f}% difference. This point assesses interpolation under the frozen fitted parameters. Mixed-composition discrepancies below are untouched predictions and cannot be counted as fitting success.

## Calibration history

{history_table}

The seed solves the two sharp-interface heat balances q/(rho r)=cp(Ts−T0)−Q, using Ts=(E/R)/ln(A/r) from each target rate. Subsequent bounded updates replace target analytic rates by target divided by the measured solver/analytic ratio. That correction is a fixed-point approximation to the numerical response, not an independently validated material-property inversion. With rho, A, E/R and T0 fixed, two distinct surface temperatures determine two coefficients cp/Q. Without those fixed assumptions the pure rates do not identify all material properties. Conductivity is absent from the steady surface balance; it is held fixed, not inferred.

## Frozen nine-case comparison

{summary_table}

`comparison.csv` also records each new analytic rate, its numerical error, the original article-fitted simulation and the stopped engineering-baseline result. The figure distinguishes the held-out pure-binder point and original simulations. Those original/stopped baseline comparators used dt ceiling scale=1, versus 0.2 here; cross-study rate changes include that timestep change as well as changed material properties. Fixed AP properties do not imply identical finite-step AP rates between studies. Chen vector digitization has graphical resolution about 0.005 cm/s; retained digits are not experimental precision. The chosen arithmetic mass-weighted density and geometric Arrhenius mixing remain unchanged. Calibrating pure HTPB does not resolve the mixture-density ambiguity described in `../REFERENCE_NOTES.md`.

## Verification and limits

All nine final cases have successful execution receipts, exact input hashes matching case metadata, the same current executable hash, matching frozen parameter dictionaries, and completed solver logs. The two completed final-stage fitting endpoints are reused, combined with the seven new frozen-parameter runs; no endpoint rerun is required. All use ell/delta=0.25, eight cells per interface width, dt ceiling scale=0.2, 32 transverse cells and six thermal relaxation times. Last raw plots reach at least 99.9% of requested duration and lie within one final solver step. Rates are independently fitted from eta=0.5 AMReX surface positions over the final 40% of each run, checked against volume loss and late-window drift. No initialized rate is used as a measurement. `execution_audit.json` retains case-level evidence, and `simulation_summary.json` combines precisely those nine measurements.

Maximum absolute late-window speed drift is {max(abs(s['late_speed_drift_percent']) for s in rows):.3f}%; surface-position and integrated-volume rates differ by at most {max(abs(100*(s['volume_r_cm_s']/s['r_cm_s']-1)) for s in rows):.3f}%. Numerical errors against the final analytic closure range from {min(s['error_eq14_19_percent'] for s in rows):+.2f}% to {max(s['error_eq14_19_percent'] for s in rows):+.2f}%. Measured incident-flux errors range from {min(s['incident_flux_error_percent'] for s in rows):+.2f}% to {max(s['incident_flux_error_percent'] for s in rows):+.2f}% of input heat, and gas storage from {min(s['gas_storage_percent_input'] for s in rows):.2f}% to {max(s['gas_storage_percent_input'] for s in rows):.2f}%. Matching calibration targets is distinct from eliminating those numerical/surrogate effects.

This fit is conditional on the executed interface width, spatial/time discretization, six-thermal-relaxation duration, and the fixed input-only gas-slab heating surrogate. Numerical-rate bias is partly absorbed into cp/Q; the final parameters should not be presented as measured pure-material constants or mesh-independent calibration. Gas chemistry is frozen, temperature advection disabled, and no Gross premixed flame/experimental-pressure validation is established. Original fitted artifacts and the stopped-study snapshot remain intact.

## Reproduction

The retained frozen parameter file is `{args.parameters}` (SHA-256 `{digest(args.parameters)}`). Its `calibration_history` records every fit input, endpoint target and measured ratio. `calibrate_pure_htpb.py` creates each stage without overwriting an existing parameter file; `run_study.py` generates retained tagged cases. Postprocess each completed tag with `compare_study.py --tag TAG --output validation_homogeneous_chen2002/analysis/calibrated_htpb/STAGE`, then run:

```bash
python validation_homogeneous_chen2002/analysis/report_calibration.py --summary {args.summary} --parameters {args.parameters}{' --frozen-summary '+str(args.frozen_summary) if args.frozen_summary else ''}
```

Use `--stage-only` for the two fitting endpoints before a frozen sweep exists. This report performs analysis only and launches no simulations.
'''
    (OUT/"STUDY_REPORT.md").write_text(report)
    print(summary_table)


if __name__=="__main__":main()
