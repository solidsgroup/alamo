#!/usr/bin/env python3
"""Independently audit A/E-only calibration stages and report a frozen sweep.

Use --stage-only for two endpoint runs. Final mode combines the two reused
final-stage endpoints with --frozen-summary containing seven unfitted runs.
No simulations or parameter edits are performed.
"""
import argparse
import csv
import json
from pathlib import Path

from compare_pure_htpb import closure, digest, table
from report_calibration import audit
from extract_runs import write_csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE=Path(__file__).resolve().parent
STUDY=HERE.parent
OUT=HERE/"arrhenius_htpb"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",type=Path,required=True)
    parser.add_argument("--parameters",type=Path,required=True)
    parser.add_argument("--frozen-summary",type=Path)
    parser.add_argument("--stage-only",action="store_true")
    args=parser.parse_args()
    rows=json.loads(args.summary.read_text())
    if args.frozen_summary: rows+=json.loads(args.frozen_summary.read_text())
    parameters=json.loads(args.parameters.read_text())
    evidence=audit(rows,parameters,arrhenius=True)
    folder=args.summary.parent if args.stage_only else OUT
    folder.mkdir(parents=True,exist_ok=True)
    (folder/"execution_audit.json").write_text(json.dumps(evidence,indent=2)+"\n")
    for s in rows:
        print(f"{s['name']}: {s['r_cm_s']:.7f} cm/s; target error {s['error_figure4_solid_percent']:+.4f}%; drift {s['late_speed_drift_percent']:+.4f}%; Ts {s['surface_temperature_K']:.3f} K; incident flux {s['incident_flux_error_percent']:+.3f}%",flush=True)
    if args.stage_only:
        assert len(rows)==2 and {(s["q_cal_cm2_s"],s["ap_volume_fraction"]) for s in rows}=={(200.,0.),(1000.,0.)}
        return
    assert len(rows)==9 and {(s["q_cal_cm2_s"],s["ap_volume_fraction"]) for s in rows}=={(q,t) for q in (200.,500.,1000.) for t in (0.,.5,1.)}
    write_csv(OUT/"simulation_summary.csv",rows)
    (OUT/"simulation_summary.json").write_text(json.dumps(rows,indent=2)+"\n")
    history=[]
    for path in sorted(OUT.glob("stage*/simulation_summary.json")):
        stage_rows=json.loads(path.read_text())
        if not stage_rows:continue
        params=json.loads((STUDY/"runs"/stage_rows[0]["name"]/"case.json").read_text())["parameters"]
        audit(stage_rows,params,arrhenius=True)
        for s in stage_rows:
            history.append(dict(stage=path.parent.name,name=s["name"],A_m_s=params["binder"]["A_m_s"],
                activation_temperature_K=params["binder"]["activation_temperature_K"],q_cal_cm2_s=s["q_cal_cm2_s"],
                r_cm_s=s["r_cm_s"],target_error_percent=s["error_figure4_solid_percent"],
                late_drift_percent=s["late_speed_drift_percent"],surface_temperature_K=s["surface_temperature_K"]))
    write_csv(OUT/"stage_history.csv",history)
    prior={(s["q_cal_cm2_s"],s["ap_volume_fraction"]):s for s in json.loads((HERE/"simulation_summary.json").read_text()) if s["name"].endswith("_nx32")}
    comparison=[]
    for s in sorted(rows,key=lambda s:(s["q_cal_cm2_s"],s["ap_volume_fraction"])):
        q,t=s["q_cal_cm2_s"],s["ap_volume_fraction"]
        comparison.append(dict(q_cal_cm2_s=q,ap_volume_fraction=t,
            role="fitted endpoint" if t==0 and q!=500 else "held-out pure binder" if t==0 else "unfitted blend" if t<1 else "fixed AP prediction",
            r_cm_s=s["r_cm_s"],analytic_r_cm_s=s["expected_r_eq14_19_cm_s"],chen_r_cm_s=s["figure4_solid_r_cm_s"],
            target_error_percent=s["error_figure4_solid_percent"],numerical_error_percent=s["error_eq14_19_percent"],
            prior_fitted_r_cm_s=prior[q,t]["r_cm_s"],surface_temperature_K=s["surface_temperature_K"]))
    write_csv(OUT/"comparison.csv",comparison)
    (OUT/"comparison.json").write_text(json.dumps(comparison,indent=2)+"\n")
    paths=[args.parameters,args.summary,HERE/"figure4_curves.csv",HERE/"simulation_summary.json",STUDY/"reference/pure_htpb.json"]
    if args.frozen_summary:paths.append(args.frozen_summary)
    (OUT/"provenance.json").write_text(json.dumps({str(p):digest(p) for p in paths},indent=2)+"\n")
    curves=list(csv.DictReader((HERE/"figure4_curves.csv").open()))
    fig,axs=plt.subplots(1,3,figsize=(11,3.9),layout="constrained")
    grid=np.linspace(0,1,201)
    analytic=[]
    for ax,q in zip(axs,(200,500,1000)):
        pts=[p for p in curves if p["kind"]=="eq15_solid" and float(p["q_cal_cm2_s"])==q]
        ax.plot([float(p["t"]) for p in pts],[float(p["r_cm_s"]) for p in pts],color="black",lw=1.2)
        values=[closure(parameters,t,q) for t in grid]
        ax.plot(grid,[v[0] for v in values],color="#0072B2",lw=1.2)
        for t,v in zip(grid,values):analytic.append(dict(q_cal_cm2_s=q,ap_volume_fraction=t,eq14_19_r_cm_s=v[0],surface_temperature_K=v[1],eq15_r_cm_s=v[2]))
        for s in [s for s in comparison if s["q_cal_cm2_s"]==q]:
            t=s["ap_volume_fraction"]
            ax.plot(t,s["r_cm_s"],marker="s" if q==500 and t==0 else "o",color="#0072B2",ls="none",ms=5)
            ax.plot(t,s["prior_fitted_r_cm_s"],"x",color="#D55E00",ms=5)
        ax.set(title=f"q={q} cal cm⁻² s⁻¹",xlabel="AP volume fraction, t",xlim=(-.025,1.025))
    axs[0].set_ylabel("Regression speed (cm s⁻¹)")
    fig.legend(handles=[Line2D([],[],color=c,ls=ls,marker=m,label=l) for c,ls,m,l in [("black","-",None,"Chen Fig. 4"),("#0072B2","-",None,"A/E-only analytic closure"),("#0072B2","none","o","New LowMach"),("#0072B2","none","s","Held-out pure binder"),("#D55E00","none","x","Original fitted study")]],loc="outside lower center",ncol=3,frameon=False,fontsize=8)
    for suffix in ("pdf","png"):fig.savefig(OUT/f"comparison.{suffix}",dpi=220)
    write_csv(OUT/"analytic_curves.csv",analytic)
    b=parameters["binder"]
    held=next(s for s in comparison if s["role"]=="held-out pure binder")
    fit_error=max(abs(s["target_error_percent"]) for s in comparison if s["role"]=="fitted endpoint")
    history_table=table(["Stage","A [m/s]","E/R [K]","q","Measured [cm/s]","Target error [%]"],[[s["stage"],f"{s['A_m_s']:.6f}",f"{s['activation_temperature_K']:.3f}",f"{s['q_cal_cm2_s']:g}",f"{s['r_cm_s']:.5f}",f"{s['target_error_percent']:+.3f}"] for s in history])
    result_table=table(["q","t","Role","LowMach [cm/s]","Chen [cm/s]","Error [%]","Surface T [K]"],[[f"{s['q_cal_cm2_s']:g}",f"{s['ap_volume_fraction']:g}",s["role"],f"{s['r_cm_s']:.5f}",f"{s['chen_r_cm_s']:.5f}",f"{s['target_error_percent']:+.2f}",f"{s['surface_temperature_K']:.2f}"] for s in comparison])
    report=f'''# Arrhenius-only pure-HTPB calibration and frozen predictions

Only binder Arrhenius A and E/R were fitted to the user-confirmed Chen Figure 4 pure-binder endpoints at q=200 and 1000 cal/(cm² s). Final A={b['A_m_s']:.9f} m/s and E/R={b['activation_temperature_K']:.6f} K. The maximum absolute fitted-endpoint error is {fit_error:.3f}%; the 1% numerical target criterion is {'met' if fit_error<1 else 'not met'}. The q=500 pure-binder point was held out and measures {held['r_cm_s']:.5f} cm/s against {held['chen_r_cm_s']:.5f}, a {held['target_error_percent']:+.3f}% discrepancy. No held-out, AP or mixed-composition result was used to refit parameters.

This is an effective calibration to article model/DNS curves, not experimental pure-HTPB measurements or validated physical kinetics. The earlier cp/Q calibration was rejected by the user and is excluded from these results.

## Fixed quantities and calibration history

Binder rho=920 kg/m³, cp=2418.29 J/(kg K), k=0.13 W/(m K), Q=−1255200 J/kg (−300 cal/g), and T0=300 K remain unchanged. Every AP property exactly equals `reference/pure_htpb.json`: rho=1950, cp=1297.90, k=0.4186, Q=−418400, A=948 and E/R=11000 in SI units. The audit verifies these equalities for every stage and final run. No thermal parameter is fitted.

{history_table}

For each fitting target, Ts=T0+[q/(rho*r)+Q]/cp follows from the fixed thermal balance. Two distinct positive Ts values determine log(A) and E/R from log(r)=log(A)−(E/R)/Ts. The initial two-endpoint solution gives A=0.452931711 m/s and E/R=1008.342770 K. Any numerical correction adjusts only those two quantities using measured solver/analytic rate ratios from the two fitting fluxes. This is a bounded fixed-point correction, conditional on the executed discretization and gas heating surrogate.

The original endpoint targets imply Ts=249.968 K at q=200 and 378.203 K at q=1000. At low flux the surface is below the 300 K deep-solid temperature; endothermic regression can receive sensible heat from cooling incoming solid. This is mathematically possible in the implemented balance, but differs from a hot-surface pyrolysis interpretation. If Ts≥300 K were required, the q=200 target would violate the fixed-Q energy bound regardless of A/E. `FEASIBILITY.md` details that condition and the initial analytic held-out prediction (+1.0244%). Numerical fits do not establish that this low-temperature kinetic regime represents real HTPB chemistry.

The standalone decks omit `temperature_cutoff` and retain the implementation's default 0 K. In contrast, the original `input.lm.ap_htpb` explicitly sets HTPB pyrolysis cutoff to 360 K; `PhaseChange.H` makes the kinetic factor zero at or below the cutoff. The q=200 and q=500 targets require surface temperatures below 360 K with the fixed thermal data. These fitted kinetics therefore cannot simply replace the original production values while preserving that cutoff and reproduce the same prescribed-flux targets. No cutoff, thermal property, gas chemistry or production input is changed here. The calibration applies to the actual standalone setup, not directly to the original production configuration.

## Frozen nine-case comparison

{result_table}

The two completed final-stage fitting endpoints are reused with seven new frozen-parameter runs. Pure AP and blend entries are unfitted predictions. `comparison.csv` separately records errors against the selected analytic closure, the digitized Chen curve, and original article-fitted rates. `comparison.pdf` / `.png` overlay these results. Original baseline comparators used dt ceiling scale=1, versus 0.2 here, so cross-study differences include timestep effects. Fixed AP properties do not imply identical finite-step AP rates between studies. The arithmetic mass-weighted density and geometric Arrhenius mixing rules are unchanged; pure-binder calibration does not settle the density-convention ambiguity in `../REFERENCE_NOTES.md`.

## Numerical evidence and limitations

All nine runs have returncode=0, matching executed-input and metadata hashes, the current executable hash, exact frozen parameter dictionaries and completed solver logs. The final raw plots reach at least 99.9% of requested duration and lie within one final solver step. Rates come from raw AMReX eta=0.5 surface-position fits over the final 40% of each run, independently checked against volume loss. All use ell/delta=0.25, eight cells per interface width, 32 transverse cells, dt ceiling scale=0.2 and six thermal relaxation times.

Maximum absolute late-window drift is {max(abs(s['late_speed_drift_percent']) for s in rows):.3f}%; volume-loss and surface-position rates differ by at most {max(abs(100*(s['volume_r_cm_s']/s['r_cm_s']-1)) for s in rows):.3f}%. Numerical errors against the final analytic closure span {min(s['error_eq14_19_percent'] for s in rows):+.2f}% to {max(s['error_eq14_19_percent'] for s in rows):+.2f}%. Measured incident-flux errors span {min(s['incident_flux_error_percent'] for s in rows):+.2f}% to {max(s['incident_flux_error_percent'] for s in rows):+.2f}% of prescribed heat input; gas storage spans {min(s['gas_storage_percent_input'] for s in rows):.2f}% to {max(s['gas_storage_percent_input'] for s in rows):.2f}%. These surrogate and discretization effects are partly absorbed by the fitted kinetics; matching endpoints is not a convergence proof.

Heating is a fixed input-only gas slab with gas k=100 W/(m K), cp=1000 J/(kg K), P=100 MPa and temperature advection disabled. Chemistry is frozen. This does not validate Gross premixed gas chemistry or experimental pressure-dependent burning. Chen curve resolution is approximately 0.005 cm/s; tiny numerical fitting residuals do not imply equivalent physical precision. Original and stopped/rejected study artifacts remain intact.

## Reproduction

Frozen parameters: `{args.parameters}`; SHA-256 `{digest(args.parameters)}`. Its `arrhenius_calibration_history` records every target, fitting input and measured ratio. Each stage is independently extracted using `compare_study.py --tag TAG --output validation_homogeneous_chen2002/analysis/arrhenius_htpb/STAGE`. To audit two fitting endpoints, add `--stage-only` to this script. Final reporting combines retained endpoints and frozen predictions:

```bash
python validation_homogeneous_chen2002/analysis/report_arrhenius.py --summary {args.summary} --parameters {args.parameters}{' --frozen-summary '+str(args.frozen_summary) if args.frozen_summary else ''}
```

`stage_history.csv`, `simulation_summary.csv` / `.json`, `execution_audit.json`, `analytic_curves.csv` and `provenance.json` preserve the measurements and reproducibility evidence. Analysis commands never launch simulations or change parameters.
'''
    (OUT/"STUDY_REPORT.md").write_text(report)
    print(result_table)


if __name__=="__main__":main()
