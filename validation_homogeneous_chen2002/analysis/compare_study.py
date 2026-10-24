#!/usr/bin/env python3
"""Compare completed standalone LowMach cases with independently read Fig. 4.

Run after the sweep: python analysis/compare_study.py [--study STUDY_DIRECTORY].
This reads case.json, run.json, and raw plotfiles. Analytic values in case.json
are reported as model targets; measured speeds always come from plotfiles.
"""
import argparse
import csv
import json
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/chen-analysis-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from extract_runs import read_case_plots, summarize, write_csv

HERE=Path(__file__).resolve().parent
COLORS={200:"#0072B2",500:"#D55E00",1000:"#009E73"}


def run_kind(case):
    if case["width_ratio"] < .249:
        return "thinner interface + smaller step", "D"
    if case.get("dt_scale",1) < .99:
        return "smaller time step", "^"
    return "baseline", "x"


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--study",type=Path,default=HERE.parent)
    p.add_argument("--late-fraction",type=float,default=.4)
    p.add_argument("--refresh",action="store_true",help="Ignore cached timeseries")
    p.add_argument("--tag",default="",help="Only process case names containing this substring")
    p.add_argument("--output",type=Path,default=HERE,help="Directory for summaries, plots, and caches")
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    curves=list(csv.DictReader((HERE/"figure4_curves.csv").open()))
    markers=list(csv.DictReader((HERE/"figure4_markers.csv").open()))
    figure={}
    for q in COLORS:
        rs=[r for r in curves if r["kind"]=="eq15_solid" and float(r["q_cal_cm2_s"])==q]
        figure[q]=np.array([[float(r["t"]),float(r["r_cm_s"])] for r in rs])
    summaries,skipped=[],[]
    for meta_path in sorted((args.study/"runs").glob("*/case.json")):
        case=json.loads(meta_path.read_text())
        if args.tag not in case["name"]:continue
        # Pilot overrides can differ from generator defaults: the executed
        # input is authoritative for whether thermal advection was enabled.
        deck=(meta_path.parent/"input").read_text()
        matches=re.findall(r"^\s*advect_temperature\s*=\s*(\d+)",deck,re.M)
        if matches:case["advect_temperature"]=bool(int(matches[-1]))
        status_path=meta_path.parent/"run.json"
        if not status_path.exists() or json.loads(status_path.read_text()).get("returncode") != 0:
            skipped.append(dict(name=case["name"],reason="No successful run.json"));continue
        plots=sorted(p.parent for p in (meta_path.parent/"output").glob("*cell/Header"))
        cache=args.output/f"timeseries_{case['name']}.csv"
        try:
            newest=max((p/"Header").stat().st_mtime for p in plots)
            if cache.exists() and cache.stat().st_mtime>=newest and not args.refresh:
                rows=[{k:(v if k=="plotfile" else float(v) if v else float("nan")) for k,v in r.items()} for r in csv.DictReader(cache.open())]
            else:
                rows=read_case_plots(plots,case)
                rows=sorted({r["time_s"]:r for r in rows}.values(),key=lambda r:r["time_s"])
                write_csv(cache,rows)
            if rows[-1]["time_s"] < .999*case["duration_s"]:
                raise ValueError("Run stopped before requested duration")
            s=summarize(rows,args.late_fraction)
            ref=case["reference"]
            q,t=ref["heat_flux_cal_cm2_s"],ref["ap_volume_fraction"]
            rf=float(np.interp(t,figure[q][:,0],figure[q][:,1]))
            r19=100*ref["regression_speed_m_s"]
            r15=100*ref["eq15_regression_speed_m_s"]
            s={k:v for k,v in case.items() if isinstance(v,(str,int,float,bool))}|s
            s.update(q_cal_cm2_s=q,ap_volume_fraction=t,
                expected_r_eq14_19_cm_s=r19,expected_r_eq15_cm_s=r15,
                figure4_solid_r_cm_s=rf,
                expected_surface_temperature_K=ref["surface_temperature_K"],
                error_eq14_19_percent=100*(s["r_cm_s"]/r19-1),
                error_eq15_percent=100*(s["r_cm_s"]/r15-1),
                error_figure4_solid_percent=100*(s["r_cm_s"]/rf-1),
                error_surface_temperature_K=s["surface_temperature_K"]-ref["surface_temperature_K"],
                thermal_length_m=ref["thermal_length_m"],
                thermal_cells=ref["thermal_length_m"]/case["dx_m"],
                final_solid_thermal_lengths=s["final_solid_depth_m"]/ref["thermal_length_m"])
            summaries.append(s)
            print(case["name"],f"{s['r_cm_s']:.5g} cm/s, {s['error_eq14_19_percent']:+.2f}% vs closure",flush=True)
        except Exception as exc:
            skipped.append(dict(name=case["name"],reason=str(exc)))
    write_csv(args.output/"simulation_summary.csv",summaries)
    (args.output/"simulation_summary.json").write_text(json.dumps(summaries,indent=2)+"\n")
    (args.output/"skipped_cases.json").write_text(json.dumps(skipped,indent=2)+"\n")
    if not summaries: raise SystemExit("No completed cases to compare")
    plt.rcParams.update({"font.size":9,"pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
    fig,axs=plt.subplots(2,1,figsize=(6.5,6.8),sharex=True,layout="constrained",gridspec_kw={"height_ratios":[2,1]})
    for q,color in COLORS.items():
        axs[0].plot(*figure[q].T,color=color,label=f"Fig. 4 Eq. 15, q={q}")
        for kind,marker in [("dns_2d_circle","o"),("dns_3d_asterisk","*")]:
            rs=[r for r in markers if r["kind"]==kind and float(r["q_cal_cm2_s"])==q]
            axs[0].plot([float(r["t"]) for r in rs],[float(r["r_cm_s"]) for r in rs],linestyle="none",marker=marker,ms=4,mfc="none",mec=color,mew=.7,alpha=.65)
        ss=[s for s in summaries if s["q_cal_cm2_s"]==q]
        for s in ss:
            _,mark=run_kind(s)
            axs[0].scatter(s["ap_volume_fraction"],s["r_cm_s"],c=color,marker=mark,s=30,zorder=4)
            axs[1].scatter(s["ap_volume_fraction"],s["error_eq14_19_percent"],c=color,marker=mark,s=30)
    axs[0].set(ylabel="Regression speed (cm s⁻¹)",title="Homogeneous LowMach regression versus Chen et al. (2002)")
    axs[0].legend(frameon=False,fontsize=8)
    axs[0].text(.015,.025,"○ paper 2D DNS   ✶ paper 3D DNS\nPrescribed flux q in cal cm⁻² s⁻¹",transform=axs[0].transAxes,fontsize=8)
    axs[1].axhline(0,color=".3",lw=.8)
    axs[1].set(xlabel="AP volume fraction, t",ylabel="Error vs chosen\nhomogeneous closure (%)",xlim=(-.025,1.025))
    kinds=sorted(set(run_kind(s) for s in summaries),key=lambda k:{"baseline":0,"smaller time step":1,"thinner interface + smaller step":2}[k[0]])
    axs[1].legend(handles=[Line2D([],[],color=".2",marker=m,linestyle="none",label=k) for k,m in kinds],frameon=False,fontsize=8,loc="best")
    for suffix in ("pdf","png"): fig.savefig(args.output/f"comparison.{suffix}",dpi=240)
    key=sorted([s for s in summaries if s["q_cal_cm2_s"]==500 and s["ap_volume_fraction"]==.5],
               key=lambda s: (s["width_ratio"]<.249,s.get("dt_scale",1)<.99))
    fig,axs=plt.subplots(1,2,figsize=(8.2,3.7),layout="constrained")
    positions=np.arange(len(key))
    labels=["Baseline" if run_kind(s)[0]=="baseline" else "Smaller\ntime step" if s["width_ratio"]>=.249 else "Thinner interface\n+ smaller step" for s in key]
    axs[0].plot(positions,[s["error_eq14_19_percent"] for s in key],"o-",color=COLORS[500],lw=1)
    for s,x in zip(key,positions):
        axs[0].annotate(f"{s['r_cm_s']:.4f} cm/s",(x,s["error_eq14_19_percent"]),xytext=(0,8),textcoords="offset points",ha="center",fontsize=8)
    for field,label,marker in [("incident_flux_error_percent","Incident flux error","o"),
                               ("gas_storage_percent_input","Gas thermal storage","s"),
                               ("thermal_balance_residual_percent_input","Snapshot thermal residual","^"),
                               ("steady_profile_balance_residual_percent_input","Steady solid profile residual","D")]:
        axs[1].plot(positions,[s[field] for s in key],marker=marker,lw=1,label=label)
    axs[0].set(ylabel="Speed error vs closure (%)",title="q = 500 cal cm⁻² s⁻¹, t = 0.5")
    axs[1].set(ylabel="Fraction of prescribed heat input (%)",title="Heating surrogate diagnostics")
    axs[1].legend(frameon=False,fontsize=8)
    for ax in axs:
        ax.axhline(0,color=".4",lw=.6)
        ax.set_xticks(positions,labels)
        ax.margins(x=.2,y=.3)
    for suffix in ("pdf","png"): fig.savefig(args.output/f"numerical_diagnostics.{suffix}",dpi=240)


if __name__ == "__main__":
    main()
