#!/usr/bin/env python3
"""Infer only identifiable heat-balance combinations from Fig. 4 endpoints.

Not independent physical parameter recovery: this calibrates to the same figure
later used for comparison. No constituent rho, cp, Q, or T0 is identified alone.
The corner perturbation range is graphical sensitivity, NOT confidence limits.
"""
import csv
import itertools
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def fit(q,r,A,E):
    temperature=E/np.log(A/r)
    alpha,intercept=np.polyfit(temperature,q/r,1)
    return alpha,-intercept/alpha,temperature


def main():
    rows=list(csv.DictReader((HERE/"figure4_curves.csv").open()))
    provenance=json.loads((HERE/"figure4_provenance.json").read_text())
    dr=provenance["conservative_graphical_resolution"]["r_cm_s"]
    results={}
    table=[]
    for name,t,A,E in [("binder",0,1036,7500),("AP",1,94800,11000)]:
        selected=sorted([r for r in rows if r["kind"]=="eq15_solid" and abs(float(r["t"])-t)<1e-8],key=lambda r:float(r["q_cal_cm2_s"]))
        q=np.array([float(r["q_cal_cm2_s"]) for r in selected])
        r=np.array([float(r["r_cm_s"]) for r in selected])
        alpha,beta,T=fit(q,r,A,E)
        perturbed=np.array([fit(q,r+dr*np.array(signs),A,E)[:2] for signs in itertools.product((-1,1),repeat=3)])
        results[name]=dict(A_cm_s=A,activation_temperature_K=E,
            inferred_rho_cp_cal_cm3_K=alpha,inferred_rho_cp_J_m3_K=alpha*4.184e6,
            inferred_T0_plus_Q_over_cp_K=beta,
            graphical_corner_sensitivity_rho_cp_cal_cm3_K=[float(perturbed[:,0].min()),float(perturbed[:,0].max())],
            graphical_corner_sensitivity_offset_K=[float(perturbed[:,1].min()),float(perturbed[:,1].max())],
            heat_balance_residual_cal_cm2_s=(alpha*r*(T-beta)-q).tolist())
        for qi,ri,Ti in zip(q,r,T):
            table.append(dict(material=name,q_cal_cm2_s=qi,r_cm_s=ri,arrhenius_implied_Ts_K=Ti))
    output=dict(source="Only user-supplied Chen et al. (2002) article: Fig.4 vector endpoints; Fig.5/6 captions give E/R and A.",
        method="At pure endpoints Eqs.12,14 give Ts=(E/R)/ln(A/r), q/r=(rho*c)*Ts-(rho*c)*(T0+Q/c). Least-squares straight line through three endpoint heat balances identifies rho*c and T0+Q/c only.",
        limitations="These are Figure-4-calibrated combinations, not tabulated physical parameters or independent validation. Individual rho,c,Q,T0 and conductivity cannot be inferred uniquely. Fig.4 gives no DNS uncertainty. One-full-stroke endpoint perturbations are an illustrative graphical sensitivity, not confidence bounds. Constituent blend-density/mass-weighting predictions require additional conventions, which must be stated and tested separately.",
        endpoint_perturbation_cm_s=dr,materials=results)
    (HERE/"article_inferred_combinations.json").write_text(json.dumps(output,indent=2)+"\n")
    with (HERE/"figure4_pure_endpoints.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=table[0].keys());w.writeheader();w.writerows(table)
    audit=[]
    b,a=results["binder"],results["AP"]
    for row in csv.DictReader((HERE/"figure4_markers.csv").open()):
        if row["kind"]!="eq19_square":continue
        t=float(row["t"]);r=float(row["r_cm_s"]);q=float(row["q_cal_cm2_s"])
        A=1036**(1-t)*94800**t;E=7500*(1-t)+11000*t
        T=E/np.log(A/r)
        alpha=(1-t)*b["inferred_rho_cp_cal_cm3_K"]+t*a["inferred_rho_cp_cal_cm3_K"]
        beta=((1-t)*b["inferred_rho_cp_cal_cm3_K"]*b["inferred_T0_plus_Q_over_cp_K"]+
              t*a["inferred_rho_cp_cal_cm3_K"]*a["inferred_T0_plus_Q_over_cp_K"])/alpha
        audit.append(dict(q_cal_cm2_s=q,t=t,r_square_cm_s=r,
            eq19_implied_Ts_K=T,volume_weighted_rho_cp_cal_cm3_K=alpha,
            capacity_weighted_offset_K=beta,
            required_capacity_multiplier=q/(r*alpha*(T-beta))))
    with (HERE/"figure4_blend_capacity_audit.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=audit[0].keys());w.writeheader();w.writerows(audit)
    print(json.dumps(results,indent=2))


if __name__=="__main__":main()
