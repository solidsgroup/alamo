#!/usr/bin/env python3
"""Root-authored temporal stationarity and correlated-block rate assessment."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.stats import t as student_t
from analyze_runs import history

def save_assessment(case,answer):
    out=case/'analysis';out.mkdir(exist_ok=True)
    temporary=out/'steady_assessment.json.tmp'
    temporary.write_text(json.dumps(answer,indent=2)+'\n')
    temporary.replace(out/'steady_assessment.json')
    return answer

def assessment(case, completed_only=False):
    meta=json.loads((case/'case.json').read_text());rows=history(case,completed_only=completed_only)
    answer=dict(case=case.name,status='insufficient_history',accepted=False,
        uncertainty_scope='Temporal variation within this pack; no inter-pack uncertainty is available.',
        note='Block confidence intervals are approximate. Stationarity and particle-scale sampling are required separately.')
    if len(rows)<4:return save_assessment(case,answer)
    t=np.array([r['time_s'] for r in rows]);h=np.array([r['mean_solid_height_m'] for r in rows])
    largest=max(float(d) for d in meta['packing']['bins'])*1e-6
    final=rows[-1];answer.update(last_time_s=float(t[-1]),recession_m=float(h[0]-h[-1]),
        final_max_temperature_K=final['temperature_max_K'],largest_particle_m=largest,
        status='not_yet_steady',candidates=[])
    if t[-1]>.001 and final['temperature_max_K']<800:
        answer.update(status='cooling_or_extinguished');return save_assessment(case,answer)
    # The initial guess is excluded even though it already contains a flame.
    startup=.00025
    for fraction in (0.,.2,.4):
        start=startup+fraction*max(0,t[-1]-startup)
        mask=t>=start
        if mask.sum()<20:continue
        ts=t[mask];hs=h[mask]
        interval=float(np.max(np.diff(ts)))
        n=int(np.floor((ts[-1]-ts[0])/interval))
        if n<19:continue
        grid=np.linspace(ts[0],ts[-1],n+1)
        heights=np.interp(grid,t,h);rates=-np.diff(heights)/np.diff(grid)*100
        rate=float(-np.polyfit(ts-ts[0],hs,1)[0]*100)
        centered=rates-rates.mean();variance=float(np.dot(centered,centered)/len(centered))
        tau=1.;acf=[]
        if variance>max(rate*rate,1e-20)*1e-20:
            for lag in range(1,len(rates)//2):
                rho=float(np.dot(centered[:-lag],centered[lag:])/(len(rates)-lag)/variance)
                acf.append(rho)
                if rho<=0:break
                tau+=2*rho
        block_size=max(4,int(np.ceil(2*tau)));count=len(rates)//block_size
        candidate=dict(start_s=float(ts[0]),end_s=float(ts[-1]),rate_cm_s=rate,
            recession_m=float(hs[0]-hs[-1]),recession_particle_diameters=float((hs[0]-hs[-1])/largest),
            correlation_time_s=float(tau*interval),block_size_samples=block_size,block_count=count,
            interval_samples=len(rates),criteria_pass=False)
        if count>=4 and rate>1e-8:
            chunks=np.array_split(rates,count);block_rates=np.array([a.mean() for a in chunks])
            mean=float(block_rates.mean());sd=float(block_rates.std(ddof=1))
            halfwidth=float(student_t.ppf(.975,count-1)*sd/np.sqrt(count))
            trend=abs(float(np.polyfit(np.arange(count),block_rates,1)[0]))*(count-1)/mean
            half=abs(float(block_rates[:count//2].mean()-block_rates[count//2:].mean()))/mean
            q=np.interp(grid,t,[r['gas_heat_release_W_m2'] for r in rows])
            qblocks=np.array([a.mean() for a in np.array_split(q,count)])
            qdrift=abs(float(qblocks[:count//2].mean()-qblocks[count//2:].mean()))/max(abs(qblocks.mean()),1e-20)
            ci=halfwidth/mean
            passed=bool(candidate['recession_particle_diameters']>=1 and ci<=.05 and trend<=.05 and half<=.05 and qdrift<=.05 and qblocks.min()>0)
            candidate.update(block_rates_cm_s=block_rates.tolist(),block_mean_cm_s=mean,
                temporal_ci95_halfwidth_cm_s=halfwidth,temporal_ci95_relative_percent=100*ci,
                block_trend_percent=100*trend,half_window_drift_percent=100*half,
                gas_heat_half_drift_percent=100*qdrift,criteria_pass=passed,
                block_duration_s=float((ts[-1]-ts[0])/count))
        answer['candidates'].append(candidate)
        if candidate['criteria_pass']:
            answer.update(status='steady_candidate',candidate=candidate);break
    out=case/'analysis';out.mkdir(exist_ok=True)
    prior_path=out/'steady_assessment.json'
    prior=json.loads(prior_path.read_text()) if prior_path.exists() else {}
    # Confirm on a later snapshot separated by at least one estimated block.
    if answer['status']=='steady_candidate':
        first=prior.get('first_passing_time_s',float(t[-1])) if prior.get('status') in ('steady_candidate','accepted_steady') else float(t[-1])
        answer['first_passing_time_s']=first
        if t[-1]-first>=answer['candidate']['block_duration_s']:
            answer.update(status='accepted_steady',accepted=True)
    return save_assessment(case,answer)

def main():
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);args=p.parse_args()
    print(json.dumps(assessment(args.case.resolve()),indent=2))

if __name__=='__main__':main()
