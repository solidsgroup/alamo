#!/usr/bin/env python3
"""Aggregate production runs only, separating missing/failed/unstable results."""
import csv
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/alamo-gross-matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

def main():
    paths=[ROOT/p for p in (HERE/'production_cases.txt').read_text().splitlines()]
    reference=list(csv.DictReader((HERE/'reference/figure10_digitized.csv').open()))
    queued={}
    for state_path in sorted((HERE/'queues').glob('*.state.json')):
        state=json.loads(state_path.read_text())
        for case_path in state.get('pending',[]):
            queued[ROOT/case_path]='queued' if state['status']=='running' else 'queue_stopped'
    records=[];groups={}
    for case in paths:
        m=json.loads((case/'case.json').read_text())
        s=json.loads((case/'analysis/steady_assessment.json').read_text()) if (case/'analysis/steady_assessment.json').exists() else {}
        r=json.loads((case/'run.json').read_text()) if (case/'run.json').exists() else {}
        timing=json.loads((case/'analysis/steady_stop.json').read_text()) if (case/'analysis/steady_stop.json').exists() else {}
        execution=r.get('status',queued.get(case,'not_run'))
        status=execution if execution in ('failed','launcher_error') else s.get('status',execution)
        accepted=(r.get('returncode')==0 and r.get('amrex_finalized') and s.get('accepted') and s.get('status')=='accepted_steady')
        candidate=s.get('candidate',{})
        record=dict(case=case.name,formulation=m['formulation'],pressure_atm=m['pressure_atm'],seed=m['seed'],
            status=status,execution_status=execution,assessment_status=s.get('status'),
            accepted=bool(accepted),last_time_s=s.get('last_time_s'),requested_stop_s=m['requested_stop_s'],
            rate_cm_s=candidate.get('block_mean_cm_s') if accepted else None,
            temporal_ci95_halfwidth_cm_s=candidate.get('temporal_ci95_halfwidth_cm_s') if accepted else None,
            half_window_drift_percent=candidate.get('half_window_drift_percent'),reference_rate_cm_s=m['reference_rate_cm_s'],
            first_passing_simulation_s=timing.get('first_passing_simulation_s'),
            time_to_steady_simulation_s=timing.get('time_to_steady_simulation_s'),
            time_to_steady_wall_s=timing.get('time_to_steady_wall_s'),
            actual_stop_simulation_s=timing.get('final_run',{}).get('final_simulation_time_s'),
            actual_run_wall_s=r.get('wall_seconds'),
            statistical_stop_requested=timing.get('stop_requested_while_running',False))
        records.append(record);groups.setdefault((m['formulation'],m['pressure_atm']),[]).append(record)
    out=HERE/'analysis';out.mkdir(exist_ok=True)
    with (out/'production_status.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    timing_fields=['case','formulation','pressure_atm','status','accepted','rate_cm_s',
                   'temporal_ci95_halfwidth_cm_s','first_passing_simulation_s',
                   'time_to_steady_simulation_s','time_to_steady_wall_s',
                   'actual_stop_simulation_s','actual_run_wall_s','statistical_stop_requested']
    with (out/'time_to_steady.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=timing_fields,extrasaction='ignore');w.writeheader();w.writerows(records)
    means=[]
    for (name,p),rows in sorted(groups.items()):
        rates=[r['rate_cm_s'] for r in rows if r['accepted']]
        complete=len(rates)==len(rows)
        mean=float(np.mean(rates)) if complete else None
        std=float(np.std(rates,ddof=1)) if complete and len(rates)>1 else None
        temporal_ci=rows[0]['temporal_ci95_halfwidth_cm_s'] if complete and len(rows)==1 else None
        gross=[r for r in reference if r['formulation']==name and r['kind']=='gross_new_solid']
        xp=np.array([float(r['pressure_atm']) for r in gross]);yp=np.array([float(r['rate_cm_s']) for r in gross])
        gross_rate=float(np.exp(np.interp(np.log(p),np.log(xp),np.log(yp)))) if xp[0]<=p<=xp[-1] else None
        means.append(dict(formulation=name,pressure_atm=p,accepted_instances=len(rates),required_instances=len(rows),
            mean_rate_cm_s=mean,packing_sd_cm_s=std,temporal_ci95_halfwidth_cm_s=temporal_ci,miller_rate_cm_s=rows[0]['reference_rate_cm_s'],
            error_vs_miller_percent=100*(mean/rows[0]['reference_rate_cm_s']-1) if complete else None,
            gross_solid_rate_cm_s=gross_rate,
            error_vs_gross_solid_percent=100*(mean/gross_rate-1) if complete and gross_rate is not None else None,
            ensemble_note='One resolved pack: temporal uncertainty only; inter-pack variability is unavailable.'))
    with (out/'figure10_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=means[0]);w.writeheader();w.writerows(means)
    accepted=sum(r['accepted'] for r in records)
    status=dict(requested_cases=len(records),accepted_cases=accepted,
        running_cases=sum(r['execution_status']=='running' for r in records),
        queued_cases=sum(r['execution_status']=='queued' for r in records),
        failed_cases=sum(r['execution_status'] in ('failed','launcher_error') for r in records),
        completed_simulations=sum(r['execution_status']=='completed' for r in records),
        complete_pressure_groups=sum(m['accepted_instances']==m['required_instances'] for m in means),
        publication_ready=False,
        note='Temporal acceptance includes particle-scale recession. Root must also audit mesh/interface-width and shortened-domain checks before accepting the final figure.')
    (out/'sweep_status.json').write_text(json.dumps(status,indent=2)+'\n')
    # Do not publish an empty "simulation comparison" if only references exist.
    if not any(m['accepted_instances']==m['required_instances'] for m in means):
        print(json.dumps(status,indent=2));return
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for ax,name in zip(axes.flat,('M03','M17','M21','M24')):
        for kind,style,label in [('gross_new_solid','-','Gross new, digitized'),('miller_diamond','D','Miller, digitized')]:
            rr=[r for r in reference if r['kind']==kind and r['formulation']==name]
            ax.plot([float(r['pressure_atm']) for r in rr],[float(r['rate_cm_s']) for r in rr],style,color='.4' if style=='-' else 'k',ms=4,label=label)
        rr=[r for r in means if r['formulation']==name and r['accepted_instances']==r['required_instances']]
        for i,r in enumerate(rr):
            ax.errorbar(r['pressure_atm'],r['mean_rate_cm_s'],yerr=r['temporal_ci95_halfwidth_cm_s'],fmt='o',color='#0072B2',ms=5,capsize=3,label='One pack: mean ± temporal 95% CI' if i==0 else None)
            ax.annotate(f"{r['error_vs_miller_percent']:+.1f}%",(r['pressure_atm'],r['mean_rate_cm_s']),xytext=(5,6),textcoords='offset points',fontsize=8,color='#0072B2')
        ax.set(xscale='log',yscale='log',xlabel='Pressure (atm)',ylabel='Regression speed (cm/s)',title=name,xlim=(5,250));ax.grid(alpha=.2,which='both')
    axes[0,0].legend(fontsize=8)
    fig.suptitle(f'Gross Figure 10 comparison: {status["complete_pressure_groups"]}/{len(means)} pressure groups\nCalibrated solids; fixed gas coefficients. Text: error versus Miller data.')
    for ext in ('png','pdf'):fig.savefig(out/f'figure10_comparison.{ext}',dpi=220)
    plt.close(fig);print(json.dumps(status,indent=2))

if __name__=='__main__':main()
