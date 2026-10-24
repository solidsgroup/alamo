#!/usr/bin/env python3
"""Root-model postprocessing of actual AMReX output (including failed pilots)."""
import argparse
import csv
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/alamo-gross-matplotlib')
import numpy as np
import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rocfire_diagnostics import heat_release
yt.set_log_level(50)
HERE=Path(__file__).resolve().parent

def snapshot(path):
    ds=yt.load(str(path));ad=ds.all_data()
    area=ad['index','dx'].d*ad['index','dy'].d
    eta=ad['boxlib','rigid_eta'].d;T=ad['boxlib','temperature'].d
    stored_heat=('boxlib','qdot') in ds.field_list
    qdot=(ad['boxlib','qdot'].d if stored_heat else
          heat_release(ad,json.loads((path.parent.parent/'case.json').read_text())))
    width=float(ds.domain_width[0]);lo_y=float(ds.domain_left_edge[1])
    # yt masks covered coarse cells: no AMR double-counting. Raw coordinates
    # are the solver's SI values; ignore yt's default cgs labels for plotfiles.
    height=lo_y+np.sum(eta*area)/width
    levels=ds.index.max_level;dims=ds.domain_dimensions*ds.refine_by**levels;dims[2:]=1
    dy=float(ds.domain_width[1])/int(dims[1])
    # Surface interpolation needs only the strip containing the diffuse front,
    # not a finest-level allocation over the entire millimetre-deep solid bed.
    transition=(eta>.05)&(eta<.95)
    left=ds.domain_left_edge.copy();strip_lo=0
    if np.any(transition):
        y=ad['index','y'].d[transition];cell_dy=ad['index','dy'].d[transition]
        strip_lo=max(0,int(np.floor((np.min(y-2*cell_dy)-lo_y)/dy)))
        strip_hi=min(int(dims[1]),int(np.ceil((np.max(y+2*cell_dy)-lo_y)/dy)))
        left[1]=lo_y+strip_lo*dy;dims[1]=strip_hi-strip_lo
    cg=ds.covering_grid(levels,left,dims)
    e=cg['boxlib','rigid_eta'].d[:,:,0];temp=cg['boxlib','temperature'].d[:,:,0]
    temps=[];positions=[];multi=0
    for ec,tc in zip(e,temp):
        cross=np.flatnonzero((ec[:-1]>=.5)!=(ec[1:]>=.5))
        multi+=len(cross)>1
        for k in cross:
            w=(.5-ec[k])/(ec[k+1]-ec[k]);temps.append(tc[k]+w*(tc[k+1]-tc[k]))
            positions.append(lo_y+(strip_lo+k+.5+w)*dy)
    P=ad['boxlib','pressure'].d
    row=dict(plotfile=path.name,time_s=float(ds.current_time),mean_solid_height_m=height,
        surface_temperature_K=float(np.mean(temps)) if temps else None,
        surface_temperature_min_K=float(np.min(temps)) if temps else None,
        surface_temperature_max_K=float(np.max(temps)) if temps else None,
        interface_min_y_m=float(np.min(positions)) if positions else None,
        interface_max_y_m=float(np.max(positions)) if positions else None,
        multivalued_column_fraction=multi/len(e),temperature_min_K=float(T.min()),temperature_max_K=float(T.max()),
        gas_heat_release_W_m2=float(np.sum(qdot*area)/width),pressure_min_Pa=float(P.min()),pressure_max_Pa=float(P.max()))
    for s in ['AP_solid','HTPB_solid']:
        row[s+'_mass_kg_m2']=float(np.sum(ad['boxlib','component_density_'+s].d*area)/width)
    uy=ad['boxlib','velocityy'].d;ux=ad['boxlib','velocityx'].d
    row['velocity_max_m_s']=float(np.max(np.hypot(ux,uy)))
    # Cell-centered top-boundary flux estimate, used as a chemistry diagnostic.
    # These are not claimed to be the solver's discrete face fluxes.
    top=np.isclose(ad['index','y'].d+.5*ad['index','dy'].d,
                   float(ds.domain_right_edge[1]),rtol=0,atol=1e-12)
    dx=ad['index','dx'].d
    total_flux=0.
    for s in ['AP_gas','HTPB_gas','Mono','Premixed','Primary','Final']:
        rho=ad['boxlib','component_density_'+s].d
        row[s+'_mass_kg_m2']=float(np.sum(rho*area)/width)
        flux=float(np.sum(rho[top]*np.maximum(uy[top],0)*dx[top])/width)
        row[s+'_outlet_flux_estimate_kg_m2_s']=flux;total_flux+=flux
    row['gas_outlet_flux_estimate_kg_m2_s']=total_flux
    row['gas_heat_method']='solver_diagnostic' if stored_heat else 'reconstructed_Rocfire_from_physical_fields'
    row['snapshot_schema']=3
    return row

def fit(t,h):
    if len(t)<3:return None
    return float(-np.polyfit(t-t[0],h,1)[0]*100)

def history(case, visited=None, completed_only=False):
    visited=set() if visited is None else visited
    if case in visited:raise ValueError(f'Restart ancestry cycle: {case}')
    visited.add(case)
    meta=json.loads((case/'case.json').read_text());out=case/'analysis';out.mkdir(exist_ok=True)
    cache_path=out/'snapshots.json';cache=json.loads(cache_path.read_text()) if cache_path.exists() else {}
    plots=sorted((case/'output').glob('*cell'))
    if completed_only:
        # The integrator appends this entry only after writing every level and
        # closing Checkpoint. Header alone can exist during an incomplete write.
        visit=case/'output/celloutput.visit'
        entries=visit.read_text().splitlines(keepends=True) if visit.exists() else []
        names={Path(line.strip()).parent.name for line in entries
               if line.endswith('\n') and line.strip().endswith('/Header')}
        plots=[p for p in plots if p.name in names and (p/'Checkpoint').exists()]
    for p in plots:
        if not (p/'Header').exists():continue
        if p.name in cache and cache[p.name].get('snapshot_schema')==3:continue
        try:cache[p.name]=snapshot(p)
        except (OSError,ValueError) as exc:print(f'{p}: skipped potentially incomplete write: {exc}',flush=True)
    cache_path.write_text(json.dumps(cache,indent=2)+'\n')
    rows=sorted((dict(r,source_case=case.name) for r in cache.values()),key=lambda r:r['time_s'])
    for row in rows:row['gas_heat_diagnostic_source']=case.name+'/'+row['plotfile']
    if rows and meta.get('restart_from'):
        checkpoint=Path(meta['restart_from'])
        if not checkpoint.is_absolute():checkpoint=HERE.parent/checkpoint
        parent=checkpoint.resolve().parent.parent
        if parent.parent != HERE/'runs':raise ValueError(f'Restart parent outside study: {parent}')
        # Exclude superseded parent evolution after the selected checkpoint.
        # Legacy restart uses float precision for time; merge near-equal times.
        cutoff=rows[0]['time_s']-1.e-10
        previous=history(parent,visited,completed_only=completed_only)
        match=[r for r in previous if r['source_case']==parent.name and r['plotfile']==checkpoint.name]
        if match:
            # A restart's first plot is written before its reference pressure
            # is synchronized. Use the checkpoint's initialized qdot diagnostic
            # for this identical physical instant, avoiding a false heat dip.
            rows[0]['gas_heat_release_W_m2']=match[-1]['gas_heat_release_W_m2']
            rows[0]['gas_heat_method']=match[-1]['gas_heat_method']
            rows[0]['gas_heat_diagnostic_source']=match[-1]['gas_heat_diagnostic_source']
        rows=[r for r in previous if r['time_s']<cutoff]+rows
    return rows

def analyze(case):
    meta=json.loads((case/'case.json').read_text());out=case/'analysis';out.mkdir(exist_ok=True)
    rows=history(case)
    if not rows:return {'case':case.name,'status':'no_readable_output'}
    with (out/'history.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    t=np.array([r['time_s'] for r in rows]);h=np.array([r['mean_solid_height_m'] for r in rows])
    start=max(meta['ignition_off_s']+.0005,.5*t[-1]);mask=t>=start
    rate=fit(t[mask],h[mask]);mid=.5*(start+t[-1])
    r1=fit(t[(t>=start)&(t<=mid)],h[(t>=start)&(t<=mid)]);r2=fit(t[t>=mid],h[t>=mid])
    drift=None if rate is None or r1 is None or r2 is None else 100*(r2-r1)/max(abs(rate),1.e-12)
    receipt=json.loads((case/'run.json').read_text()) if (case/'run.json').exists() else {}
    reached=bool(t[-1]>=meta['requested_stop_s']*(1-1.e-6))
    status='in_progress'
    returncode=receipt.get('returncode',receipt.get('exit_code'))
    if returncode is not None and returncode!=0:status='failed'
    elif reached:
        status='stationary_candidate' if drift is not None and abs(drift)<5 else 'not_stationary'
    if t[-1]>meta['ignition_off_s']+.0005 and rows[-1]['temperature_max_K']<800:
        status='cooling_or_extinguished'
    fit_recession=float(h[mask][0]-h[mask][-1]) if mask.sum()>1 else None
    diameters=[float(d) for d in meta.get('packing',{}).get('bins',{})]
    largest_diameter=max(diameters,default=0)*1e-6
    summary=dict(case=case.name,status=status,pressure_atm=meta['pressure_atm'],formulation=meta['formulation'],
        seed=meta['seed'],pilot=meta['pilot'],last_time_s=float(t[-1]),requested_stop_s=meta['requested_stop_s'],
        completed_duration=reached,fit_start_s=start,fit_rate_cm_s=rate,fit_half_drift_percent=drift,
        reference_rate_cm_s=meta['reference_rate_cm_s'],
        comparison_error_percent=(100*(rate/meta['reference_rate_cm_s']-1) if rate is not None and not meta['pilot'] else None),
        late_surface_temperature_K=float(np.mean([r['surface_temperature_K'] for r in rows if r['time_s']>=start and r['surface_temperature_K'] is not None])) if np.any(mask) else None,
        max_temperature_K=max(r['temperature_max_K'] for r in rows),
        final_max_temperature_K=rows[-1]['temperature_max_K'],
        recession_m=float(h[0]-h[-1]),fit_samples=int(mask.sum()),
        fit_recession_m=fit_recession,
        fit_recession_largest_particle_diameters=fit_recession/largest_diameter if fit_recession is not None and largest_diameter>0 else None,
        history_start_s=float(t[0]),
        note='Window drift is a stationarity diagnostic, not a statistical confidence interval. Pilot rates are not validation data.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    fig,ax=plt.subplots(2,2,figsize=(10,6),layout='constrained')
    ax[0,0].plot(t*1e3,(h[0]-h)*1e6);ax[0,0].set(ylabel='Mean recession (µm)')
    ax[0,1].plot(t*1e3,[r['surface_temperature_K'] for r in rows],label='Surface')
    ax[0,1].plot(t*1e3,[r['temperature_max_K'] for r in rows],label='Domain maximum')
    ax[0,1].legend();ax[0,1].set(ylabel='Temperature (K)')
    if len(t)>2:ax[1,0].plot(.5*(t[:-1]+t[1:])*1e3,-np.diff(h)/np.diff(t)*100)
    ax[1,0].set(ylabel='Interval regression speed (cm/s)')
    ax[1,1].plot(t*1e3,[r['gas_heat_release_W_m2']/1e6 for r in rows]);ax[1,1].set(ylabel='Integrated gas heat release (MW/m²)')
    for a in ax.flat:a.axvline(meta['ignition_off_s']*1e3,color='.5',ls='--');a.set(xlabel='Time (ms)');a.grid(alpha=.2)
    fig.suptitle(f'{case.name}\n{status}');fig.savefig(out/'history.png',dpi=160);plt.close(fig)
    return summary

def main():
    parser=argparse.ArgumentParser();parser.add_argument('cases',nargs='+',type=Path);args=parser.parse_args()
    summaries=[analyze(c.resolve()) for c in args.cases]
    print(json.dumps(summaries,indent=2))

if __name__=='__main__':main()
