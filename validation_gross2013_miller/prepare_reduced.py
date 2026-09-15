#!/usr/bin/env python3
"""Prepare 12 flame-initialized cases and four zero-step preview cases."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
os.environ.setdefault('MPLCONFIGDIR','/tmp/alamo-gross-matplotlib')
import numpy as np
import yt
from prepare_runs import generate, HERE, ROOT
yt.set_log_level(50)
SPECIES=('AP_gas','HTPB_gas','Mono','Premixed','Primary','Final')

def source_profile():
    source=HERE/'runs/M03_p68.0000_seed101_pilot_dx1_w4_v5_restart/output/19013cell'
    ds=yt.load(str(source));cg=ds.covering_grid(0,ds.domain_left_edge,ds.domain_dimensions)
    y=float(ds.domain_left_edge[1])+(np.arange(ds.domain_dimensions[1])+.5)*float(ds.domain_width[1])/ds.domain_dimensions[1]
    mean=lambda f:cg['boxlib',f].d[:,:,0].mean(axis=0)
    eta=mean('rigid_eta');i=np.flatnonzero((eta[:-1]>=.5)&(eta[1:]<.5))[0]
    surface=y[i]+(.5-eta[i])*(y[i+1]-y[i])/(eta[i+1]-eta[i]);z=y-surface
    T=mean('temperature');u=np.maximum(mean('velocityy'),0)
    rho=np.column_stack([mean('component_density_'+s) for s in SPECIES])
    total=rho.sum(axis=1);Y=np.maximum(rho,0)/np.maximum(total[:,None],1e-100)
    Y[(total<1e-8)|(z<-8e-6)]=[.8737,.1263,0,0,0,0]
    Y/=Y.sum(axis=1)[:,None]
    Ts=float(np.interp(0,z,T));Tf=float(np.mean(T[z>100e-6]))
    # Tail values are held constant outside this sampled flame/preheat layer.
    keep=(z>=-80e-6)&(z<=200e-6)
    z0=z[keep];v0=np.column_stack((T,u,Y))[keep]
    zz=np.sort(np.r_[z0,0]);vv=np.column_stack([np.interp(zz,z0,v0[:,j]) for j in range(8)])
    # Retain shared knots so the six piecewise-linear fractions sum to one.
    selected={0,len(zz)-1,int(np.searchsorted(zz,0))};tolerance=np.array([1.,.003]+[.001]*6)
    while True:
        idx=sorted(selected);reconstructed=np.column_stack([np.interp(zz,zz[idx],vv[idx,j]) for j in range(8)])
        error=np.max(np.abs(reconstructed-vv)/tolerance,axis=1);k=int(np.argmax(error))
        if error[k]<=1:break
        selected.add(k)
    idx=sorted(selected);zz=zz[idx];vv=vv[idx]
    profile=dict(source_plotfile=str(source.relative_to(ROOT)),source_time_s=float(ds.current_time),
        source_header_sha256=hashlib.sha256((source/'Header').read_bytes()).hexdigest(),
        source_pressure_atm=68.,source_surface_y_m=float(surface),surface_temperature_K=Ts,
        product_temperature_K=Tf,coordinate_m=zz.tolist(),values=vv.tolist(),
        columns=['temperature_K','velocity_y_m_s']+list(SPECIES),
        note='Coupled planar seed extracted 0.9 ms after external heating ended. It is an initial guess, not a measured steady result for any new pack.')
    out=HERE/'initial_conditions';out.mkdir(exist_ok=True)
    (out/'flame_seed_source.json').write_text(json.dumps(profile,indent=2)+'\n')
    return profile

def expr(knots,values,variable='y'):
    """Balanced conditional tree for exact piecewise-linear interpolation."""
    def branch(lo,hi):
        if hi==lo+1:
            a,b=values[lo],values[hi]
            if abs(b-a)<1e-16:return f'({a:.17g})'
            return f'({a:.17g}+({(b-a)/(knots[hi]-knots[lo]):.17g})*({variable}-({knots[lo]:.17g})))'
        mid=(lo+hi)//2
        return f'if({variable}<({knots[mid]:.17g}),{branch(lo,mid)},{branch(mid,hi)})'
    return f'if({variable}<({knots[0]:.17g}),({values[0]:.17g}),if({variable}>({knots[-1]:.17g}),({values[-1]:.17g}),{branch(0,len(knots)-1)}))'

def set_flame(case,profile):
    m=json.loads((case/'case.json').read_text());P=m['pressure_Pa'];ratio=m['pressure_atm']/68.
    # Dimensional seed scaling from first-order gas rates ~P^1.6/P^1.7:
    # flame thickness ~P^-0.825, gas speed ~P^-0.175. No Miller rate fitting.
    exponent=.825;scale=ratio**(-exponent)
    z=np.array(profile['coordinate_m']);v=np.array(profile['values']);y=z*scale
    Ts0=profile['surface_temperature_K'];Tf=profile['product_temperature_K']
    activation=9526.52460524
    Ts=activation/(activation/Ts0-exponent*np.log(ratio))
    T=np.where(z<=0,300+(v[:,0]-300)*(Ts-300)/(Ts0-300),Ts+(v[:,0]-Ts0)*(Tf-Ts)/(Tf-Ts0))
    T=np.maximum(T,300);u=v[:,1]*ratio**(exponent-1)*T/v[:,0];Y=v[:,2:]
    w=m['interface_width_m'];temperature=expr(y,T);velocity=expr(y,u)
    deck={}
    for line in (case/'input').read_text().splitlines():
        if '=' in line and not line.lstrip().startswith('#'):
            k,value=line.split('=',1);deck[k.strip()]=value.strip()
    for key in list(deck):
        if key.startswith(('temperature.ic.','velocity.ic.','Final.density.ic.','heat_source.ic.')):del deck[key]
    deck['temperature.ic.type']='expression';deck['temperature.ic.expression.region0']=f'"{temperature}"'
    deck['velocity.ic.type']='expression';deck['velocity.ic.expression.region0']='"0"'
    deck['velocity.ic.expression.region1']=f'"{velocity}"'
    # Seed gas density is EOS-consistent with local T and diffuse gas occupancy.
    R=8.31446261815324/.026
    for j,s in enumerate(SPECIES):
        fraction=expr(y,Y[:,j]);prefix=s+'.density.ic.'
        deck[prefix+'type']='expression'
        deck[prefix+'expression.region0']=f'"temp=({temperature}); frac=({fraction}); (0.5+0.5*tanh(2*y/{w:.17g}))*{P:.17g}/({R:.17g}*temp)*max(0,frac)"'
    # First segment is deliberately short. Root evaluates correlated block
    # statistics and extends only when duration/recession is insufficient.
    initial_segment=.0005
    deck['stop_time']=f'{initial_segment:.17g}';deck['amr.plot_dt']='2.5e-5'
    deck['run_control.stop_file']=str(case.relative_to(ROOT)/'STOP')
    deck['run_control.poll_interval']='100'
    text='# Reduced study: developed flame at t=0; no external heating.\n'
    text+='\n'.join(f'{k} = {value}' for k,value in deck.items())+'\n'
    (case/'input').write_text(text)
    m.update(initial_condition='developed_planar_flame',ignition_off_s=0.,requested_stop_s=initial_segment,
        plot_dt_s=2.5e-5,input_sha256=hashlib.sha256(text.encode()).hexdigest(),
        resolve_M03_20um=True,initial_flame=dict(source=profile['source_plotfile'],source_time_s=profile['source_time_s'],
        source_header_sha256=profile['source_header_sha256'],pressure_scaling_exponent=exponent,
        length_scale=scale,surface_temperature_seed_K=Ts,product_temperature_seed_K=Tf,
        note='Free initial field; surface temperature is not prescribed after t=0. Same planar gas seed across formulations; resolved flames develop during relaxation.'),
        duration_policy=dict(initial_segment_s=initial_segment,maximum_reference_duration_s=.0015+2*max(float(d) for d in m['packing']['bins'])*1e-6/(m['reference_rate_cm_s']*.01),
        criterion='Discard startup; require four correlated-time blocks, <5% drift, and relative 95% temporal mean uncertainty <=5%, with particle-scale sampling. Extend if unmet; never accept merely because the segment ends.',
        packing_uncertainty='One realization: inter-pack variance is unavailable.'))
    (case/'case.json').write_text(json.dumps(m,indent=2)+'\n')
    np.savez_compressed(HERE/'initial_conditions'/f"{m['case']}_profile.npz",y_m=y,temperature_K=T,velocity_y_m_s=u,mass_fractions=Y)
    return m

def main():
    profile=source_profile();references=list(csv.DictReader((HERE/'reference/figure10_digitized.csv').open()))
    production=[];previews=[];cases=[]
    for name in ('M03','M17','M21','M24'):
        rows=[r for r in references if r['formulation']==name and r['kind']=='miller_diamond']
        for position in (0,2,6):
            row=rows[position]
            path=generate(name,101,float(row['pressure_atm']),float(row['rate_cm_s']),tag='_warm',ignition_off=0.,pack_dir=HERE/'packs_reduced',resolve_m03=True)
            case=ROOT/path;m=set_flame(case,profile);production.append(path);cases.append(m)
            if position==2:
                preview=case.with_name(case.name+'_initial');preview.mkdir(exist_ok=True)
                if (preview/'run.json').exists():raise FileExistsError(preview)
                text=(case/'input').read_text().replace(str(case.relative_to(ROOT)),str(preview.relative_to(ROOT))).replace('max_step = 2147483647','max_step = 0')
                (preview/'input').write_text(text)
                pm=dict(m,case=preview.name,pilot=True,initialization_only=True,production_case=path,input_sha256=hashlib.sha256(text.encode()).hexdigest())
                (preview/'case.json').write_text(json.dumps(pm,indent=2)+'\n');previews.append(str(preview.relative_to(ROOT)))
    old=HERE/'production_cases.txt';archive=HERE/'production_cases_84_archived.txt'
    if not archive.exists():shutil.copy2(old,archive)
    old.write_text('\n'.join(production)+'\n')
    (HERE/'preview_cases.txt').write_text('\n'.join(previews)+'\n')
    (HERE/'initial_conditions/reduced_cases.json').write_text(json.dumps(cases,indent=2)+'\n')
    print(json.dumps(dict(production_cases=len(production),preview_cases=previews,source_knots=len(profile['coordinate_m'])),indent=2))

if __name__=='__main__':main()
