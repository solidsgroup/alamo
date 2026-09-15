#!/usr/bin/env python3
"""Compact domains and material-aware diffusion-flame initial guesses.

This is an initialization model, not a second combustion solver. Surface
streams mix through a constant-D convective diffusion approximation. Progress
fractions from the saved developed flame are mapped onto the local conserved
AP/binder mixture, with all six species kept nonnegative. No data in Figure 10
are used to tune this initial field.
"""
import hashlib
import json
import math
from pathlib import Path
import shutil
import numpy as np
from prepare_reduced import expr, HERE, ROOT, SPECIES

BETA=7.974
GAMMA=3.19
R=8.31446261815324/.026

def deck_read(path):
    return dict((a.strip(),b.strip()) for line in path.read_text().splitlines()
                if '=' in line and not line.lstrip().startswith('#')
                for a,b in [line.split('=',1)])

def save_case(case,deck,meta):
    case.mkdir(exist_ok=True)
    if (case/'run.json').exists():raise FileExistsError(case)
    text='# Compact domain; approximate local diffusion-flame seed; no external heating.\n'
    text+='\n'.join(f'{k} = {v}' for k,v in deck.items())+'\n'
    (case/'input').write_text(text)
    meta=dict(meta,case=case.name,input_sha256=hashlib.sha256(text.encode()).hexdigest())
    (case/'case.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta

def surface_modes(meta):
    width=meta['packing']['width_mm']*1e-3
    # Oversample the disk intersections; the Fourier filter smooths on the
    # resolved diffuse material interface, not on a whole particle diameter.
    n=8192;x=np.arange(n)*width/n;mask=np.zeros(n)
    disk_path=ROOT/meta['packing_file'] if meta.get('packing_file') else HERE/'packs_reduced'/f"{meta['formulation']}_seed101.xyzr"
    disks=np.loadtxt(disk_path)*1e-3
    for xx,yy,_,radius in disks:
        if abs(yy)>=radius:continue
        half=np.sqrt(radius*radius-yy*yy)
        distance=(x-xx+.5*width)%width-.5*width
        mask=np.maximum(mask,abs(distance)<half)
    modes=np.fft.rfft(mask)/n
    count=min(192,int(np.ceil(width/(4*meta['dx_m']))))
    return modes[:count+1],width,float(mask.mean())

def fourier(modes,width,sigma,decay):
    terms=[f'({modes[0].real:.17g})']
    for n,c in enumerate(modes[1:],1):
        k=2*np.pi*n/width;filter_=np.exp(-.5*(k*sigma)**2)
        if abs(c)*filter_<2e-5:continue
        trig=f'({2*c.real*filter_:.17g}*cos({k:.17g}*x)+{-2*c.imag*filter_:.17g}*sin({k:.17g}*x))'
        terms.append(f'{trig}*exp(-({decay(k):.17g})*max(y,0))')
    return '('+'+'.join(terms)+')'

def remap_progress(profile):
    v=np.array(profile['values']);Y=v[:,2:]
    e2=Y[:,4]/(BETA+1);e3=Y[:,5]/(GAMMA+1)
    e0=Y[:,2]+GAMMA*e3;e1=Y[:,3]+e3
    Z=Y[:,0]+e0+BETA*e2
    divide=lambda a,b:np.divide(a,b,out=np.zeros_like(a),where=b>1e-12)
    g2=np.clip(divide(e2,np.minimum(Z/BETA,1-Z)),0,1)
    g0=np.clip(divide(e0,Z-BETA*e2),0,1)
    g1=np.clip(divide(e1,1-Z-e2),0,1)
    g3=np.clip(divide(e3,np.minimum(e0/GAMMA,e1)),0,1)
    return np.column_stack((g0,g1,g2,g3))

def adiabatic_expr(end):
    # The final extents are piecewise linear in mixture fraction. Compress
    # their enthalpy into a small table to stay within the native parser's
    # 16-slot stack, without changing the initialization model.
    Z=np.linspace(0,1,10001)
    e2=np.minimum(Z/BETA,1-Z)*end[2]
    e0=(Z-BETA*e2)*end[0];e1=(1-Z-e2)*end[1]
    e3=np.minimum(e0/GAMMA,e1)*end[3]
    T=300+(430*e0+447*e1+8365*e2+2127*e3-66-34*Z)/.3
    selected={0,len(Z)-1}
    while True:
        ix=sorted(selected);error=abs(np.interp(Z,Z[ix],T[ix])-T);i=int(error.argmax())
        if error[i]<1e-4:break
        selected.add(i)
    ix=sorted(selected)
    return expr(Z[ix],T[ix],variable='zeta')

def field_expressions(meta,profile):
    name=meta['formulation'];P=meta['pressure_Pa'];scale=meta['initial_flame']['length_scale']
    z=np.array(profile['coordinate_m']);y=z*scale;v=np.array(profile['values'])
    progress=remap_progress(profile);ts0=profile['surface_temperature_K'];tf0=profile['product_temperature_K']
    TsM=meta['initial_flame']['surface_temperature_seed_K'];w=meta['material_properties']['matrix_AP_mass_fraction']
    # Gross §4.3 gives 830 K for AP at 60 atm. Use it only as an initial
    # guess, scaled with the calibrated AP E/R; never impose a surface BC.
    TsAP=10739.15322108806/(10739.15322108806/830-.825*np.log(meta['pressure_atm']/60))
    prop=meta['material_properties'];jAP=1950*1179.5847619198041*np.exp(-10739.15322108806/TsAP)
    jM=prop['matrix_density_kg_m3']*prop['matrix_A_m_s']*np.exp(-prop['matrix_activation_temperature_K']/TsM)
    modes,width,area=surface_modes(meta);mean_j=area*jAP+(1-area)*jM
    U=mean_j*R*1800/P
    conductivity=(2.13037e-7*1800+5.32743e-6)*418.4
    D=conductivity/((P/(R*1800))*(.3*4184))
    sigma=2*meta['dx_m']
    alpha=f'min(1,max(0,{fourier(modes,width,sigma,lambda k:k*k*D/U)}))'
    surface=f'min(1,max(0,{fourier(modes,width,sigma,lambda k:0)}))'
    # A divergence-free mass-flux extension of the heterogeneous surface
    # source; nonzero Fourier modes decay above the surface.
    mass=modes*(jAP-jM);mass[0]+=jM
    jy=fourier(mass,width,sigma,lambda k:abs(k))
    transverse=-1j*mass;transverse[0]=0
    jx=fourier(transverse,width,sigma,lambda k:abs(k))
    s=expr(y,np.clip((v[:,0]-300)/(ts0-300),0,1))
    g=expr(y,np.clip((v[:,0]-ts0)/(tf0-ts0),0,1))
    # Final reaction extents determine an adiabatic enthalpy estimate for
    # each local mixture, including the unchanged mass-specific condensed Q.
    end=progress[-1]
    prefix=f'a=({alpha}); zeta=(a*{jAP:.17g}+(1-a)*{jM*w:.17g})/(a*{jAP:.17g}+(1-a)*{jM:.17g}); '
    prefix+=f'tf=({adiabatic_expr(end)}); '
    prefix+=f'ts={TsM:.17g}+({TsAP-TsM:.17g})*({surface}); '
    prefix+=f'temp=if(y<0,300+(ts-300)*({s}),ts+(max(tf,ts)-ts)*({g})); '
    occupancy=f'(0.5+0.5*tanh(2*y/{meta["interface_width_m"]:.17g}))'
    temperature=prefix+'temp'
    velocities=[prefix+f'{occupancy}*({j})*{R:.17g}*temp/{P:.17g}' for j in (jx,jy)]
    prefix+=f'g3=({expr(y,progress[:,3])}); '
    prefix+=f'e2=min(zeta/{BETA:.17g},1-zeta)*({expr(y,progress[:,2])}); '
    prefix+=f'e0=(zeta-{BETA:.17g}*e2)*({expr(y,progress[:,0])}); '
    prefix+=f'e1=(1-zeta-e2)*({expr(y,progress[:,1])}); e3=min(e0/{GAMMA:.17g},e1)*g3; '
    fractions=(f'zeta-e0-{BETA:.17g}*e2','1-zeta-e1-e2',f'e0-{GAMMA:.17g}*e3','e1-e3',f'{BETA+1:.17g}*e2',f'{GAMMA+1:.17g}*e3')
    densities=[prefix+f'{occupancy}*{P:.17g}/({R:.17g}*temp)*max(0,{f})' for f in fractions]
    info=dict(method='Conserved-mixture diffusion approximation; mapped four-reaction progress; EOS-consistent density and potential mass flux. Not a solved steady flame.',
        source_AP_surface_area_fraction=area,AP_surface_temperature_seed_K=float(TsAP),matrix_surface_temperature_seed_K=float(TsM),
        AP_surface_mass_flux_seed_kg_m2_s=float(jAP),matrix_surface_mass_flux_seed_kg_m2_s=float(jM),
        mean_mass_flux_seed_kg_m2_s=float(mean_j),mixing_diffusivity_m2_s=float(D),mixing_speed_m_s=float(U),
        surface_AP_temperature_source='Gross 2013 section 4.3: 830 K at 60 atm; initial guess only. All solid kinetic parameters remain calibrated values.')
    return temperature,velocities,densities,info

def main():
    profile=json.loads((HERE/'initial_conditions/flame_seed_source.json').read_text())
    archive=HERE/'production_cases_12_planar_archived.txt'
    if not archive.exists():shutil.copy2(HERE/'production_cases.txt',archive)
    preview_archive=HERE/'preview_cases_planar_archived.txt'
    if not preview_archive.exists():shutil.copy2(HERE/'preview_cases.txt',preview_archive)
    originals=(HERE/'initial_conditions/planar_full_height');originals.mkdir(exist_ok=True)
    for p in (HERE/'initial_conditions').glob('*.png'):
        if not (originals/p.name).exists():shutil.copy2(p,originals/p.name)
    production=[];previews=[];records=[]
    for path in archive.read_text().splitlines():
        source=ROOT/path;m=json.loads((source/'case.json').read_text());d=deck_read(source/'input')
        case=source.with_name(source.name.replace('_warm','_compact_diffusion'))
        nx=int(d['amr.n_cell'].split()[0]);coarse=m['packing']['width_mm']*.001/nx;block=4*coarse
        bed={'M03':150e-6,'M17':350e-6,'M21':1200e-6,'M24':650e-6}[m['formulation']]
        gas=max(100e-6,150e-6*(m['pressure_atm']/34)**(-.5))
        bed=math.ceil(bed/block)*block;gas=math.ceil(gas/block)*block
        d['amr.n_cell']=f'{nx} {round((bed+gas)/coarse)}';d['geometry.prob_lo']=f'0 {-bed:.17g} 0'
        d['geometry.prob_hi']=f'{nx*coarse:.17g} {gas:.17g} 0'
        d['plot_file']=str(case.relative_to(ROOT)/'output');d['run_control.stop_file']=str(case.relative_to(ROOT)/'STOP')
        d['diagnostics.interval']='100'
        temp,vel,density,info=field_expressions(m,profile)
        d['temperature.ic.expression.region0']=f'"{temp}"'
        for i,f in enumerate(vel):d[f'velocity.ic.expression.region{i}']=f'"{f}"'
        for species,f in zip(SPECIES,density):d[f'{species}.density.ic.expression.region0']=f'"{f}"'
        m.update(initial_condition='approximate_local_diffusion_flames',initialization_model=info,
            shortened_domain=dict(initial_solid_depth_m=bed,initial_gas_height_m=gas,
                preserved_horizontal_width_m=nx*coarse,preserved_finest_dx_m=m['dx_m'],
                geometry='Crop the same seed-101 disk field; no radius, horizontal position, or mixture-property change.',
                boundary_validation='Pending: require cold solid buffer below deepest interface and negligible reacting outflow; compare enlarged domain at steady state.'),
            predecessor_case=source.name,rate_status='not_run')
        m=save_case(case,d,m);production.append(str(case.relative_to(ROOT)));records.append(m)
        if 20<m['pressure_atm']<50:
            preview=case.with_name(case.name+'_initial_v2');pd=dict(d)
            pd['max_step']='0';pd['plot_file']=str(preview.relative_to(ROOT)/'output');pd['run_control.stop_file']=str(preview.relative_to(ROOT)/'STOP')
            save_case(preview,pd,dict(m,pilot=True,initialization_only=True,production_case=str(case.relative_to(ROOT))))
            previews.append(str(preview.relative_to(ROOT)))
    (HERE/'production_cases.txt').write_text('\n'.join(production)+'\n')
    (HERE/'preview_cases.txt').write_text('\n'.join(previews)+'\n')
    (HERE/'initial_conditions/compact_cases.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps([dict(case=m['case'],**m['shortened_domain']) for m in records],indent=2))

if __name__=='__main__':main()
