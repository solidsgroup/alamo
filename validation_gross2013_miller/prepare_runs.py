#!/usr/bin/env python3
"""Generate reactive input decks; never launch simulations or fit their rates."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from generate_packs import HERE, PARAMETERS, properties

ROOT=HERE.parent
GAS_MODEL=dict(
    source='Existing Alamo Rocfire implementation / Gross (2013) Table 1',
    A_g_cm3_s=[27000,25300,5000,400],pressure_exponent=[1.6,1.7,1.7,1.7],
    activation_energy_kcal_mol=[15.9,14.5,11.3,14.5],
    qgas_cal_g=[430,447,8365,2127],qsolid_cal_g=[0,0],beta=7.974,gamma=3.19,
    cp_cal_g_K=.3,molecular_weight_g_mol=26,pressure_reference_Pa=100000,
    approximation='Fixed gas coefficients for every formulation/pressure. Gross composition-dependent adiabatic temperature inputs and fine-particle correction table are unavailable.',
    homogeneous_gas_products='Mass split into AP_gas and HTPB_gas, as implemented by the existing homogeneous phase-change model. This differs from Gross composition-specific pseudo-binder gas.',
    heat_accounting='Calibrated condensed Q applied once on mass transfer; legacy gas qsolid disabled.')

def generate(name,seed,p_atm,reference_rate,pilot=False,dx_um=1.,width_cells=4,tag='',ignition_off=.0003,pack_dir=None,resolve_m03=False):
    if pilot and name != 'M03':
        raise ValueError('Narrow planar pilots require M03; clone a full pack deck for a packed pilot')
    pack_dir=HERE/'packs' if pack_dir is None else Path(pack_dir)
    prop=properties(name,resolve_m03);pack=json.loads((pack_dir/f'{name}_seed{seed}.json').read_text())
    resolved_pack=pack['disk_count']>0
    b,ap=PARAMETERS['binder'],PARAMETERS['AP'];P=p_atm*101325
    if pilot:
        nx=round(4/dx_um);dx=dx_um*1.e-6;ny=round(512/dx_um)
        lo_y=-128.e-6;hi_y=384.e-6;levels=0
        width=width_cells*dx;stop=.003;plot_dt=stop/60
    else:
        nx=128 if name in ('M21','M24') else (64 if resolved_pack else 4)
        coarse_dx=pack['width_mm']*.001/nx
        levels=4 if resolved_pack else 0
        dx=coarse_dx/2**levels;width=4*dx
        lo_y=-math.ceil(pack['bed_height_mm']*.001/coarse_dx/8)*8*coarse_dx
        hi_y=math.ceil(.0005/coarse_dx/8)*8*coarse_dx
        ny=round((hi_y-lo_y)/coarse_dx)
        target={'M03':.0002,'M17':.0008,'M21':.0016,'M24':.0012}[name]
        stop=.0015+target/(reference_rate*.01);plot_dt=stop/100
        if not resolved_pack:
            dx=1.e-6;width=4*dx;lo_y=-.0004;hi_y=.0004;ny=800
    v0=1.5*width
    case_name=f'{name}_p{p_atm:.4f}_seed{seed}'+(f'_pilot_dx{dx_um:g}_w{width_cells:g}' if pilot else '')+tag
    case=HERE/'runs'/case_name;case.mkdir(parents=True,exist_ok=True)
    if (case/'run.json').exists() or (case/'stdout.log').exists():
        raise FileExistsError(f'Preserving launched case {case}; use a new tag')
    d={
      'plot_file':str(case.relative_to(ROOT)/'output'),'amr.plot_int':2147483647,'amr.plot_dt':plot_dt,
      'amr.max_level':levels,'amr.n_cell':f'{nx} {ny}','amr.max_grid_size':64,'amr.blocking_factor':4,
      'amr.regrid_int':5 if levels else -1,'amr.grid_eff':.7,'amr.nsubsteps':1,
      'amr.reinitialize_condensed_composition':1,'amr.reinitialize_condensed_composition_eta_min':.99,
      'eta_refinement_criterion':.1,'reaction_refinement_criterion':'1.e4_1/s','temperature_refinement_criterion':100,
      'max_step':2147483647,'stop_time':stop,'timestep':'1.e-8_s',
      'dynamictimestep.on':1,'dynamictimestep.min':1.e-14,'dynamictimestep.max':5.e-7,
      'integration.type':'RungeKutta','integration.rk.type':3,
      'geometry.prob_lo':f'0 {lo_y:.17g} 0','geometry.prob_hi':f'{nx*dx*2**levels:.17g} {hi_y:.17g} 0',
      'geometry.is_periodic':'1 0 0','system.amount':'kmol',
      'gas.mw':' '.join(['26_g/mol']*6),'gas.thermo.type':'rocfire','gas.thermo.rocfire.mw':'26_g/mol',
      'gas.thermo.rocfire.cp_mass':'.3_cal/g/K','gas.thermo.rocfire.Tref':'300_K',
      'gas.transport.type':'rocfire','gas.transport.rocfire.mw':'26_g/mol',
      'gas.transport.rocfire.prandtl':.69,'gas.transport.rocfire.lewis':1,
      'gas.transport.rocfire.lambda_a':'2.13037e-7_cal/cm/s/K^2',
      'gas.transport.rocfire.lambda_b':'5.32743e-6_cal/cm/s/K','gas.eos.type':'tpg',
      'species.names':'AP_gas HTPB_gas Mono Premixed Primary Final AP_solid HTPB_solid',
      'chemistry.model.type':'rocfire','chemistry.model.rocfire.A':' '.join(f'{a}_g/cm^3/s' for a in GAS_MODEL['A_g_cm3_s']),
      'chemistry.model.rocfire.pressure_exponent':'1.6 1.7 1.7 1.7',
      'chemistry.model.rocfire.activation_energy':'15.9_kcal/mol 14.5_kcal/mol 11.3_kcal/mol 14.5_kcal/mol',
      'chemistry.model.rocfire.qgas':'430_cal/g 447_cal/g 8365_cal/g 2127_cal/g',
      'chemistry.model.rocfire.qsolid':'0_cal/g 0_cal/g','chemistry.model.rocfire.beta':7.974,
      'chemistry.model.rocfire.gamma':3.19,'chemistry.model.rocfire.pressure_reference':'1_bar',
      'chemistry.solver.type':'backward_euler','chemistry.solver.backward_euler.max_iter':100,
      'chemistry.solver.backward_euler.reltol':1.e-8,'chemistry.solver.backward_euler.abstol':1.e-12,
      'mechanisms.names':'AP_decomposition HTPB_pyrolysis','rigid.relaxation_time':'1.e-8_s','rigid.velocity':'0 0',
      'advection.type':'muscl','advection.muscl.limiter.type':'minmod','cfl':.4,'phase_field.cfl':.5,
      'include_viscosity':1,'implicit_viscosity':1,'include_conduction':1,'advect_temperature':1,
      'projection.enabled':1,'projection.tol_rel':1.e-9,'projection.tol_abs':1.e-12,
      'projection.verbose':0,'projection.update_pressure':1,
      'projection.bottom_solver':'bicgstab','projection.bottom_max_iter':1000,
      'diffusion.tol_rel':1.e-10,'diffusion.tol_abs':1.e-12,'diffusion.verbose':0,
      'diagnostics.interval':1000,'diagnostics.extended_fields':1,
      'velocity.ic.type':'constant','velocity.ic.constant.value':'0 0',
      'temperature.ic.type':'constant','temperature.ic.constant.value':300,
      'pressure.ic.type':'constant','pressure.ic.constant.value':P,
      'Final.density.ic.type':'expression',
      'Final.density.ic.expression.region0':f'"(0.5+0.5*tanh(2*y/{width:.17g}))*{P:.17g}/(319.787*300)"',
      'heat_source.ic.type':'expression','heat_source.ic.expression.unit':'W/m^3',
      'heat_source.ic.expression.region0':f'"3.e7*exp(-(y/{max(width,6.25e-6):.17g})^2)/(1.772453850905516*{max(width,6.25e-6):.17g})*(t<{ignition_off:.17g})"',
    }
    for s in ('AP_gas','HTPB_gas','Mono','Premixed','Primary','Final'):d[s+'.mechanics']='fluid'
    for s,mat in [('AP_solid',ap),('HTPB_solid',b)]:
        d[s+'.mechanics']='rigid_solid';d[s+'.reference_density']=mat['density_kg_m3']
        d[s+'.specific_heat']=mat['cp_J_kg_K'];d[s+'.thermal_conductivity']=mat['conductivity_W_m_K']
    for mech,s,g,mat in [('AP_decomposition','AP_solid','AP_gas',ap),('HTPB_pyrolysis','HTPB_solid','HTPB_gas',b)]:
        d[mech+'.type']='phase_change';prefix=mech+'.phase_change.'
        d[prefix+'in']=s;d[prefix+'out']=g;d[prefix+'rate_multiplier']=mat['A_m_s']/v0
        d[prefix+'activation_temperature']=mat['activation_temperature_K'];d[prefix+'temperature_cutoff']=0
        d[prefix+'phase_field.type']='allencahn'
        for key,value in dict(mobility=1,**{'lambda':1},kappa=3*width**2,w0=0,w12=2,w1=1).items():
            d[prefix+'phase_field.allencahn.'+key]=value
    d['AP_decomposition.phase_change.heat_release']=ap['heat_release_J_kg']
    prefix='HTPB_pyrolysis.phase_change.'
    for key,value in {'homogeneous':'true','homogeneous.ap_solid':'AP_solid','homogeneous.ap_gas':'AP_gas',
        'homogeneous.mass_fraction':prop['matrix_AP_mass_fraction'],
        'homogeneous.ap_rate_multiplier':ap['A_m_s']/v0,
        'homogeneous.ap_activation_temperature':ap['activation_temperature_K'],
        'homogeneous.binder_heat_release':b['heat_release_J_kg'],
        'homogeneous.ap_heat_release':ap['heat_release_J_kg']}.items():d[prefix+key]=value
    for s,rho in [('AP_solid',ap['density_kg_m3']),('HTPB_solid',prop['matrix_density_kg_m3'])]:
        if not resolved_pack:
            d[s+'.density.ic.type']='expression'
            d[s+'.density.ic.expression.region0']=f'"{rho if s=="HTPB_solid" else 0:.17g}*(0.5-0.5*tanh(2*y/{width:.17g}))"'
        else:
            prefix=s+'.density.ic.psread.';d[s+'.density.ic.type']='psread'
            d[prefix+'file.name']=str((pack_dir/f'{name}_seed{seed}.xyzr').relative_to(ROOT))
            d[prefix+'file.unit']='mm';d[prefix+'eps']=2*dx;d[prefix+'value']=rho
            d[prefix+'invert']=int(s=='HTPB_solid');d[prefix+'plane.point']='0 0 0'
            d[prefix+'plane.normal']='0 -1 0';d[prefix+'plane.eps']=width
    for field in ('component_density','velocity','temperature','pressure'):
        d[field+'.bc.type']='constant';prefix=field+'.bc.constant.'
        d[prefix+'type.xlo']='periodic';d[prefix+'type.xhi']='periodic'
        d[prefix+'type.ylo']='neumann';d[prefix+'type.yhi']='neumann'
    d['velocity.bc.constant.type.ylo']='dirichlet';d['velocity.bc.constant.val.ylo']='0 0'
    d['temperature.bc.constant.type.ylo']='dirichlet';d['temperature.bc.constant.val.ylo']=300
    d['pressure.bc.constant.type.yhi']='dirichlet';d['pressure.bc.constant.val.yhi']=P
    text=f'# Generated by validation_gross2013_miller/prepare_runs.py.\n# Reactive fixed-coefficient four-flame approximation; heating is exactly off at {ignition_off*1.e3:g} ms.\n'
    text+='\n'.join(f'{k} = {v:.17g}' if isinstance(v,float) else f'{k} = {v}' for k,v in d.items())+'\n'
    (case/'input').write_text(text)
    meta=dict(case=case_name,formulation=name,seed=seed,pressure_atm=p_atm,pressure_Pa=P,
        reference_rate_cm_s=reference_rate,pilot=pilot,dx_m=dx,interface_width_m=width,
        ignition_off_s=ignition_off,requested_stop_s=stop,plot_dt_s=plot_dt,
        rate_status='not_run',material_properties=prop,packing=pack,
        gas_model=GAS_MODEL,input_sha256=hashlib.sha256(text.encode()).hexdigest(),
        solid_parameters_sha256=hashlib.sha256((HERE/'parameters_solid_frozen.json').read_bytes()).hexdigest())
    (case/'case.json').write_text(json.dumps(meta,indent=2)+'\n')
    return str(case.relative_to(ROOT))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--pilot',action='store_true')
    parser.add_argument('--dx-um',type=float,default=1.);parser.add_argument('--width-cells',type=int,default=4)
    parser.add_argument('--formulations',nargs='+',default=['M03']);parser.add_argument('--tag',default='')
    parser.add_argument('--ignition-off-ms',type=float,default=.3)
    args=parser.parse_args();cases=[]
    if args.pilot:
        for name in args.formulations:cases.append(generate(name,101,68.,1.,True,args.dx_um,args.width_cells,args.tag,args.ignition_off_ms*.001))
    else:
        rows=list(csv.DictReader((HERE/'reference/figure10_digitized.csv').open()))
        for row in rows:
            if row['kind']!='miller_diamond' or row['formulation'] not in args.formulations:continue
            for seed in (101,202,303):cases.append(generate(row['formulation'],seed,float(row['pressure_atm']),float(row['rate_cm_s']),ignition_off=args.ignition_off_ms*.001))
    (HERE/'gas_parameters.json').write_text(json.dumps(GAS_MODEL,indent=2)+'\n')
    print('\n'.join(cases))

if __name__=='__main__':main()
