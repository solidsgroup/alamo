#!/usr/bin/env python3
"""Use an entire smaller M21 tile to preserve its three AP populations."""
import json
from pathlib import Path
import numpy as np
from generate_packs import generate,properties,PARAMETERS
from prepare_compact import HERE,ROOT,deck_read,save_case,field_expressions,SPECIES

def main():
    counts={400:5,200:20,50:107};nx=128;nsolid=72
    area=sum(n*np.pi*(d*.0005)**2 for d,n in counts.items())
    width_mm=np.sqrt(area*nx/(nsolid*properties('M21')['resolved_AP_area_fraction']))
    directory=HERE/'packs_compact';saved=directory/'M21_seed101.json'
    if saved.exists():pack=json.loads(saved.read_text())
    else:pack,_,_=generate('M21',101,output_dir=directory,counts=counts,width_mm=width_mm)
    assert abs(pack['bed_height_mm']/pack['width_mm']-nsolid/nx)<1e-12
    profile=json.loads((HERE/'initial_conditions/flame_seed_source.json').read_text())
    paths=(HERE/'production_cases.txt').read_text().splitlines();new_previews=[]
    for path in paths:
        case=ROOT/path;m=json.loads((case/'case.json').read_text())
        if m['formulation']!='M21':continue
        assert not (case/'run.json').exists()
        d=deck_read(case/'input');coarse=pack['width_mm']*.001/nx;dx=coarse/16;width=4*dx
        m.update(packing=pack,packing_file=str((directory/'M21_seed101.xyzr').relative_to(ROOT)),dx_m=dx,interface_width_m=width)
        bed=pack['bed_height_mm']*.001;ngas=int(np.ceil(max(100e-6,150e-6*(m['pressure_atm']/34)**(-.5))/(4*coarse)))*4;gas=ngas*coarse
        d['amr.n_cell']=f'{nx} {nsolid+ngas}';d['geometry.prob_lo']=f'0 {-bed:.17g} 0';d['geometry.prob_hi']=f'{nx*coarse:.17g} {gas:.17g} 0'
        for species in ('AP_solid','HTPB_solid'):
            prefix=species+'.density.ic.psread.'
            d[prefix+'file.name']=m['packing_file'];d[prefix+'eps']=str(2*dx);d[prefix+'plane.eps']=str(width)
        for mech,mat in [('AP_decomposition',PARAMETERS['AP']),('HTPB_pyrolysis',PARAMETERS['binder'])]:
            pre=mech+'.phase_change.';d[pre+'rate_multiplier']=str(mat['A_m_s']/(1.5*width));d[pre+'phase_field.allencahn.kappa']=str(3*width*width)
        d['HTPB_pyrolysis.phase_change.homogeneous.ap_rate_multiplier']=str(PARAMETERS['AP']['A_m_s']/(1.5*width))
        temp,vel,density,info=field_expressions(m,profile)
        d['temperature.ic.expression.region0']=f'"{temp}"'
        for i,f in enumerate(vel):d[f'velocity.ic.expression.region{i}']=f'"{f}"'
        for s,f in zip(SPECIES,density):d[f'{s}.density.ic.expression.region0']=f'"{f}"'
        m.update(initialization_model=info,shortened_domain=dict(initial_solid_depth_m=bed,initial_gas_height_m=gas,
            horizontal_width_m=nx*coarse,finest_dx_m=dx,
            geometry='Entire smaller periodic seed-101 tile: 5x400um, 20x200um, 107x50um disks. Width adjusted by <1% to fit the tile exactly on a square-cell grid.',
            boundary_validation='Pending: check cold solid buffer and reacting outflow; enlarged-domain steady comparison required.'))
        save_case(case,d,m)
        if 20<m['pressure_atm']<50:
            preview=case.with_name(case.name+'_initial_v4');pd=dict(d,max_step='0')
            pd['plot_file']=str(preview.relative_to(ROOT)/'output');pd['run_control.stop_file']=str(preview.relative_to(ROOT)/'STOP')
            save_case(preview,pd,dict(m,pilot=True,initialization_only=True,production_case=path))
            new_previews.append(str(preview.relative_to(ROOT)))
    manifest=HERE/'preview_cases.txt';rows=manifest.read_text().splitlines()
    rows=[new_previews[0] if Path(p).name.startswith('M21_') else p for p in rows];manifest.write_text('\n'.join(rows)+'\n')
    records=[json.loads((ROOT/p/'case.json').read_text()) for p in paths]
    (HERE/'initial_conditions/compact_cases.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps(dict(width_mm=pack['width_mm'],depth_mm=pack['bed_height_mm'],bins=pack['bins'],preview=new_previews),indent=2))

if __name__=='__main__':main()
