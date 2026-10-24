#!/usr/bin/env python3
"""Reproducible periodic 2D disks with Table 2 mass fractions.

Disk counts approximate each resolved size-bin mass ratio. Box height is set
from the exact disk area and requested resolved area fraction. Disks are
relaxed in a doubly periodic box; vertical images are explicitly included in
the solver file because its vertical simulation boundary is not periodic.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/alamo-gross-matplotlib')
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

HERE=Path(__file__).resolve().parent
PARAMETERS=json.loads((HERE/'parameters_solid_frozen.json').read_text())
TABLE={
 'M03':{20:.5579,.7:.3158},
 'M17':{90:.3158,20:.5579},
 'M21':{400:.3158,200:.3158,50:.1053,20:.1368},
 'M24':{200:.3158,50:.4211,20:.1368},
}
COUNTS={'M03':{},'M17':{90:80},'M21':{400:12,200:48,50:256},'M24':{200:48,50:1024}}
WIDTH_MM={'M03':.032,'M17':.8,'M21':2.,'M24':2.}

def properties(name,resolve_m03=False):
    from scipy.optimize import brentq
    b,ap=PARAMETERS['binder'],PARAMETERS['AP']
    resolved=lambda d: d>20 or (resolve_m03 and name=='M03' and d==20)
    fine=sum(v for d,v in TABLE[name].items() if not resolved(d))
    coarse=sum(v for d,v in TABLE[name].items() if resolved(d))
    binder=1-sum(TABLE[name].values()); matrix=fine+binder;w=fine/matrix
    va=fine/ap['density_kg_m3'];vb=binder/b['density_kg_m3']
    t=va/(va+vb);rho=matrix/(va+vb)
    totalv=va+vb+coarse/ap['density_kg_m3']
    k0,k1=b['conductivity_W_m_K'],ap['conductivity_W_m_K']
    k=brentq(lambda k:k-k1-(1-t)*(k0-k1)*np.sqrt(k/k0),k0,k1)
    return dict(formulation=name,binder_mass_fraction=binder,fine_AP_mass_fraction=fine,
        resolved_AP_mass_fraction=coarse,matrix_AP_mass_fraction=w,matrix_AP_volume_fraction=t,
        resolved_AP_area_fraction=(coarse/ap['density_kg_m3'])/totalv,
        propellant_density_kg_m3=1/totalv,matrix_density_kg_m3=rho,
        matrix_cp_J_kg_K=(1-w)*b['cp_J_kg_K']+w*ap['cp_J_kg_K'],
        matrix_Q_J_kg=(1-w)*b['heat_release_J_kg']+w*ap['heat_release_J_kg'],
        matrix_k_W_m_K=k,matrix_A_m_s=np.exp((1-t)*np.log(b['A_m_s'])+t*np.log(ap['A_m_s'])),
        matrix_activation_temperature_K=(1-t)*b['activation_temperature_K']+t*ap['activation_temperature_K'])

def objective(flat,radii,box):
    x=flat.reshape(-1,2)%box
    pairs=cKDTree(x,boxsize=box).query_pairs(2*radii.max()+.0001,output_type='ndarray')
    grad=np.zeros_like(x)
    if len(pairs)==0:return 0.,grad.ravel()
    a,b=pairs.T;d=x[a]-x[b];d-=box*np.round(d/box)
    dist=np.sqrt(np.sum(d*d,axis=1));over=radii[a]+radii[b]+.0001-dist
    active=over>0;a,b,d,dist,over=a[active],b[active],d[active],dist[active],over[active]
    force=-over[:,None]*d/np.maximum(dist[:,None],1.e-12)
    np.add.at(grad,a,force);np.add.at(grad,b,-force)
    return .5*np.dot(over,over),grad.ravel()

def generate(name,seed,resolve_m03=False,output_dir=None,counts=None,width_mm=None):
    props=properties(name,resolve_m03);rng=np.random.default_rng(seed)
    counts=COUNTS[name] if counts is None else counts
    radii=np.array([d*.0005 for d,n in counts.items() for _ in range(n)]) # mm
    width=WIDTH_MM[name] if width_mm is None else width_mm
    if len(radii):
        height=np.sum(np.pi*radii**2)/(props['resolved_AP_area_fraction']*width)
        box=np.array([width,height]);x=rng.random((len(radii),2))*box
        # Inflate slowly to avoid trapping initially coincident particles.
        for scale in np.linspace(.25,1,16):
            result=minimize(objective,x.ravel(),args=(radii*scale,box),jac=True,method='L-BFGS-B',
                options=dict(maxiter=3000,ftol=1.e-18,gtol=1.e-10,maxcor=20))
            x=result.x.reshape(-1,2)%box
        energy,_=objective(x.ravel(),radii,box)
        a,b=np.triu_indices(len(radii),1);d=x[a]-x[b];d-=box*np.round(d/box)
        min_gap=float(np.min(np.sqrt(np.sum(d*d,axis=1))-radii[a]-radii[b]))
        if min_gap < -1.e-6:raise RuntimeError(f'{name} seed={seed}: residual overlap {min_gap} mm')
        x[:,1]-=height
    else:
        height=.6;x=np.empty((0,2));energy=0.;min_gap=None
    directory=HERE/'packs' if output_dir is None else Path(output_dir)
    directory.mkdir(exist_ok=True)
    stem=directory/f'{name}_seed{seed}'
    np.savetxt(stem.with_suffix('.csv'),np.column_stack((x,radii)),delimiter=',',header='x_mm,y_mm,radius_mm',comments='')
    # Tiling in y allows an unbiased planar cut through disks at y=0 and
    # makes the initialized composition available throughout the solid bed.
    diskrows=[]
    for shift in (-height,0,height):
        diskrows.extend((xx,yy+shift,0.,r) for (xx,yy),r in zip(x,radii))
    if diskrows:np.savetxt(stem.with_suffix('.xyzr'),diskrows,fmt='%.15g')
    bins={}
    for diameter,n in counts.items():
        area=n*np.pi*(diameter*.0005)**2/(width*height)
        mass=area*PARAMETERS['AP']['density_kg_m3']/props['propellant_density_kg_m3']
        bins[str(diameter)]=dict(count=n,diameter_um=diameter,realized_mass_fraction=mass,
                               target_mass_fraction=TABLE[name][diameter],error_mass_fraction=mass-TABLE[name][diameter])
    metadata=dict(**props,seed=seed,width_mm=width,bed_height_mm=height,
        disk_count=len(radii),bins=bins,minimum_gap_mm=min_gap,overlap_energy_mm2=energy,
        geometry='2D disks; periodic x; y-periodic packing tile explicitly repeated',
        independence='Independent random geometry' if len(radii) else 'Identical planar control: no packing variance exists',
        parameter_sha256=hashlib.sha256((HERE/'parameters_solid_frozen.json').read_bytes()).hexdigest())
    stem.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps({k:metadata[k] for k in ('formulation','seed','disk_count','width_mm','bed_height_mm','minimum_gap_mm')}),flush=True)
    return metadata,x,radii

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--formulations',nargs='+',default=list(TABLE))
    parser.add_argument('--seeds',nargs='+',type=int,default=[101,202,303]);args=parser.parse_args()
    (HERE/'packs').mkdir(exist_ok=True)
    rows=[properties(name) for name in TABLE]
    with (HERE/'reference/formulations.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    fig,axes=plt.subplots(len(args.seeds),len(args.formulations),figsize=(12,9),squeeze=False,layout='constrained')
    for j,name in enumerate(args.formulations):
        for i,seed in enumerate(args.seeds):
            meta,x,r=generate(name,seed);ax=axes[i,j]
            ax.set_facecolor('#efd2a2')
            for offset in (-meta['width_mm'],0,meta['width_mm']):
                for (xx,yy),rr in zip(x,r):ax.add_patch(Circle((xx+offset,yy),rr,facecolor='#406ca0',lw=.2,edgecolor='white'))
            ax.set(xlim=(0,meta['width_mm']),ylim=(-meta['bed_height_mm'],0),title=f'{name}, seed {seed}',xlabel='x (mm)',ylabel='y (mm)')
            ax.set_aspect('equal')
            if not len(r):ax.text(.5,.5,'Fully homogenized',ha='center',rotation=90,transform=ax.transAxes)
    fig.savefig(HERE/'packs/pack_overview.png',dpi=180);plt.close(fig)

if __name__=='__main__':main()
