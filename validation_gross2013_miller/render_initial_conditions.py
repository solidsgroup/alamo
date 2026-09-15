#!/usr/bin/env python3
"""Render actual zero-step AMReX flame/packing fields for user review."""
import gc
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/alamo-gross-matplotlib')
import numpy as np
import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch,Circle,Rectangle
from PIL import Image,ImageOps,ImageDraw
yt.set_log_level(50)
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SPECIES=('AP_gas','HTPB_gas','Mono','Premixed','Primary','Final')
COLORS=('#0072b2','#d55e00','#009e73','#cc79a7','#8c564b','#555555')
AP=np.array([.24,.43,.65]);MATRIX=np.array([.92,.80,.60])
norm=Normalize(300,3400);cmap=plt.get_cmap('inferno')

def covering(ds,level,ybounds=None):
    dims=ds.domain_dimensions*ds.refine_by**level;dims[2:]=1
    left=ds.domain_left_edge.copy();dy=float(ds.domain_width[1])/dims[1]
    if ybounds is not None:
        lo=max(0,int(np.floor((ybounds[0]-float(left[1]))/dy)))
        hi=min(int(dims[1]),int(np.ceil((ybounds[1]-float(left[1]))/dy)))
        left[1]=float(left[1])+lo*dy;dims[1]=hi-lo
    cg=ds.covering_grid(level,left,dims)
    dx=float(ds.domain_width[0])/dims[0]
    extent=(float(left[0])*1e6,(float(left[0])+dims[0]*dx)*1e6,float(left[1])*1e6,(float(left[1])+dims[1]*dy)*1e6)
    return cg,extent,dy

def render(case):
    meta=json.loads((case/'case.json').read_text());receipt=json.loads((case/'run.json').read_text())
    assert receipt['returncode']==0
    source=case/'output/00000cell';ds=yt.load(str(source))
    # Read only exact domain cells; bypass a several-ulp yt right-edge guard.
    ds.force_periodicity()
    ad=ds.all_data();area=ad['index','dx'].d*ad['index','dy'].d
    gas=sum(ad['boxlib','component_density_'+s].d for s in SPECIES)
    T=ad['boxlib','temperature'].d;eta=ad['boxlib','rigid_eta'].d
    P=meta['pressure_Pa'];R=8.31446261815324/.026
    residual=eta+gas*R*T/P-1
    assert np.all(np.isfinite(T)) and T.min()>=299.9
    assert np.all(np.isfinite(gas)) and gas.min()>=-1e-12
    assert float(ds.current_time)==0
    width=float(ds.domain_width[0]);rho_ap=ad['boxlib','component_density_AP_solid'].d
    audit=dict(formulation=meta['formulation'],case=case.name,pressure_atm=meta['pressure_atm'],
        time_s=0.,max_AMR_level=ds.index.max_level,valid_cells=len(area),
        temperature_min_K=float(T.min()),temperature_max_K=float(T.max()),
        max_abs_volume_constraint_residual=float(np.max(abs(residual))),
        resolved_AP_solid_area_fraction=float(np.sum(rho_ap/1950*area)/np.sum(eta*area)),
        tile_target_resolved_AP_area_fraction=meta['packing']['resolved_AP_area_fraction'],
        gas_heat_release_W_m2=float(np.sum(ad['boxlib','qdot'].d*area)/width),
        external_heat_source_present='heat_source.ic' in (case/'input').read_text(),
        reference_note='Field values from actual zero-step solver output. Full-bed area fraction includes the rounded simulation bed and coarse-cell disk sampling.')
    assert not audit['external_heat_source_present']
    cg,extent,_=covering(ds,min(2,ds.index.max_level))
    e=cg['boxlib','rigid_eta'].d[:,:,0]
    a=cg['boxlib','component_density_AP_solid'].d[:,:,0]/1950
    temp=cg['boxlib','temperature'].d[:,:,0]
    fraction=np.clip(a/np.maximum(e,1e-15),0,1)
    solid_rgb=fraction[:,:,None]*AP+(1-fraction[:,:,None])*MATRIX
    rgb=np.clip(e,0,1)[:,:,None]*solid_rgb+(1-np.clip(e,0,1))[:,:,None]*cmap(norm(temp))[:,:,:3]
    zoom,zoom_extent,dy=covering(ds,ds.index.max_level,(-100e-6,160e-6))
    ze=zoom['boxlib','rigid_eta'].d[:,:,0];zt=zoom['boxlib','temperature'].d[:,:,0]
    za=zoom['boxlib','component_density_AP_solid'].d[:,:,0]/1950
    zq=zoom['boxlib','qdot'].d[:,:,0]
    y=np.linspace(zoom_extent[2]+dy*5e5,zoom_extent[3]-dy*5e5,zt.shape[1])
    x=np.linspace(zoom_extent[0],zoom_extent[1],zt.shape[0])
    fig=plt.figure(figsize=(12,7.3),layout='constrained');grid=fig.add_gridspec(2,3,width_ratios=(1.12,1,1))
    full=fig.add_subplot(grid[:,0]);zoom_ax=fig.add_subplot(grid[0,1:]);thermal=fig.add_subplot(grid[1,1]);species=fig.add_subplot(grid[1,2])
    # Show the exact input disks, since cold solids intentionally remain coarse
    # in AMR and are reconstructed from these disks as the front approaches.
    geometry_rgb=cmap(norm(temp))[:,:,:3]
    yc=np.linspace(extent[2],extent[3],temp.shape[1])
    geometry_rgb[:,yc<0]=MATRIX
    full.imshow(np.transpose(geometry_rgb,(1,0,2)),origin='lower',extent=extent,aspect='equal',interpolation='nearest')
    disk_path=ROOT/meta['packing_file'] if meta.get('packing_file') else HERE/'packs_reduced'/f"{meta['formulation']}_seed101.xyzr"
    disks=np.loadtxt(disk_path)*1000
    def draw_disks(ax,filled):
        bounds=extent if filled else zoom_extent
        bottom,top=bounds[2],min(0,bounds[3])
        clip=Rectangle((bounds[0],bottom),bounds[1]-bounds[0],top-bottom,transform=ax.transData)
        ax.autoscale(False)
        for shift in (-(extent[1]-extent[0]),0,extent[1]-extent[0]):
            for xx,yy,_,radius in disks:
                if yy-radius>top or yy+radius<bottom:continue
                if xx+shift+radius<bounds[0] or xx+shift-radius>bounds[1]:continue
                c=Circle((xx+shift,yy),radius,facecolor=AP if filled else 'none',edgecolor='white' if filled else '#a8d7ff',lw=.25 if filled else .5)
                ax.add_patch(c);c.set_clip_path(clip)
    draw_disks(full,True)
    full.axhline(0,color='white',lw=.7)
    full.set(xlabel='x (µm)',ylabel='y (µm)',title='Input disk geometry + gas temperature')
    full.legend(handles=[Patch(facecolor=AP,label='Resolved AP'),Patch(facecolor=MATRIX,label='Binder + fine AP')],loc='lower left',fontsize=8,framealpha=.95)
    im=zoom_ax.imshow(zt.T,origin='lower',extent=zoom_extent,aspect='auto',cmap=cmap,norm=norm,interpolation='nearest')
    zoom_ax.contour(x,y,ze.T,levels=[.5],colors='white',linewidths=.8)
    draw_disks(zoom_ax,False)
    zoom_ax.set_xlim(zoom_extent[:2]);zoom_ax.set_ylim(zoom_extent[2:])
    if zq.max()>0:zoom_ax.contour(x,y,zq.T,levels=[.1*zq.max(),.5*zq.max()],colors='#54ffbc',linewidths=.65)
    zoom_ax.set(xlabel='x (µm)',ylabel='y (µm)',title='Surface detail: temperature; green contours mark gas reaction zone')
    fig.colorbar(im,ax=zoom_ax,label='Temperature (K)',fraction=.03,pad=.02)
    thermal.fill_between(y,zt.min(0),zt.max(0),color='#b34012',alpha=.15)
    thermal.plot(y,zt.mean(0),color='#b34012',label='Mean temperature; band: range across x')
    thermal.axvline(0,color='.4',lw=.8,ls='--');thermal.set(xlabel='Height above initial surface (µm)',ylabel='Temperature (K)',xlim=(-70,min(140,zoom_extent[3])));thermal.grid(alpha=.2)
    thermal.legend(fontsize=6,loc='lower right')
    velocity=zoom['boxlib','velocityy'].d[:,:,0].mean(0)
    other=thermal.twinx();other.plot(y,velocity,color='#176c98',ls='--');other.set(ylabel='Gas velocity y (m/s)')
    yy=np.stack([zoom['boxlib','component_density_'+s].d[:,:,0].mean(0) for s in SPECIES]);yy/=np.maximum(yy.sum(0),1e-100)
    for s,color,f in zip(SPECIES,COLORS,yy):species.plot(y,np.where(y>=0,f,np.nan),label=s,color=color,lw=1.3)
    species.set(xlabel='Height above initial surface (µm)',ylabel='Gas mass fraction',xlim=(0,100),ylim=(0,1));species.grid(alpha=.2);species.legend(fontsize=7,ncol=2)
    sizes=', '.join(f'{int(float(d))}' for d in sorted(meta['packing']['bins'],key=float))
    fine='0.7 µm' if meta['formulation']=='M03' else '≤20 µm'
    seed='Approximate local diffusion-flame seed' if meta['initial_condition']=='approximate_local_diffusion_flames' else 'Initialized planar flame'
    fig.suptitle(f"{meta['formulation']} · {meta['pressure_atm']:.1f} atm · seed 101 · t = 0\nResolved AP: {sizes} µm; homogenized AP: {fine}. {seed}; no external heating.",fontsize=11)
    out=HERE/'initial_conditions';out.mkdir(exist_ok=True)
    fig.savefig(out/f"{meta['formulation']}_initial_condition.png",dpi=190);plt.close(fig)
    (out/f"{meta['formulation']}_initial_audit.json").write_text(json.dumps(audit,indent=2)+'\n')
    return audit

def main():
    audits=[]
    for path in (HERE/'preview_cases.txt').read_text().splitlines():
        audits.append(render(ROOT/path));gc.collect()
    out=HERE/'initial_conditions';tile=(1140,695);sheet=Image.new('RGB',(2*tile[0],2*tile[1]),'white')
    for i,name in enumerate(('M03','M17','M21','M24')):
        with Image.open(out/f'{name}_initial_condition.png') as im:
            thumb=ImageOps.contain(im.convert('RGB'),tile,Image.Resampling.LANCZOS)
            sheet.paste(thumb,((i%2)*tile[0],(i//2)*tile[1]))
    sheet.save(out/'initial_conditions_overview.png')
    (out/'initialization_audit.json').write_text(json.dumps(audits,indent=2)+'\n')
    print(json.dumps(audits,indent=2))

if __name__=='__main__':main()
