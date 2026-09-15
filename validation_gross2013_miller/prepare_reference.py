#!/usr/bin/env python3
"""Digitize the original embedded raster of Gross (2013), Fig. 10.

Coordinates are manually audited marker centers in the 1153 x 1068 image,
not unpublished experimental measurements. No burning-rate model is fitted.
"""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

HERE = Path(__file__).resolve().parent
PDF = Path.home() / 'Downloads/gross_fourflame_model.pdf'
IMAGE = HERE / 'reference/figure10_source-000.jpg'
# x=1,10,100 atm; y=1,10 cm/s (M21: y=.1,1 cm/s).
AXES = {
    'M03': dict(x1=80., x100=464., y1=274., y10=31.),
    'M17': dict(x1=675., x100=1059., y1=273., y10=31.),
    'M21': dict(x1=79., x100=464., y1=737., y10=520.),
    'M24': dict(x1=676., x100=1061., y1=827., y10=585.),
}
CENTERS = {
    'M03': [(240,333),(330,245),(375,190),(404,174),(432,138),(490,90),(524,67)],
    'M17': [(835,353),(925,261),(970,230),(1000,211),(1027,195),(1085,160),(1119,135)],
    'M21': [(239,855),(330,803),(374,784),(404,769),(432,754),(489,727),(523,705)],
    'M24': [(836,926),(926,858),(970,828),(1000,818),(1028,796),(1086,763),(1120,747)],
}
# Vertices on the solid "New" line, avoiding diamonds where possible.
SOLID = {
    'M03': [(240,329),(330,231),(375,190),(404,166),(432,141),(490,88),(524,60)],
    'M17': [(835,348),(925,266),(970,234),(1000,214),(1027,195),(1085,160),(1119,136)],
    'M21': [(239,859),(330,806),(374,778),(404,765),(432,753),(489,728)],
    'M24': [(836,928),(926,865),(970,831),(1000,808),(1028,789),(1086,752),(1120,731)],
}

def convert(name, points):
    a = AXES[name]
    p = np.array(points, dtype=float)
    return np.column_stack((10**(2*(p[:,0]-a['x1'])/(a['x100']-a['x1'])),
                            10**((a['y1']-p[:,1])/(a['y1']-a['y10']))))

def main():
    rows=[]
    for kind, groups in [('miller_diamond', CENTERS), ('gross_new_solid', SOLID)]:
        for name, points in groups.items():
            for (x,y),(p,r) in zip(points,convert(name,points)):
                rows.append(dict(formulation=name,kind=kind,pressure_atm=p,rate_cm_s=r,pixel_x=x,pixel_y=y))
    with (HERE/'reference/figure10_digitized.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]); w.writeheader(); w.writerows(rows)
    metadata=dict(source_pdf=str(PDF),source_pdf_sha256=hashlib.sha256(PDF.read_bytes()).hexdigest(),
        source_image_sha256=hashlib.sha256(IMAGE.read_bytes()).hexdigest(),pdf_page=9,printed_page=990,
        image_size=list(Image.open(IMAGE).size),method='Manual centers/solid-line knots in original embedded JPEG; logarithmic axis transform',
        axis_calibration_pixels=AXES,conservative_coordinate_uncertainty_pixels=2,
        rate_graphical_resolution_percent={k:100*(10**(2/(a['y1']-a['y10']))-1) for k,a in AXES.items()},
        notes=['Graphical resolution is not experimental uncertainty.',
               'Caption says M02; panel, Table 2, and text identify M03.',
               'Do not extrapolate M21 solid line beyond its last plotted point.',
               'Use extracted experimental pressures for the sweep; decimals do not imply measurement precision.'])
    (HERE/'reference/figure10_provenance.json').write_text(json.dumps(metadata,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(11.53,10.68))
    ax.imshow(Image.open(IMAGE))
    for name,points in CENTERS.items():
        x,y=np.array(points).T;ax.scatter(x,y,s=150,facecolors='none',edgecolors='#009E73',lw=1)
    ax.axis('off');fig.tight_layout(pad=0);fig.savefig(HERE/'reference/figure10_digitization_audit.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(9,7),layout='constrained')
    for ax,name in zip(axes.flat,AXES):
        p,r=convert(name,CENTERS[name]).T;ax.loglog(p,r,'D',color='k',ms=4,label='Miller, digitized')
        p,r=convert(name,SOLID[name]).T;ax.loglog(p,r,'-',color='.45',label='Gross new, digitized')
        ax.set(xlabel='Pressure (atm)',ylabel='Regression speed (cm/s)',title=name,xlim=(5,250))
        ax.grid(True,which='both',alpha=.2)
    axes[0,0].legend(fontsize=9)
    fig.savefig(HERE/'reference/figure10_reference.png',dpi=200);plt.close(fig)
    print(json.dumps({name:convert(name,CENTERS[name]).round(4).tolist() for name in AXES},indent=2))

if __name__=='__main__':main()
