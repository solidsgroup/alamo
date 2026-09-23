#!/usr/bin/env python3
"""Load-depth curve and checks for this input's loading schedule (3D only)."""
import csv
import json
import re
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yt

root = Path(sys.argv[1])
log = (root/'out.log').read_text()
assert 'finalized' in log and 'SIGABRT' not in log, 'Run did not finish successfully'
matches = re.findall(r'Contact time=(\S+) stable=1 active_nodes=(\d+) load=(\S+) penetration=(\S+) tension=(\S+)', log)
states = {float(t): (int(n),float(p),float(g),float(q)) for t,n,p,g,q in matches}
assert states, 'No converged contact results'
assert max(states) >= .25-1e-10, 'The loading/unloading schedule is incomplete'
rows = []
for t, (n, load, penetration, tension) in sorted(states.items()):
    depth = .5*t if t < .1 else (.05 if t < .15 else .05-.6*(t-.15))
    rows.append((t,depth*1000,load,n,penetration,tension))
with (root/'load_depth.csv').open('w') as stream:
    writer = csv.writer(stream)
    writer.writerow(['solve_time_s','tip_depth_nm','load_uN','active_nodes','penetration_um','tensile_traction_MPa'])
    writer.writerows(rows)
fig, ax = plt.subplots(figsize=(5,4), layout='constrained')
ax.plot([r[1] for r in rows], [r[2] for r in rows], '.-')
ax.set(xlabel='Tip depth (nm)', ylabel='Compressive load (µN)', title='Small-slope spherical indentation')
fig.savefig(root/'load_depth.png', dpi=180)
yt.set_log_level(50)
plots = sorted(root.glob('*cell'))
assert plots
ds = yt.load(str(plots[-1]))
assert ds.dimensionality == 3
data = ds.all_data()
slips = [name for kind,name in ds.field_list if kind == 'boxlib' and 'gamma' in name]
max_slip = max((float(np.max(np.abs(data['boxlib',name].d))) for name in slips), default=0.)
level = int(ds.index.max_level)
dims = ds.domain_dimensions * ds.refine_by**level
grid = ds.covering_grid(level, ds.domain_left_edge, dims)
section = np.zeros((int(dims[0]),int(dims[2])))
for name in slips:
    section += np.abs(grid['boxlib',name].d[:,int(dims[1])//2,:])
    grid.clear_data()
fig, ax = plt.subplots(figsize=(6,3.5), layout='constrained')
im = ax.imshow(section.T, origin='lower', aspect='equal',
               extent=[float(ds.domain_left_edge[0]),float(ds.domain_right_edge[0]),
                       float(ds.domain_left_edge[2]),float(ds.domain_right_edge[2])])
ax.set(xlabel='x (µm)', ylabel='z (µm)', title='Final central section: sum of absolute signed slips')
fig.colorbar(im, ax=ax, label='Σ |γᵅ| (not accumulated absolute slip)')
fig.savefig(root/'slip_section.png', dpi=180)
report = dict(grid=ds.domain_dimensions.tolist(), solved_times=len(rows),
              finest_level=level, finest_equivalent_grid=dims.tolist(),
              stored_cells=sum(int(np.prod(g.ActiveDimensions)) for g in ds.index.grids),
              last_solve_time_s=rows[-1][0], peak_load_uN=max(r[2] for r in rows),
              final_load_uN=rows[-1][2], final_active_nodes=rows[-1][3],
              max_penetration_um=max(r[4] for r in rows),
              max_tensile_traction_MPa=max(r[5] for r in rows),
              final_max_abs_signed_slip=max_slip,
              mesh_converged=False)
(root/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
