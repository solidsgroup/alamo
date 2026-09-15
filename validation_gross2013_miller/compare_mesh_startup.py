#!/usr/bin/env python3
"""Compare common-time startup histories; these are not steady rate errors."""
import argparse
import json
from pathlib import Path
import numpy as np
from analyze_runs import HERE,history


def compare(fine,coarse):
    histories=[history(case) for case in (fine,coarse)]
    end=min(rows[-1]['time_s'] for rows in histories)
    if end<=0:raise ValueError('Both cases must advance')
    rows=[]
    for case,values in zip((fine,coarse),histories):
        meta=json.loads((case/'case.json').read_text());receipt=json.loads((case/'run.json').read_text())
        times=np.array([v['time_s'] for v in values]);height=np.array([v['mean_solid_height_m'] for v in values])
        row=dict(case=case.name,dx_um=meta['dx_m']*1e6,receipt_status=receipt['status'],
                 comparison_time_us=end*1e6,
                 recession_um=float((height[0]-np.interp(end,times,height))*1e6))
        for key in ('surface_temperature_K','temperature_max_K','gas_heat_release_W_m2'):
            row[key]=float(np.interp(end,times,[v[key] for v in values]))
        outputs=sorted((case/'output').glob('*cell'))
        if receipt.get('returncode')==0 and len(outputs)>1:
            # Filesystem timestamps bracket advancement, excluding costly IC
            # expression parsing; include output cost and concurrent-load effects.
            elapsed=(outputs[-1]/'Header').stat().st_mtime-(outputs[0]/'Header').stat().st_mtime
            interval=values[-1]['time_s']-meta.get('continuation_start_s',0)
            row.update(wall_s_after_initial_header=elapsed,
                       approximate_wall_s_per_simulated_us=elapsed/(interval*1e6),
                       total_receipt_wall_s=receipt['wall_seconds'])
        rows.append(row)
    a,b=rows
    return dict(formulation=json.loads((fine/'case.json').read_text())['formulation'],records=rows,
                recession_difference_percent=100*(b['recession_um']/a['recession_um']-1),
                surface_temperature_difference_K=b['surface_temperature_K']-a['surface_temperature_K'],
                gas_heat_difference_percent=100*(b['gas_heat_release_W_m2']/a['gas_heat_release_W_m2']-1),
                accepted_production_mesh=False,
                note='Common-time interpolation of short startup histories. These differences are not errors versus Miller or evidence of steady mesh convergence. Timing is approximate and load-dependent.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fine',type=Path);parser.add_argument('coarse',type=Path);args=parser.parse_args()
    result=compare(args.fine.resolve(),args.coarse.resolve())
    (HERE/'analysis'/(result['formulation']+'_mesh_startup_comparison.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
