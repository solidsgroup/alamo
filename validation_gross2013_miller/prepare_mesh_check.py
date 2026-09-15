#!/usr/bin/env python3
"""Prepare a cheaper mesh pilot with the identical physical initial expressions."""
import argparse
import json
from pathlib import Path
from prepare_compact import ROOT,HERE,deck_read,save_case


def prepare(source):
    source=source.resolve();d=deck_read(source/'input');m=json.loads((source/'case.json').read_text())
    old_level=int(d['amr.max_level']);new_level=old_level-1
    if new_level<1:raise ValueError('Insufficient AMR levels for this check')
    case=source.with_name(source.name.replace('_compact_diffusion',f'_mesh{new_level}_relax2us'))
    if case.exists():raise FileExistsError(case)
    d.update({'amr.max_level':str(new_level),'max_step':'1000','stop_time':'2e-6',
              'amr.plot_dt':'2.5e-7','plot_file':str(case.relative_to(ROOT)/'output'),
              'run_control.stop_file':str(case.relative_to(ROOT)/'STOP'),
              'run_control.poll_interval':'10','diagnostics.interval':'50'})
    d.pop('restart_cell',None);d.pop('restart.in_place',None)
    reference_dx=m['dx_m'];m['dx_m']=reference_dx*2
    m.update(pilot=True,initialization_only=False,requested_stop_s=2e-6,plot_dt_s=2.5e-7,
             rate_status='not_run',check_scope='Startup discretization/cost check only; cannot establish a steady regression rate.',
             mesh_check=dict(reference_case=str(source.relative_to(ROOT)),reference_dx_m=reference_dx,
                             candidate_dx_m=m['dx_m'],physical_interface_width_unchanged=True,
                             cells_per_interface_width=m['interface_width_m']/m['dx_m'],
                             initial_fields='Same physical temperature, velocity, species and packing definitions; sampled on a coarser mesh.',
                             acceptance='Requires comparison to the fine mesh. No automatic production-mesh change.'))
    m['shortened_domain'].pop('preserved_finest_dx_m',None)
    m['shortened_domain']['candidate_finest_dx_m']=m['dx_m']
    save_case(case,d,m)
    return case


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('source',type=Path);args=parser.parse_args()
    print(prepare(args.source).relative_to(ROOT))
