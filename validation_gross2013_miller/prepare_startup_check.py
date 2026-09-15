#!/usr/bin/env python3
"""Prepare, but do not launch, a 100-step check from an audited t=0 field."""
import argparse
import hashlib
import json
from pathlib import Path
from prepare_compact import ROOT, HERE, deck_read, save_case, SPECIES

def prepare(preview):
    meta=json.loads((preview/'case.json').read_text())
    receipt=json.loads((preview/'run.json').read_text())
    assert receipt['returncode']==0 and receipt['amrex_finalized']
    assert meta.get('initialization_only')
    checkpoint=preview/'output/00000cell'
    assert (checkpoint/'Header').exists()
    base=ROOT/meta['production_case'];case=base.with_name(base.name+'_startup100')
    d=deck_read(preview/'input')
    # The saved field, not a second expression evaluation, supplies the entire
    # initial state. Keep the solid geometry IC for subsequent AMR refinement.
    for key in list(d):
        if key.startswith(('temperature.ic.','velocity.ic.')) or any(key.startswith(s+'.density.ic.') for s in SPECIES):del d[key]
    d['temperature.ic.type']='constant';d['temperature.ic.constant.value']='300'
    d['velocity.ic.type']='constant';d['velocity.ic.constant.value']='0 0'
    for s in SPECIES:d[s+'.density.ic.type']='constant';d[s+'.density.ic.constant.value']='0'
    d['restart_cell']=str(checkpoint.relative_to(ROOT));d['restart.in_place']='0'
    d['max_step']='100';d['plot_file']=str(case.relative_to(ROOT)/'output')
    d['run_control.stop_file']=str(case.relative_to(ROOT)/'STOP');d['diagnostics.interval']='10'
    meta.update(pilot=True,initialization_only=False,restart_from=str(checkpoint.relative_to(ROOT)),
        checkpoint_header_sha256=hashlib.sha256((checkpoint/'Header').read_bytes()).hexdigest(),
        check_scope='100-step startup/cost diagnostic; never an accepted regression rate.',requested_stop_s=.0005)
    save_case(case,d,meta)
    return case

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('preview',type=Path);a=p.parse_args()
    print(prepare(a.preview.resolve()).relative_to(ROOT))
