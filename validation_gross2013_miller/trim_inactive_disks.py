#!/usr/bin/env python3
"""Cull disk entries whose compact PSRead support cannot touch a case's domain.

Original line text and ordering are preserved. A four-base-cell halo exceeds
the three ghost cells used by this study's MUSCL setup. Periodic x copies are
unaffected; y must be nonperiodic. This changes geometry-search cost only.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from prepare_compact import ROOT, HERE, deck_read, save_case


def trim(case):
    case=case.resolve()
    if (case/'run.json').exists():raise ValueError('Never alter a launched input')
    deck=deck_read(case/'input');meta=json.loads((case/'case.json').read_text())
    if meta.get('packing_cull'):return meta['packing_cull']
    if deck['geometry.is_periodic'].split()[1]!='0' or deck['advection.type']!='muscl':
        raise ValueError('This culling bound is for nonperiodic-y MUSCL cases')
    keys=[s+'.density.ic.psread.' for s in ('AP_solid','HTPB_solid')]
    source=Path(deck[keys[0]+'file.name'].strip('"'))
    if not source.is_absolute():source=ROOT/source
    for key in keys:
        if deck[key+'file.unit']!='mm' or Path(deck[key+'file.name'].strip('"')).name!=source.name:
            raise ValueError('Expected shared millimetre disk file')
        if any(float(x.strip('"'))!=0 for x in deck.get(key+'x0','0 0 0').split()):
            raise ValueError('Translated disk input needs a transformed culling bound')
        if float(deck.get(key+'mult','1').strip('"'))!=1:
            raise ValueError('Scaled disk input needs a transformed culling bound')
    lines=source.read_text().splitlines(keepends=True)
    disks=np.array([[float(x) for x in line.split()] for line in lines])*.001
    lo=float(deck['geometry.prob_lo'].split()[1]);hi=float(deck['geometry.prob_hi'].split()[1])
    dy=(hi-lo)/int(deck['amr.n_cell'].split()[1]);halo=4*dy
    eps=max(float(deck[key+'eps']) for key in keys)
    # PSRead has exactly zero disk support outside radius + eps.
    clearance=np.maximum(lo-halo-(disks[:,1]+disks[:,3]+eps),
                         (disks[:,1]-disks[:,3]-eps)-(hi+halo))
    keep=clearance<=0
    out=HERE/'packs_domain';out.mkdir(exist_ok=True)
    target=out/(case.name+'.xyzr')
    text=''.join(line for line,retained in zip(lines,keep) if retained)
    target.write_text(text)
    proof=dict(original_file=str(source.relative_to(ROOT)),
               original_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
               original_entries=len(disks),retained_entries=int(keep.sum()),
               omitted_entries=int((~keep).sum()),halo_m=halo,
               minimum_omitted_support_clearance_m=float(clearance[~keep].min()) if (~keep).any() else None,
               retained_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
               note='Only entries with zero support throughout the domain and ghost halo were removed; retained lines/order and material parameters are unchanged.')
    for key in keys:deck[key+'file.name']=str(target.relative_to(ROOT))
    meta.update(packing_file=str(target.relative_to(ROOT)),packing_cull=proof)
    save_case(case,deck,meta)
    return proof


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cases',type=Path,nargs='+');args=parser.parse_args()
    for case in args.cases:print(json.dumps(dict(case=case.name,**trim(case))))
