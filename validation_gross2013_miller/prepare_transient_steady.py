#!/usr/bin/env python3
"""Prepare transient cases whose statistical monitor requests their early stop."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from prepare_compact import HERE, ROOT, SPECIES, deck_read, save_case

RESTARTS = {
    'M03': 'M03_p34.3921_seed101_compact_diffusion_relax2us_v2/output/00307cell',
    'M17': 'M17_p34.3921_seed101_compact_diffusion_relax2us/output/00222cell',
    'M21': 'M21_p34.0775_seed101_compact_diffusion_relax2us/output/00225cell',
    'M24': 'M24_p33.6723_seed101_compact_diffusion_relax2us/output/00235cell',
}


def main():
    listing=HERE/'production_cases.txt'
    backup=HERE/'production_cases_before_transient_steady.txt'
    if backup.exists():
        raise RuntimeError('Transient sweep has already been prepared')
    originals=[ROOT/p for p in listing.read_text().splitlines() if p]
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    binary=digest(ROOT/'bin/lowmach-2d-clang++')
    if binary!='f2bc201ab63641fc1e2c95cb4076c9401f16dcb35c253075484c9a9be34dce6c':
        raise RuntimeError('Expected the restored pre-feature transient executable')
    cases=[];groups={};prepared=[]
    for original in originals:
        meta=json.loads((original/'case.json').read_text());deck=deck_read(original/'input')
        case=original.with_name(original.name+'_transient_steady')
        if case.exists():raise FileExistsError(case)
        if digest(original/'input')!=meta['input_sha256']:raise RuntimeError('Input hash mismatch')
        if any('quasi_steady' in k for k in deck):raise RuntimeError('Transient inputs required')
        # Reference data set only output cadence and a ceiling, never a rate or
        # the acceptance decision. About forty samples per particle transit
        # avoid thousands of unnecessarily large plotfiles for the coarse packs.
        particle_time=max(float(d) for d in meta['packing']['bins'])*1e-6/(meta['reference_rate_cm_s']*.01)
        cadence=max(25e-6,particle_time/40)
        ceiling=max(meta['duration_policy']['maximum_reference_duration_s'],.00025+60*cadence)
        if 30<meta['pressure_atm']<40:
            checkpoint=HERE/'runs'/RESTARTS[meta['formulation']]
            parent=checkpoint.parent.parent
            receipt=json.loads((parent/'run.json').read_text())
            parent_meta=json.loads((parent/'case.json').read_text())
            assert receipt['returncode']==0 and receipt['amrex_finalized']
            assert parent_meta['pressure_Pa']==meta['pressure_Pa']
            assert parent_meta['material_properties']==meta['material_properties']
            parent_deck=deck_read(parent/'input')
            parent_pack=ROOT/parent_deck['AP_solid.density.ic.psread.file.name']
            assert digest(parent_pack) in (meta['packing_cull']['original_sha256'],meta['packing_cull']['retained_sha256'])
            for key in ('geometry.prob_lo','geometry.prob_hi','amr.n_cell','amr.max_level',
                        'AP_solid.density.ic.psread.eps','HTPB_solid.density.ic.psread.eps'):
                assert parent_deck[key]==deck[key]
            for key in list(deck):
                if key.startswith(('temperature.ic.','velocity.ic.')) or any(key.startswith(s+'.density.ic.') for s in SPECIES):
                    del deck[key]
            deck.update({'temperature.ic.type':'constant','temperature.ic.constant.value':'300',
                         'velocity.ic.type':'constant','velocity.ic.constant.value':'0 0'})
            for s in SPECIES:
                deck[s+'.density.ic.type']='constant';deck[s+'.density.ic.constant.value']='0'
            deck.update({'restart_cell':str(checkpoint.relative_to(ROOT)),'restart.in_place':'0'})
            header=(checkpoint/'Header').read_text().splitlines()
            meta.update(restart_from=str(checkpoint.relative_to(ROOT)),
                        checkpoint_header_sha256=digest(checkpoint/'Header'),
                        continuation_start_s=float(header[int(header[1])+3]))
        deck.update({'plot_file':str(case.relative_to(ROOT)/'output'),
                     'run_control.stop_file':str(case.relative_to(ROOT)/'STOP'),
                     'run_control.poll_interval':'100','max_step':'2147483647',
                     'stop_time':f'{ceiling:.17g}','amr.plot_dt':f'{cadence:.17g}'})
        meta.update(pilot=False,requested_stop_s=ceiling,plot_dt_s=cadence,rate_status='not_run',
                    transient_steady_stop=True,restored_input_source=str(original.relative_to(ROOT)),
                    stop_control='Root monitor writes STOP after confirmed statistical acceptance; elapsed physical and wall times recorded.',
                    duration_policy=dict(meta['duration_policy'],mode='stop_on_confirmed_statistical_stationarity',
                                         safety_ceiling_s=ceiling,output_interval_s=cadence,
                                         reference_use='Scheduling only; never rate acceptance.'))
        prepared.append((case,deck,meta));cases.append(case)
        groups.setdefault(meta['formulation'],[]).append((meta['pressure_atm'],case))
    for case,deck,meta in prepared:save_case(case,deck,meta)
    for rows in groups.values():rows.sort()
    (HERE/'queues').mkdir(exist_ok=True)
    for label,forms in [('small',('M03','M17')),('large',('M21','M24'))]:
        ordered=[groups[f][i][1] for i in (1,0,2) for f in forms]
        manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),binary_sha256=binary,mpi_ranks=4,
                      monitor_heartbeat_max_age_s=300,
                      cases=[dict(case=str(c.relative_to(ROOT)),input_sha256=digest(c/'input')) for c in ordered])
        (HERE/'queues'/f'fig10_transient_{label}.json').write_text(json.dumps(manifest,indent=2)+'\n')
    backup.write_bytes(listing.read_bytes())
    listing.write_text(''.join(str(c.relative_to(ROOT))+'\n' for c in cases))
    print(f'Prepared {len(cases)} transient cases with statistical early stopping.')


if __name__=='__main__':main()
