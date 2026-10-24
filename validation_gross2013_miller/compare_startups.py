#!/usr/bin/env python3
"""Report actual startup fields/cost; never interpret them as steady rates."""
import json
from pathlib import Path
import numpy as np
from analyze_runs import snapshot

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent

def main():
    records=[]
    for line in (HERE/'preview_cases.txt').read_text().splitlines():
        preview=ROOT/line;m=json.loads((preview/'case.json').read_text())
        old=HERE/'runs'/m['predecessor_case']
        new=ROOT/m['production_case'];new=new.with_name(new.name+'_startup100')
        for label,case in [('planar_full_height',old),('local_flame_compact',new)]:
            receipt=case/'run.json';final=case/'output/00100cell'
            if not receipt.exists() or not (final/'Header').exists():continue
            r=json.loads(receipt.read_text())
            if r.get('returncode')!=0:continue
            cm=json.loads((case/'case.json').read_text())
            initial=ROOT/cm['restart_from'] if cm.get('restart_from') else case/'output/00000cell'
            # Restart's initial diagnostic write precedes reference-pressure
            # synchronization and can contain qdot=0. The checkpoint holds
            # the same physical state with initialized chemistry diagnostics.
            first=snapshot(initial);last=snapshot(final)
            dt=last['time_s']-first['time_s']
            records.append(dict(formulation=m['formulation'],initialization=label,case=case.name,
                initial_diagnostic_plotfile=str(initial.relative_to(ROOT)),
                physical_time_s=dt,mean_dt_s=dt/100,wall_seconds=r['wall_seconds'],
                wall_seconds_per_simulated_us=r['wall_seconds']/(dt*1e6),
                initial_temperature_max_K=first['temperature_max_K'],final_temperature_max_K=last['temperature_max_K'],
                initial_velocity_max_m_s=first['velocity_max_m_s'],final_velocity_max_m_s=last['velocity_max_m_s'],
                initial_gas_heat_MW_m2=first['gas_heat_release_W_m2']/1e6,
                final_gas_heat_MW_m2=last['gas_heat_release_W_m2']/1e6,
                initial_mean_surface_temperature_K=first['surface_temperature_K'],
                final_mean_surface_temperature_K=last['surface_temperature_K'],
                recession_um=(first['mean_solid_height_m']-last['mean_solid_height_m'])*1e6))
    answer=dict(check='100-step startup diagnostics only; no statistically steady rates.',
        cost_note='Timings include initialization/output and may reflect different concurrent loads; not a controlled speedup benchmark.',records=records)
    (HERE/'analysis/startup_comparison.json').write_text(json.dumps(answer,indent=2)+'\n')
    print(json.dumps(answer,indent=2))

if __name__=='__main__':main()
