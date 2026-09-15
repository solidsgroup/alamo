#!/usr/bin/env python3
"""Diagnose startup cost, flame evolution, and boundary proximity from saved fields."""
import argparse
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/alamo-gross-matplotlib')
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import yt
from analyze_runs import HERE, history
from prepare_compact import deck_read
from rocfire_diagnostics import stored_or_reconstructed_heat


def field_diagnostics(case, plotfile):
    meta = json.loads((case / 'case.json').read_text())
    deck = deck_read(case / 'input')
    ds = yt.load(str(case / 'output' / plotfile)); ad = ds.all_data()
    area = ad['index', 'dx'].d * ad['index', 'dy'].d
    x = ad['index', 'x'].d; y = ad['index', 'y'].d
    eta = ad['boxlib', 'rigid_eta'].d
    temp = ad['boxlib', 'temperature'].d
    speed = np.hypot(ad['boxlib', 'velocityx'].d, ad['boxlib', 'velocityy'].d)
    gas = sum(ad['boxlib', 'component_density_' + s].d for s in
              ('AP_gas', 'HTPB_gas', 'Mono', 'Premixed', 'Primary', 'Final'))
    residual = eta + gas * (8.31446261815324 / .026) * temp / meta['pressure_Pa'] - 1
    peak = int(np.argmax(speed / np.minimum(ad['index', 'dx'].d, ad['index', 'dy'].d)))
    valid_rate = float(np.max(speed / np.minimum(ad['index', 'dx'].d, ad['index', 'dy'].d)))
    # Match the solver's CFL reduction, including covered coarse cells.
    all_rate = 0.
    for grid in ds.index.grids:
        velocity = np.hypot(grid['boxlib', 'velocityx'].d, grid['boxlib', 'velocityy'].d)
        all_rate = max(all_rate, float(velocity.max() / min(grid.dds.d[:2])))
    q = np.maximum(stored_or_reconstructed_heat(ad, meta), 0)
    top_y = float(ds.domain_right_edge[1]); bottom_y = float(ds.domain_left_edge[1])
    top_band = y > .9 * top_y
    bottom_band = y < bottom_y + .1 * abs(bottom_y)
    answer = dict(valid_cells=len(area), cfl=float(deck['cfl']),
                  predicted_advection_dt_s=float(deck['cfl']) / max(all_rate, 1e-100),
                  covered_to_valid_CFL_rate_ratio=all_rate / max(valid_rate, 1e-100),
                  limiting_valid_cell=dict(x_um=float(x[peak]*1e6), y_um=float(y[peak]*1e6),
                                          temperature_K=float(temp[peak]), eta=float(eta[peak]),
                                          speed_m_s=float(speed[peak])),
                  maximum_volume_residual=float(np.max(np.abs(residual))),
                  rms_volume_residual=float(np.sqrt(np.sum(residual**2 * area) / area.sum())),
                  heat_fraction_in_top_tenth_of_initial_gas=float(np.sum(q[top_band]*area[top_band]) / max(np.sum(q*area), 1e-100)),
                  bottom_tenth_temperature_max_K=float(temp[bottom_band].max()),
                  boundary_note='Screening diagnostics only; domain independence requires an enlarged-domain comparison.')
    return answer


def report(cases):
    records=[]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout='constrained')
    for case in cases:
        rows=history(case)
        if not rows:continue
        receipt=json.loads((case/'run.json').read_text())
        meta=json.loads((case/'case.json').read_text())
        final=rows[-1]
        # Restart's first derived qdot may be zero before pressure synchronization.
        # The selected checkpoint has the identical physical state and valid qdot.
        diagnostic_path=HERE/'runs'/final['gas_heat_diagnostic_source']
        diagnostic=field_diagnostics(diagnostic_path.parent, diagnostic_path.name)
        diagnostic['diagnostic_plotfile']=str((diagnostic_path.parent/'output'/diagnostic_path.name).relative_to(HERE.parent))
        t=np.array([r['time_s'] for r in rows]); h=np.array([r['mean_solid_height_m'] for r in rows])
        entry=dict(case=case.name, formulation=meta['formulation'], receipt_status=receipt['status'],
                   accepted_steady_rate=False, saved_time_s=float(t[-1]),
                   recession_um=float((h[0]-h[-1])*1e6),
                   final_surface_temperature_K=final['surface_temperature_K'],
                   final_gas_heat_MW_m2=final['gas_heat_release_W_m2']/1e6,
                   final_temperature_max_K=final['temperature_max_K'], **diagnostic)
        if receipt.get('wall_seconds') and meta.get('continuation_start_s') is not None:
            interval=t[-1]-meta['continuation_start_s']
            if interval>0:entry['wall_seconds_per_simulated_us']=receipt['wall_seconds']/(interval*1e6)
        records.append(entry)
        axes[0,0].plot(t*1e6, (h[0]-h)*1e6, '.-', label=meta['formulation'])
        axes[0,1].plot(t*1e6, [r['surface_temperature_K'] for r in rows], '.-')
        axes[1,0].plot(t*1e6, [r['gas_heat_release_W_m2']/1e6 for r in rows], '.-')
        axes[1,1].plot(t*1e6, [r['velocity_max_m_s'] for r in rows], '.-')
    for ax,label in zip(axes.flat, ('Mean recession (µm)', 'Mean surface temperature (K)',
                                   'Gas heat release (MW/m²)', 'Maximum velocity (m/s)')):
        ax.set(xlabel='Physical time (µs)', ylabel=label);ax.grid(alpha=.2)
    if records:axes[0,0].legend()
    fig.suptitle('Startup relaxation diagnostics — no accepted regression rates')
    fig.savefig(HERE/'analysis/startup_relaxation.png', dpi=180);plt.close(fig)
    answer=dict(note='Short startup records; not Figure 10 validation rates.', records=records)
    (HERE/'analysis/startup_relaxation.json').write_text(json.dumps(answer,indent=2)+'\n')
    return answer


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cases',type=Path,nargs='+');args=parser.parse_args()
    print(json.dumps(report([p.resolve() for p in args.cases]),indent=2))
