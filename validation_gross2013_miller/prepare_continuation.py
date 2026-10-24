#!/usr/bin/env python3
"""Prepare a distinct continuation from a completed, explicitly selected checkpoint."""
import argparse
import hashlib
import json
from pathlib import Path
from prepare_compact import ROOT, HERE, SPECIES, deck_read, save_case


def prepare(checkpoint, name, stop_s, plot_dt_s, max_steps):
    checkpoint = checkpoint.resolve()
    parent = checkpoint.parent.parent
    if parent.parent != HERE / 'runs' or checkpoint.parent.name != 'output':
        raise ValueError('Checkpoint must belong to this study')
    receipt = json.loads((parent / 'run.json').read_text())
    if receipt.get('returncode') != 0 or not receipt.get('amrex_finalized'):
        raise ValueError('Use a finalized parent run, not a possibly incomplete write')
    header = (checkpoint / 'Header').read_text().splitlines()
    start_s = float(header[int(header[1]) + 3])
    start_step = int(header[int(header[1]) + 9].split()[0])
    if not stop_s > start_s or plot_dt_s <= 0 or max_steps < 1:
        raise ValueError('Invalid continuation interval')
    if Path(name).name != name:
        raise ValueError('Case name must be a basename')
    case = HERE / 'runs' / name
    if case.exists():
        raise FileExistsError(case)
    meta = json.loads((parent / 'case.json').read_text())
    deck = deck_read(parent / 'input')
    for key in list(deck):
        if key.startswith(('temperature.ic.', 'velocity.ic.')) or any(
                key.startswith(s + '.density.ic.') for s in SPECIES):
            del deck[key]
    deck.update({'temperature.ic.type': 'constant',
                 'temperature.ic.constant.value': '300',
                 'velocity.ic.type': 'constant',
                 'velocity.ic.constant.value': '0 0'})
    for species in SPECIES:
        deck[species + '.density.ic.type'] = 'constant'
        deck[species + '.density.ic.constant.value'] = '0'
    deck.update({'restart_cell': str(checkpoint.relative_to(ROOT)),
                 'restart.in_place': '0', 'stop_time': f'{stop_s:.17g}',
                 'max_step': str(start_step + max_steps), 'amr.plot_dt': f'{plot_dt_s:.17g}',
                 'plot_file': str(case.relative_to(ROOT) / 'output'),
                 'run_control.stop_file': str(case.relative_to(ROOT) / 'STOP'),
                 'run_control.poll_interval': '10', 'diagnostics.interval': '50'})
    meta.update(pilot=True, initialization_only=False,
                restart_from=str(checkpoint.relative_to(ROOT)),
                checkpoint_header_sha256=hashlib.sha256((checkpoint / 'Header').read_bytes()).hexdigest(),
                parent_binary_sha256=receipt['binary_sha256'],
                requested_stop_s=stop_s, plot_dt_s=plot_dt_s,
                continuation_start_s=start_s, continuation_start_step=start_step,
                maximum_additional_steps=max_steps, absolute_max_step=start_step+max_steps,
                check_scope='Bounded startup relaxation and cost diagnostic; not an accepted regression rate.',
                rate_status='not_run')
    save_case(case, deck, meta)
    return case


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path)
    parser.add_argument('name')
    parser.add_argument('--stop', type=float, required=True)
    parser.add_argument('--plot-dt', type=float, required=True)
    parser.add_argument('--max-steps', type=int, default=1000)
    args = parser.parse_args()
    print(prepare(args.checkpoint, args.name, args.stop, args.plot_dt, args.max_steps).relative_to(ROOT))
