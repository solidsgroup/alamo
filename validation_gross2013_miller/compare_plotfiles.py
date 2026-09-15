#!/usr/bin/env python3
"""Compare valid AMR cell fields, independent of MPI output ordering."""
import argparse
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/alamo-gross-matplotlib')
import numpy as np
import yt
yt.set_log_level(50)

def ordered(path):
    ds = yt.load(str(path))
    ad = ds.all_data()
    coords = np.column_stack([ad['index', key].d for key in ('dx', 'x', 'y')])
    order = np.lexsort(coords.T[::-1])
    return ds, ad, coords[order], order

def compare(reference, candidate):
    a, aa, ac, ai = ordered(reference)
    b, ba, bc, bi = ordered(candidate)
    if ac.shape != bc.shape or not np.allclose(ac, bc, rtol=1e-13, atol=1e-18):
        raise ValueError('Comparisons require the same valid AMR cells')
    if not np.isclose(float(a.current_time), float(b.current_time), rtol=1e-12, atol=1e-15):
        raise ValueError('Plot times differ')
    result = dict(reference=str(reference), candidate=str(candidate),
                  time_s=float(a.current_time), valid_cells=len(ai), fields={})
    for field in sorted(set(a.field_list) & set(b.field_list)):
        if field[0] != 'boxlib':
            continue
        av = aa[field].d[ai]; bv = ba[field].d[bi]
        if not np.all(np.isfinite(av)) or not np.all(np.isfinite(bv)):
            raise ValueError(f'Nonfinite field: {field}')
        delta = bv-av
        result['fields'][field[1]] = dict(
            max_absolute_difference=float(np.max(np.abs(delta))),
            relative_L2=float(np.linalg.norm(delta)/max(np.linalg.norm(av), 1e-300)))
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidates', nargs='+', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    results = [compare(args.reference, candidate) for candidate in args.candidates]
    args.output.write_text(json.dumps(results, indent=2)+'\n')
    for r in results:
        print(r['candidate'])
        for name in ('temperature', 'rigid_eta', 'pressure', 'velocityx', 'velocityy'):
            print(name, r['fields'][name])

if __name__ == '__main__':
    main()
