#!/usr/bin/env python3
"""Fit deepest-front recession and audit sensitivity to startup truncation."""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from analyze_runs import HERE, history, plt


def linear_fit(time, recession, cutoff, minimum_samples=4):
    """Ordinary least squares with a free intercept; coordinates are SI."""
    keep = (time >= cutoff) & np.isfinite(time) & np.isfinite(recession)
    t, depth = time[keep], recession[keep]
    if len(t) < minimum_samples or np.any(np.diff(t) <= 0):
        return None
    slope, intercept = np.polyfit(t - t[0], depth, 1)
    residual = depth - (intercept + slope * (t - t[0]))
    sse = float(residual @ residual)
    sst = float(np.sum((depth - depth.mean()) ** 2))
    return dict(requested_cutoff_s=float(cutoff), start_s=float(t[0]),
                end_s=float(t[-1]), samples=len(t), rate_cm_s=float(slope * 100),
                intercept_recession_m=float(intercept),
                rmse_m=float(np.sqrt(sse / len(t))),
                r_squared=1 - sse / sst if sst > 0 else None,
                measured_recession_m=float(depth[-1] - depth[0]))


def analyze_front(case, discard_fraction):
    rows = history(case, completed_only=True)
    if not rows:
        return None
    meta = json.loads((case / 'case.json').read_text())
    t = np.array([r['time_s'] for r in rows])
    y = np.array([r.get('interface_min_y_m') for r in rows], dtype=float)
    if not np.isfinite(y[0]):
        raise ValueError(f'{case.name}: initial interface position unavailable')
    # The solid is below the gas, so decreasing minimum y means deeper burnback.
    depth = y[0] - y
    h = np.array([r['mean_solid_height_m'] for r in rows])
    mean_depth = h[0] - h
    cutoff = float(t[0] + discard_fraction * (t[-1] - t[0]))
    full = linear_fit(t, depth, t[0])
    truncated = linear_fit(t, depth, cutoff)
    scan = [fit for start in t
            if (fit := linear_fit(t, depth, start)) is not None]
    records = [dict(time_s=r['time_s'], deepest_y_m=r.get('interface_min_y_m'),
                    deepest_recession_m=float(d), mean_recession_m=float(m),
                    multivalued_column_fraction=r['multivalued_column_fraction'],
                    source_case=r['source_case'], plotfile=r['plotfile'])
               for r, d, m in zip(rows, depth, mean_depth)]
    result = dict(case=case.name, formulation=meta['formulation'],
                  pressure_atm=meta['pressure_atm'], latest_saved_time_s=float(t[-1]),
                  initial_deepest_y_m=float(y[0]), snapshot_count=len(rows),
                  discard_fraction_of_elapsed_time=discard_fraction,
                  full_fit=full, truncated_fit=truncated, cutoff_scan=scan,
                  mean_front_full_fit=linear_fit(t, mean_depth, t[0]),
                  mean_front_truncated_fit=linear_fit(t, mean_depth, cutoff),
                  max_multivalued_column_fraction=max(r['multivalued_column_fraction'] for r in rows),
                  interpretation='Diagnostic fits only. Removing startup samples or a high R-squared does not establish statistical stationarity.')
    return result, records


def plot_results(results, output):
    columns = min(2, len(results))
    nrows = (len(results) + columns - 1) // columns
    fig, axes = plt.subplots(nrows, columns, figsize=(6.2 * columns, 4.4 * nrows),
                             squeeze=False, layout='constrained')
    for axis, (result, records) in zip(axes.flat, results):
        t = np.array([r['time_s'] for r in records])
        scale, unit = (1e6, 'µs') if t[-1] < .001 else (1e3, 'ms')
        axis.plot(t * scale, [r['deepest_recession_m'] * 1e6 for r in records],
                  'o', color='#192b41', markersize=4, label='Deepest surface point')
        axis.plot(t * scale, [r['mean_recession_m'] * 1e6 for r in records],
                  '.-', color='#888888', alpha=.8, label='Mean recession (volume)')
        for key, label, color, style in [('full_fit', 'Full fit', '#2268a2', '--'),
                                          ('truncated_fit', 'Truncated fit', '#bc4c15', '-')]:
            fit = result[key]
            if fit is None:
                continue
            tf = np.array([fit['start_s'], fit['end_s']])
            yf = fit['intercept_recession_m'] + fit['rate_cm_s'] / 100 * (tf - tf[0])
            axis.plot(tf * scale, yf * 1e6, style, color=color, linewidth=2,
                      label=f"{label}: {fit['rate_cm_s']:.3f} cm/s (n={fit['samples']})")
        cut = t[0] + result['discard_fraction_of_elapsed_time'] * (t[-1] - t[0])
        axis.axvspan(t[0] * scale, cut * scale, color='#cfcfcf', alpha=.22)
        axis.axvline(cut * scale, color='#999999', linestyle=':', linewidth=1)
        axis.set(title=f"{result['formulation']} · {result['pressure_atm']:.2f} atm",
                 xlabel=f'Physical time ({unit})', ylabel='Recession from initial surface (µm)')
        axis.grid(alpha=.2)
        axis.legend(fontsize=8, loc='upper left')
    for axis in axes.flat[len(results):]:
        axis.set_visible(False)
    percent = results[0][0]['discard_fraction_of_elapsed_time'] * 100
    fig.suptitle(f'Deepest-front least-squares fits · shaded first {percent:g}% excluded from truncated fit\n'
                 'Saved transient data; fitted slopes do not establish steady burning', fontsize=12)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cases', nargs='*', type=Path)
    parser.add_argument('--discard-fraction', type=float, default=.5,
                        help='Fraction of the available elapsed time excluded from the displayed truncated fit.')
    parser.add_argument('--output-dir', type=Path, default=HERE / 'analysis/front_regression')
    args = parser.parse_args()
    if not 0 <= args.discard_fraction < 1:
        parser.error('--discard-fraction must be in [0, 1)')
    cases = args.cases or [HERE.parent / line.strip()
                          for line in (HERE / 'production_cases.txt').read_text().splitlines()
                          if line.strip()]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = []
    for case in cases:
        case = case.resolve()
        if not (case / 'output/celloutput.visit').exists():
            continue
        analyzed = analyze_front(case, args.discard_fraction)
        if analyzed is None:
            continue
        result, records = analyzed
        results.append(analyzed)
        with (output / f'{case.name}_history.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=records[0])
            writer.writeheader()
            writer.writerows(records)
    if not results:
        parser.error('No completed field history available.')
    report = dict(computed_utc=datetime.now(timezone.utc).isoformat(),
                  metric='min_x y at eta=0.5; recession = initial minimum y minus current minimum y',
                  fit_method='Unweighted ordinary least squares, degree 1, free intercept; each saved time has equal weight.',
                  cutoff_policy='Display a prescribed elapsed-time cutoff and record every cutoff retaining at least four samples; no cutoff chosen by agreement with reference data.',
                  limitations='The extremum can change horizontal location. It describes the deepest front, not bulk regression. Linear slope is an interval average, not an instantaneous derivative.',
                  results=[result for result, records in results])
    (output / 'fits.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    plot_results(results, output / 'front_regression.png')
    print(json.dumps(report, indent=2))
    print(f"Plot: {output / 'front_regression.png'}")


if __name__ == '__main__':
    main()
