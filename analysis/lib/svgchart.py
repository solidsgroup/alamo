#!/usr/bin/env python3
"""Zero-dependency SVG chart helpers (no matplotlib required).

We deliberately avoid a plotting dependency: the suite must run on a bare HPC
login node.  These produce clean, self-contained SVGs that embed directly in the
HTML report.

Exposed:
    grouped_bars(path, title, categories, series, ylabel, log=False)
    timeline(path, title, x, series, xlabel, ylabel)
"""
from __future__ import annotations
import html
import math

_PALETTE = ["#2563eb", "#dc2626", "#16a34a", "#d97706", "#7c3aed", "#0891b2"]


def _fmt(v: float) -> str:
    if v == 0:
        return "0"
    a = abs(v)
    if a >= 1000 or a < 0.01:
        return f"{v:.2e}"
    if a >= 1:
        return f"{v:.2f}"
    return f"{v:.3f}"


def grouped_bars(path, title, categories, series, ylabel="", log=False):
    """series = {name: [values aligned with categories]}."""
    W, H = 900, 460
    ml, mr, mt, mb = 70, 180, 56, 70
    pw, ph = W - ml - mr, H - mt - mb
    names = list(series.keys())
    allvals = [v for vs in series.values() for v in vs if v is not None]
    vmax = max(allvals) if allvals else 1.0
    if log:
        vmax = math.log10(max(vmax, 1e-9)) + 0.05
        def sy(v): return 0 if v in (None, 0) else max(0, math.log10(max(v, 1e-9)))
    else:
        vmax = vmax * 1.12 or 1.0
        def sy(v): return v or 0.0

    def Y(v): return mt + ph - (sy(v) / vmax) * ph

    s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'font-family="Segoe UI,Helvetica,Arial,sans-serif" font-size="13">']
    s.append(f'<rect width="{W}" height="{H}" fill="white"/>')
    s.append(f'<text x="{W/2}" y="28" text-anchor="middle" font-size="17" '
             f'font-weight="600">{html.escape(title)}</text>')
    # gridlines
    for i in range(6):
        gy = mt + ph - i / 5 * ph
        s.append(f'<line x1="{ml}" y1="{gy:.1f}" x2="{ml+pw}" y2="{gy:.1f}" '
                 f'stroke="#e5e7eb"/>')
        gv = (i / 5 * vmax)
        gv = 10 ** gv if log else gv
        s.append(f'<text x="{ml-8}" y="{gy+4:.1f}" text-anchor="end" '
                 f'fill="#6b7280">{_fmt(gv)}</text>')
    if ylabel:
        s.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" fill="#374151" '
                 f'transform="rotate(-90 16 {mt+ph/2})">{html.escape(ylabel)}</text>')

    ncat = len(categories)
    gw = pw / max(ncat, 1)
    nb = len(names)
    bw = gw * 0.7 / max(nb, 1)
    for ci, cat in enumerate(categories):
        gx = ml + ci * gw
        for bi, name in enumerate(names):
            v = series[name][ci] if ci < len(series[name]) else None
            x = gx + gw * 0.15 + bi * bw
            y = Y(v)
            h = mt + ph - y
            col = _PALETTE[bi % len(_PALETTE)]
            s.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bw*0.92:.1f}" '
                     f'height="{max(h,0):.1f}" fill="{col}" rx="2"/>')
            if v is not None:
                s.append(f'<text x="{x+bw*0.46:.1f}" y="{y-4:.1f}" '
                         f'text-anchor="middle" font-size="11" fill="#374151">'
                         f'{_fmt(v)}</text>')
        s.append(f'<text x="{gx+gw/2:.1f}" y="{mt+ph+22}" text-anchor="middle" '
                 f'font-weight="600">{html.escape(str(cat))}</text>')
    # legend
    for bi, name in enumerate(names):
        ly = mt + bi * 22
        col = _PALETTE[bi % len(_PALETTE)]
        s.append(f'<rect x="{ml+pw+24}" y="{ly}" width="14" height="14" '
                 f'fill="{col}" rx="2"/>')
        s.append(f'<text x="{ml+pw+44}" y="{ly+12}">{html.escape(name)}</text>')
    s.append('</svg>')
    with open(path, "w") as fh:
        fh.write("\n".join(s))
    return path


def timeline(path, title, x, series, xlabel="", ylabel=""):
    """series = {name: [y aligned with x]}.  Auto-scales to combined max."""
    W, H = 980, 420
    ml, mr, mt, mb = 70, 170, 56, 60
    pw, ph = W - ml - mr, H - mt - mb
    xs = [float(v) for v in x] or [0.0]
    xmin, xmax = min(xs), max(xs) or 1.0
    xspan = (xmax - xmin) or 1.0
    allvals = [v for vs in series.values() for v in vs if v is not None]
    ymax = (max(allvals) * 1.1) if allvals else 1.0
    ymax = ymax or 1.0

    def X(v): return ml + (float(v) - xmin) / xspan * pw
    def Y(v): return mt + ph - (float(v) / ymax) * ph

    s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'font-family="Segoe UI,Helvetica,Arial,sans-serif" font-size="13">']
    s.append(f'<rect width="{W}" height="{H}" fill="white"/>')
    s.append(f'<text x="{W/2}" y="28" text-anchor="middle" font-size="17" '
             f'font-weight="600">{html.escape(title)}</text>')
    for i in range(6):
        gy = mt + ph - i / 5 * ph
        s.append(f'<line x1="{ml}" y1="{gy:.1f}" x2="{ml+pw}" y2="{gy:.1f}" stroke="#eee"/>')
        s.append(f'<text x="{ml-8}" y="{gy+4:.1f}" text-anchor="end" fill="#6b7280">'
                 f'{_fmt(i/5*ymax)}</text>')
    for i in range(7):
        gx = ml + i / 6 * pw
        xv = xmin + i / 6 * xspan
        s.append(f'<text x="{gx:.1f}" y="{mt+ph+20}" text-anchor="middle" '
                 f'fill="#6b7280">{_fmt(xv)}</text>')
    if xlabel:
        s.append(f'<text x="{ml+pw/2}" y="{H-12}" text-anchor="middle" '
                 f'fill="#374151">{html.escape(xlabel)}</text>')
    if ylabel:
        s.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" fill="#374151" '
                 f'transform="rotate(-90 16 {mt+ph/2})">{html.escape(ylabel)}</text>')
    for ni, (name, ys) in enumerate(series.items()):
        col = _PALETTE[ni % len(_PALETTE)]
        pts = " ".join(f"{X(xx):.1f},{Y(yy):.1f}"
                       for xx, yy in zip(x, ys) if yy is not None)
        if pts:
            s.append(f'<polyline points="{pts}" fill="none" stroke="{col}" '
                     f'stroke-width="2"/>')
        ly = mt + ni * 22
        s.append(f'<rect x="{ml+pw+24}" y="{ly}" width="14" height="14" fill="{col}" rx="2"/>')
        s.append(f'<text x="{ml+pw+44}" y="{ly+12}">{html.escape(name)}</text>')
    s.append('</svg>')
    with open(path, "w") as fh:
        fh.write("\n".join(s))
    return path
