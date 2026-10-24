#!/usr/bin/env python3
"""Extract published Fig. 4 vectors, without tracing raster pixels or fitting data.

Requires the system `mutool` executable and numpy/matplotlib. Coordinate bounds
and path selections are specific to the user-supplied Chen et al. (2002) PDF.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

import numpy as np

HERE = Path(__file__).resolve().parent
X0, X1, Y0, Y1 = 73.394, 238.013, 590.031, 719.816


def subpaths(path):
    result, current = [], []
    for item in path:
        if item.tag == "moveto" and current:
            result.append(current)
            current = []
        current.append(item)
    if current:
        result.append(current)
    return result


def xy(item):
    return np.array([float(item.attrib[k]) for k in ("x", "y")])


def transform(points):
    return (np.asarray(points) - [X0, Y0]) * [1 / (X1-X0), 3.5 / (Y1-Y0)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    trace = subprocess.check_output(
        ["mutool", "draw", "-F", "trace", str(args.pdf), "5"],
        stderr=subprocess.DEVNULL)
    paths = ET.fromstring(trace).findall(".//stroke_path")
    assert len(paths) == 49, "PDF differs from inspected document; re-audit paths"
    assert np.allclose(xy(paths[30][0]), [X0, Y0])
    curves, rows = {}, []
    for path_id, kind, qs in [(47, "eq15_solid", (500,1000,200)),
                             (48, "single_temperature_dashed", (200,500,1000))]:
        assert len(subpaths(paths[path_id])) == 3
        for q, segment in zip(qs, subpaths(paths[path_id])):
            points = transform([xy(p) for p in segment])
            if kind == "eq15_solid":
                curves[q] = points
            for t, r in points:
                rows.append(dict(kind=kind, q_cal_cm2_s=q, t=t, r_cm_s=r,
                                 pdf_stroke_path=path_id))
    for segment in subpaths(paths[46]):
        if len(segment) == 5 and segment[1].tag == "curveto":
            # Circle's four cubic Bezier endpoints give its exact center.
            endpoints = [xy(segment[0])] + [
                [float(p.attrib["x3"]),float(p.attrib["y3"])] for p in segment[1:]]
            center = (np.min(endpoints,axis=0)+np.max(endpoints,axis=0))/2
            kind = "dns_2d_circle"
        elif len(segment) == 5 and segment[-1].tag == "closepath":
            center = np.mean([xy(p) for p in segment[:-1]],axis=0)
            kind = "eq19_square"
        elif len(segment) == 2 and abs(xy(segment[0])[1]-xy(segment[1])[1]) < .002:
            # Each asterisk consists of horizontal/vertical/diagonal lines.
            # Read only its horizontal stroke, once per asterisk.
            center = (xy(segment[0])+xy(segment[1]))/2
            kind = "dns_3d_asterisk"
        else:
            continue
        t, r = transform(center)
        q = min(curves, key=lambda q: abs(r-np.interp(t,curves[q][:,0],curves[q][:,1])))
        rows.append(dict(kind=kind,q_cal_cm2_s=q,t=t,r_cm_s=r,pdf_stroke_path=46))
    for filename, subset in [("figure4_curves.csv",[r for r in rows if r["pdf_stroke_path"]!=46]),
                             ("figure4_markers.csv",[r for r in rows if r["pdf_stroke_path"]==46])]:
        with (args.out/filename).open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(subset)
    provenance = dict(
        source=str(args.pdf),sha256=hashlib.sha256(args.pdf.read_bytes()).hexdigest(),
        reference="Chen, Buckmaster, Jackson and Massa (2002), Proc. Combust. Inst. 29, 2923–2929, Fig. 4, printed p.2927 / PDF page 5",
        extraction="mutool draw -F trace; exact vector path vertices / marker centers; no fit or pixel trace",
        pdf_bounds=dict(x0=X0,x1=X1,y0=Y0,y1=Y1,t_min=0,t_max=1,r_min=0,r_max=3.5),
        path_ids=dict(markers=46,solid=47,dashed=48),
        coordinate_precision_pdf_points=.001,line_width_pdf_points=.184,
        conservative_graphical_resolution=dict(t=.184/(X1-X0),r_cm_s=.184*3.5/(Y1-Y0)),
        uncertainty_note="Stored decimals preserve vector extraction, not simulation/experimental precision. One full printed stroke width is a conservative graphical resolution, not a statistical confidence interval. Source DNS sampling/numerical uncertainties are unavailable. Endpoint DNS markers overlap and should not be treated as independent validation observations.",
        marker_counts={k:sum(r["kind"]==k for r in rows) for k in ("eq19_square","dns_2d_circle","dns_3d_asterisk")})
    (args.out/"figure4_provenance.json").write_text(json.dumps(provenance,indent=2)+"\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":9,"pdf.fonttype":42})
    fig,ax=plt.subplots(figsize=(5.8,3.8),layout="constrained")
    colors={200:"#0072B2",500:"#D55E00",1000:"#009E73"}
    for q in (200,500,1000):
        pts=curves[q]
        ax.plot(*pts.T,color=colors[q],label=f"q = {q} cal cm⁻² s⁻¹")
        for kind,marker in [("dns_2d_circle","o"),("dns_3d_asterisk","*"),("eq19_square","s")]:
            subset=[r for r in rows if r["kind"]==kind and r["q_cal_cm2_s"]==q]
            ax.plot([r["t"] for r in subset],[r["r_cm_s"] for r in subset],linestyle="none",marker=marker,ms=4,mfc="none",mec=colors[q],mew=.7)
    ax.set(xlim=(0,1),ylim=(0,3.5),xlabel="AP volume fraction, t",ylabel="Regression speed (cm s⁻¹)",title="Chen et al. (2002), Figure 4 — extracted vectors")
    ax.legend(frameon=False)
    ax.text(.02,.02,"Lines: Eq. 15   □ Eq. 19   ○ 2D DNS   ✶ 3D DNS",transform=ax.transAxes,fontsize=8)
    for suffix in ("png","pdf"):
        fig.savefig(args.out/f"figure4_extracted.{suffix}",dpi=240)
    print(json.dumps(provenance["marker_counts"]))


if __name__ == "__main__":
    main()
