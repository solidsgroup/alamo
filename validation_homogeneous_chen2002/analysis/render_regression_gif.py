#!/usr/bin/env python3
"""Animate saved LowMach solid/temperature fields without rerunning the solver."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/chen-regression-gif-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image, ImageDraw
import yt

from extract_runs import crossing

yt.set_log_level(40)
STUDY = Path(__file__).resolve().parent.parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--frame-ms", type=int, default=150)
    args = parser.parse_args()
    assert args.frame_ms >= 20
    folder = args.case.resolve()
    case = json.loads((folder / "case.json").read_text())
    receipt = json.loads((folder / "run.json").read_text())
    assert receipt["returncode"] == 0
    assert receipt["input_sha256"] == case["input_sha256"] == digest(folder / "input")
    selected_path = STUDY / "parameters_current.json"
    assert case["parameters"] == json.loads(selected_path.read_text()), "Choose the selected parameter set"
    assert "finalized" in (folder / "output/out.log").read_text()
    plotfiles = sorted(path.parent for path in (folder / "output").glob("*cell/Header"))
    assert len(plotfiles) >= 3
    frames = []
    for path in plotfiles:
        ds = yt.load(str(path))
        lev = ds.index.max_level
        dims = ds.domain_dimensions * ds.refine_by**lev
        if ds.dimensionality < 3:
            dims[ds.dimensionality:] = 1
        grid = ds.covering_grid(lev, ds.domain_left_edge, dims)
        eta = np.asarray(grid["boxlib", "rigid_eta"]).squeeze(axis=2)
        temperature = np.asarray(grid["boxlib", "temperature"]).squeeze(axis=2)
        left, right = np.asarray(ds.domain_left_edge), np.asarray(ds.domain_right_edge)
        x = left[0]+(np.arange(dims[0])+.5)*(right[0]-left[0])/dims[0]
        y = left[1]+(np.arange(dims[1])+.5)*(right[1]-left[1])/dims[1]
        values = np.array([crossing(y, e, t, .5) for e, t in zip(eta, temperature)])
        frames.append(dict(time_s=float(ds.current_time), eta=eta.copy(), T=temperature.copy(),
            surface_y_m=values[:, 0], surface_T_K=float(values[:, 1].mean()),
            x_m=x, y_m=y, extent_um=1.e6*np.array([left[0], right[0], left[1], right[1]]),
            plotfile=str(path), header_sha256=digest(path / "Header")))
        print("Read", path.name, flush=True)
    assert np.all(np.diff([f["time_s"] for f in frames]) > 0)
    assert frames[-1]["time_s"] >= .999*case["duration_s"]
    assert all(np.array_equal(f["extent_um"], frames[0]["extent_um"]) for f in frames)
    initial_y = float(frames[0]["surface_y_m"].mean())
    positions = np.array([f["surface_y_m"].mean() for f in frames])
    assert np.all(np.diff(positions) < 0), "Surface should recede monotonically"

    q = case["reference"]["heat_flux_cal_cm2_s"]
    fraction = case["reference"]["ap_volume_fraction"]
    tmin = case["T0_K"]
    tmax = 50*np.ceil(max(float(f["T"].max()) for f in frames)/50)
    norm = Normalize(tmin, tmax)
    cmap = plt.get_cmap("inferno")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig = plt.figure(figsize=(8.3, 7.6), dpi=135, facecolor="white")
    ax = fig.add_axes([.13, .18, .65, .62], facecolor="#eaf1f6")
    cax = fig.add_axes([.82, .18, .024, .62])
    bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    bar.set_label("Propellant temperature (K)", labelpad=12)
    bar.set_ticks(np.arange(tmin, tmax+1, 100))
    f = frames[0]
    rgba = cmap(norm(f["T"].T))
    rgba[..., 3] = np.clip(f["eta"].T, 0, 1)
    im = ax.imshow(rgba, extent=f["extent_um"], origin="lower", aspect="auto",
                   interpolation="nearest", zorder=1)
    xmin, xmax, ymin, ymax = f["extent_um"]
    width, height = xmax-xmin, ymax-ymin
    initial_line = ax.axhline(initial_y*1.e6, color="#397da6", lw=1.5, ls="--", zorder=4)
    front, = ax.plot(f["x_m"]*1.e6, f["surface_y_m"]*1.e6,
                     color="#22c9c3", lw=2.6, zorder=5)
    ax.text(.5, .07, "PROPELLANT", transform=ax.transAxes, ha="center", va="center",
            color="white", fontsize=13, weight="bold", zorder=6)
    ax.text(xmin+.5*width, initial_y*1.e6+.05*height, "GAS", color="#41566b",
            ha="center", fontsize=12, weight="bold", zorder=6)
    heat_y0, heat_y1 = np.array(case["source_y_m"])*1.e6
    ax.axhspan(heat_y0, heat_y1, facecolor="#edab58", alpha=.23, zorder=2)
    ax.text(xmin+.5*width, (heat_y0+heat_y1)/2, "Prescribed heat input",
            ha="center", va="center", color="#875318", fontsize=10, zorder=6)
    for xfrac in (.22, .5, .78):
        xp = xmin+xfrac*width
        ax.annotate("", (xp, heat_y0-.072*height), (xp, heat_y0-.008*height),
                    arrowprops=dict(arrowstyle="-|>", lw=1.7, color="#a76a24"), zorder=6)
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), xlabel="Transverse position (µm)",
           ylabel="Vertical position (µm)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(.13, .948, "Propellant regression", fontsize=20, weight="bold", color="#192f44")
    fig.text(.13, .91, f"Homogeneous mixture  •  {100*fraction:g}% AP by volume  •  q = {q:g} cal/(cm²·s)",
             fontsize=11, color="#4d6277")
    time_label = fig.text(.13, .856, "", fontsize=12, color="#192f44", family="DejaVu Sans Mono")
    progress_background = fig.add_axes([.13, .823, .65, .007])
    progress_background.set(xlim=(0, 1), ylim=(0, 1))
    progress_background.axis("off")
    progress_background.axhspan(0, 1, color="#e2e9f0")
    progress_line, = progress_background.plot([0, 0], [.5, .5], color="#25828b", lw=5)
    fig.legend(handles=[initial_line, Line2D([], [], color="#22c9c3", lw=2.6)],
               labels=["Initial surface", "Moving surface (solid fraction = 0.5)"],
               loc="lower center", bbox_to_anchor=(.5, .075), ncol=2, frameon=False, fontsize=10)
    fig.text(.5, .038, "Transverse scale stretched for visibility  •  Saved simulation fields  •  Frozen gas chemistry",
             ha="center", fontsize=9, color="#637386")

    images = []
    frame_data = []
    for f in frames:
        rgba = cmap(norm(f["T"].T))
        rgba[..., 3] = np.clip(f["eta"].T, 0, 1)
        im.set_data(rgba)
        front.set_data(f["x_m"]*1.e6, f["surface_y_m"]*1.e6)
        depth = 1.e6*(initial_y-f["surface_y_m"].mean())
        time_label.set_text(f"Time {f['time_s']*1000:5.3f} ms     Recession {depth:5.1f} µm")
        progress_line.set_data([0, f["time_s"]/frames[-1]["time_s"]], [.5, .5])
        fig.canvas.draw()
        images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba()).copy()).convert("RGB"))
        frame_data.append(dict(time_ms=1000*f["time_s"], recession_um=float(depth),
            surface_y_um=float(f["surface_y_m"].mean()*1.e6), surface_temperature_K=f["surface_T_K"],
            plotfile=f["plotfile"], header_sha256=f["header_sha256"]))
    plt.close(fig)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    durations = np.rint(np.diff([f["time_s"] for f in frames])/np.median(
        np.diff([f["time_s"] for f in frames]))*args.frame_ms/10).astype(int)*10
    durations = [int(t) for t in durations]+[1200]
    durations[0] += 500
    images[0].save(args.output, save_all=True, append_images=images[1:], duration=durations,
                   loop=0, optimize=True, disposal=2)
    # Still previews allow inspection of the GIF's beginning, middle and end.
    preview_path = args.output.with_name(args.output.stem+"_preview.png")
    thumbs = []
    for idx in (0, len(images)//2, len(images)-1):
        thumb = images[idx].copy()
        thumb.thumbnail((560, 520), Image.Resampling.LANCZOS)
        thumbs.append(thumb)
    contact = Image.new("RGB", (sum(t.width for t in thumbs), max(t.height for t in thumbs)), "white")
    offset = 0
    for thumb in thumbs:
        contact.paste(thumb, (offset, 0))
        offset += thumb.width
    contact.save(preview_path)
    images[-1].save(args.output.with_name(args.output.stem+"_final.png"))
    with args.output.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(frame_data[0]))
        writer.writeheader()
        writer.writerows(frame_data)
    with Image.open(args.output) as gif:
        assert gif.n_frames == len(images)
        assert gif.size == images[0].size
        assert gif.info["loop"] == 0
    manifest = dict(case=folder.name, parameters_sha256=digest(selected_path),
        input_sha256=digest(folder / "input"), recorded_binary_sha256=receipt["binary_sha256"],
        script_sha256=digest(Path(__file__)), frame_count=len(images),
        duration_ms=sum(durations), simulation_duration_ms=1000*frames[-1]["time_s"],
        final_recession_um=frame_data[-1]["recession_um"], temperature_scale_K=[tmin, tmax],
        frame_durations_ms=durations, transverse_scale_stretched=True,
        interface_definition="rigid_eta=0.5, linearly interpolated in each transverse column",
        temperature_display="Raw temperature colored with fixed scale; opacity equals solid fraction",
        gif_sha256=digest(args.output), gif_bytes=args.output.stat().st_size)
    args.output.with_suffix(".json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps(manifest, indent=2))
    print("GIF:", args.output)
    print("Preview:", preview_path)


if __name__ == "__main__":
    main()
