#!/usr/bin/env python3
"""Render LowMach aluminum/alumina phase fields as an animated GIF."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from PIL import Image
import yt


FIELD_TYPE = "boxlib"
AL_FIELD = (FIELD_TYPE, "eta_liquid_Al_liquid")
OXIDE_FIELD = (FIELD_TYPE, "eta_liquid_Al2O3_liquid")
PRESSURE_FIELD = (FIELD_TYPE, "pressure")
VELOCITY_X_FIELD = (FIELD_TYPE, "velocityx")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render molten aluminum and liquid alumina plotfiles."
    )
    parser.add_argument("plotfile_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--frames-dir", type=Path)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--fps", type=float, default=5.0)
    parser.add_argument("--dpi", type=int, default=140)
    parser.add_argument(
        "--background",
        choices=("phase", "x_velocity", "pressure"),
        default="phase",
        help="field drawn behind the liquid phases (default: phase)",
    )
    return parser.parse_args()


def plotfiles(root: Path, stride: int) -> list[Path]:
    visit = root / "celloutput.visit"
    paths: list[Path] = []
    if visit.exists():
        for raw in visit.read_text(encoding="utf-8").splitlines():
            entry = raw.strip()
            if entry.endswith("/Header"):
                entry = entry[: -len("/Header")]
            candidate = root / entry
            if entry and candidate.exists():
                paths.append(candidate)
    if not paths:
        paths = sorted(root.glob("*cell"))
    if not paths:
        raise FileNotFoundError(f"No *cell plotfiles found under {root}")
    return paths[::stride]


def load_frame(path: Path) -> dict[str, object]:
    ds = yt.load(str(path))
    grid = ds.covering_grid(0, ds.domain_left_edge, ds.domain_dimensions)
    aluminum = np.asarray(grid[AL_FIELD])[:, :, 0]
    oxide = np.asarray(grid[OXIDE_FIELD])[:, :, 0]
    pressure = np.asarray(grid[PRESSURE_FIELD])[:, :, 0]
    velocity_x = np.asarray(grid[VELOCITY_X_FIELD])[:, :, 0]
    aluminum = np.clip(aluminum, 0.0, 1.0)
    oxide = np.clip(oxide, 0.0, 1.0)
    nx, ny = aluminum.shape
    lo = np.asarray(ds.domain_left_edge.d[:2], dtype=float)
    hi = np.asarray(ds.domain_right_edge.d[:2], dtype=float)
    x = lo[0] + (np.arange(nx) + 0.5) * (hi[0] - lo[0]) / nx
    y = lo[1] + (np.arange(ny) + 0.5) * (hi[1] - lo[1]) / ny
    xx, yy = np.meshgrid(x, y, indexing="ij")
    oxide_mass = float(oxide.sum())
    aluminum_mass = float(aluminum.sum())
    oxide_centroid = np.array(
        [(oxide * xx).sum(), (oxide * yy).sum()]
    ) / oxide_mass
    aluminum_centroid = np.array(
        [(aluminum * xx).sum(), (aluminum * yy).sum()]
    ) / aluminum_mass
    return {
        "path": path,
        "time": float(ds.current_time),
        "extent": (lo[0], hi[0], lo[1], hi[1]),
        "x": x,
        "y": y,
        "aluminum": aluminum,
        "oxide": oxide,
        "pressure": pressure,
        "velocity_x": velocity_x,
        "oxide_mass": oxide_mass,
        "oxide_centroid": oxide_centroid,
        "aluminum_centroid": aluminum_centroid,
    }


def render(
    frame: dict[str, object],
    relative_history: list[np.ndarray],
    initial_relative_centroid: np.ndarray,
    output: Path,
    dpi: int,
    pressure_limit: float,
    velocity_limit: float,
    background: str,
) -> None:
    aluminum = frame["aluminum"]
    oxide = frame["oxide"]
    total = np.clip(aluminum + oxide, 0.0, 1.0)
    al_color = np.array([0.45, 0.72, 0.82])
    oxide_color = np.array([1.0, 0.28, 0.025])
    pressure = np.asarray(frame["pressure"])
    pressure_delta = pressure - np.median(pressure)
    velocity_x = np.asarray(frame["velocity_x"])

    extent = np.asarray(frame["extent"], dtype=float) * 1.0e6
    x = np.asarray(frame["x"], dtype=float) * 1.0e6
    y = np.asarray(frame["y"], dtype=float) * 1.0e6
    fig, ax = plt.subplots(figsize=(8.4, 5.45), constrained_layout=True)
    flow_image = None
    if background == "pressure":
        flow_image = ax.imshow(
            (pressure_delta * 1.0e-3).T,
            origin="lower",
            extent=extent,
            interpolation="bilinear",
            cmap="RdBu_r",
            vmin=-pressure_limit * 1.0e-3,
            vmax=pressure_limit * 1.0e-3,
        )
        flow_label = "pressure perturbation [kPa]"
        title_field = "pressure perturbation"
    elif background == "x_velocity":
        velocity_cmap = matplotlib.colormaps["RdBu_r"].copy()
        flow_image = ax.imshow(
            velocity_x.T,
            origin="lower",
            extent=extent,
            interpolation="bilinear",
            cmap=velocity_cmap,
            vmin=-velocity_limit,
            vmax=velocity_limit,
        )
        flow_label = "x-velocity, $u_x$ [m/s]"
        title_field = "x-velocity"
    else:
        ax.set_facecolor("#17202a")
        title_field = "phase fields"
    if flow_image is not None:
        fig.colorbar(
            flow_image, ax=ax, fraction=0.045, pad=0.025, label=flow_label
        )

    al_overlay = np.empty((aluminum.shape[1], aluminum.shape[0], 4))
    al_overlay[:, :, :3] = al_color
    al_overlay[:, :, 3] = 0.88 * aluminum.T
    oxide_overlay = np.empty((oxide.shape[1], oxide.shape[0], 4))
    oxide_overlay[:, :, :3] = oxide_color
    oxide_overlay[:, :, 3] = 0.94 * oxide.T
    ax.imshow(al_overlay, origin="lower", extent=extent, interpolation="bicubic")
    ax.imshow(
        oxide_overlay, origin="lower", extent=extent, interpolation="bicubic"
    )
    ax.contour(
        x,
        y,
        total.T,
        levels=[0.5],
        colors=["white"],
        linewidths=1.2,
    )
    ax.contour(
        x,
        y,
        oxide.T,
        levels=[0.5],
        colors=["#ffd166"],
        linewidths=1.1,
    )

    aluminum_centroid = np.asarray(frame["aluminum_centroid"])
    oxide_centroid = np.asarray(frame["oxide_centroid"])
    relative_centroid = oxide_centroid - aluminum_centroid
    relative_displacement = relative_centroid - initial_relative_centroid
    # Anchor the entire oxide history at the current aluminum centroid.  This
    # removes aluminum translation while retaining changes in the oxide/Al
    # separation vector.
    trail = (
        aluminum_centroid + np.asarray(relative_history)
    ) * 1.0e6
    reference_centroid = (
        aluminum_centroid + initial_relative_centroid
    ) * 1.0e6
    aluminum_centroid_um = aluminum_centroid * 1.0e6
    oxide_centroid_um = oxide_centroid * 1.0e6
    ax.plot(
        trail[:, 0], trail[:, 1], color="#ffd166", linewidth=1.3,
        marker="o", markersize=2.8, alpha=0.9,
    )
    ax.scatter(
        aluminum_centroid_um[0], aluminum_centroid_um[1],
        marker="+", s=65, linewidths=1.7, color="#7fdbff", zorder=6,
    )
    ax.scatter(
        reference_centroid[0], reference_centroid[1], marker="o", s=55,
        linewidths=1.4, facecolors="none", edgecolors="white", zorder=6,
    )
    ax.scatter(
        oxide_centroid_um[0], oxide_centroid_um[1], marker="x", s=50,
        linewidths=1.7, color="#ffd166", zorder=7,
    )
    ax.annotate(
        "", xy=oxide_centroid_um, xytext=reference_centroid,
        arrowprops={"arrowstyle": "->", "color": "#ffd166", "lw": 1.8},
        zorder=6,
    )

    displacement_um = relative_displacement * 1.0e6
    displacement_norm_um = np.linalg.norm(relative_displacement) * 1.0e6
    ax.set_title(
        f"Oxide motion relative to aluminum over {title_field}\n"
        f"t = {float(frame['time']) * 1.0e6:6.2f} µs, "
        f"Δr$_{{oxide/Al}}$ = ({displacement_um[0]:5.2f}, "
        f"{displacement_um[1]:5.2f}) µm, "
        f"|Δr| = {displacement_norm_um:5.2f} µm"
    )
    ax.set_xlabel("x [µm]")
    ax.set_ylabel("y [µm]")
    ax.set_aspect("equal")
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.legend(
        handles=[
            Patch(facecolor=al_color, edgecolor="white", label="Al liquid"),
            Patch(facecolor=oxide_color, edgecolor="#ffd166", label="Al₂O₃ liquid"),
            Line2D([], [], color="#ffd166", marker="o", markersize=3,
                   label="oxide path in Al frame"),
            Line2D([], [], color="white", marker="o", markerfacecolor="none",
                   linestyle="none", label="initial relative position"),
        ],
        loc="upper left",
        framealpha=0.88,
    )
    fig.savefig(output, dpi=dpi, facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if args.stride < 1 or args.fps <= 0.0:
        raise ValueError("--stride must be positive and --fps must exceed zero")
    paths = plotfiles(args.plotfile_root, args.stride)
    frames_dir = args.frames_dir or args.output.with_suffix("")
    frames_dir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    yt.set_log_level(50)

    data = [load_frame(path) for path in paths]
    pressure_samples = np.concatenate([
        np.abs(
            np.asarray(frame["pressure"]) -
            np.median(np.asarray(frame["pressure"]))
        ).ravel()
        for frame in data
    ])
    pressure_limit = max(float(np.percentile(pressure_samples, 99.0)), 1.0)
    gas_velocity_x_samples = np.concatenate([
        np.abs(np.asarray(frame["velocity_x"]))[
            np.asarray(frame["aluminum"]) + np.asarray(frame["oxide"]) < 0.20
        ]
        for frame in data
    ])
    velocity_limit = max(
        float(np.percentile(gas_velocity_x_samples, 99.0)), 1.0e-6
    )
    initial_centroid = np.asarray(data[0]["oxide_centroid"])
    initial_relative_centroid = (
        initial_centroid - np.asarray(data[0]["aluminum_centroid"])
    )
    relative_history: list[np.ndarray] = []
    pngs: list[Path] = []
    for number, frame in enumerate(data):
        relative_history.append(
            np.asarray(frame["oxide_centroid"]) -
            np.asarray(frame["aluminum_centroid"])
        )
        png = frames_dir / f"frame_{number:04d}.png"
        render(
            frame, relative_history, initial_relative_centroid, png, args.dpi,
            pressure_limit, velocity_limit, args.background
        )
        pngs.append(png)

    images = [Image.open(path).convert("RGB") for path in pngs]
    duration_ms = round(1000.0 / args.fps)
    images[0].save(
        args.output,
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0,
        optimize=True,
    )
    for image in images:
        image.close()

    final_displacement = np.linalg.norm(
        np.asarray(data[-1]["oxide_centroid"]) - initial_centroid
    )
    relative_displacement = np.linalg.norm(
        (np.asarray(data[-1]["oxide_centroid"]) -
         np.asarray(data[-1]["aluminum_centroid"])) -
        (initial_centroid - np.asarray(data[0]["aluminum_centroid"]))
    )
    mass_drift = (
        float(data[-1]["oxide_mass"]) / float(data[0]["oxide_mass"]) - 1.0
    )
    print(f"Rendered {len(pngs)} frames to {args.output}")
    print(f"Final oxide displacement: {final_displacement * 1.0e6:.6g} um")
    print(f"Motion relative to aluminum: {relative_displacement * 1.0e6:.6g} um")
    print(f"Relative oxide phase-volume drift: {mass_drift:.6g}")


if __name__ == "__main__":
    main()
