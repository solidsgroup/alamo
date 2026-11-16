#!/usr/bin/env python3
"""Export 1D ALAMO AMReX plotfiles to CSV when VisIt lacks a 1D reader."""

import argparse
import csv
from pathlib import Path

import numpy as np
import yt


BASE_FIELDS = (
    "rigid_eta",
    "temperature",
    "pressure",
    "velocity",
    "component_density_AP_solid",
    "component_density_AP_gas",
    "component_density_Mono",
    "component_density_Premixed",
    "component_density_Primary",
    "component_density_Final",
)


def snapshots(plot_dir):
    visit_file = plot_dir / "celloutput.visit"
    if not visit_file.is_file():
        if (plot_dir / "Header").is_file():
            return [plot_dir]
        raise FileNotFoundError(f"No celloutput.visit or plotfile Header in {plot_dir}")

    result = []
    for line in visit_file.read_text().splitlines():
        entry = line.strip()
        if not entry:
            continue
        path = Path(entry)
        if not path.is_absolute():
            path = plot_dir / path
        if path.name == "Header":
            path = path.parent
        if not (path / "Header").is_file():
            raise FileNotFoundError(f"Invalid plotfile entry in {visit_file}: {entry}")
        result.append(path)
    if not result:
        raise ValueError(f"No plotfile entries in {visit_file}")
    return result


def leaf_profile(dataset, fields):
    grid_info = []
    for grid in dataset.index.grids:
        grid_info.append(
            (
                int(grid.Level),
                float(grid.LeftEdge[0].d),
                float(grid.RightEdge[0].d),
                grid,
            )
        )

    pieces = []
    for level, left, right, grid in grid_info:
        ncell = int(grid.ActiveDimensions[0])
        dx = (right - left) / ncell
        x = left + (np.arange(ncell) + 0.5) * dx
        covered = np.zeros(ncell, dtype=bool)
        for fine_level, fine_left, fine_right, _ in grid_info:
            if fine_level > level:
                covered |= (x >= fine_left) & (x < fine_right)
        keep = ~covered
        piece = {
            "x_m": x[keep],
            "dx_m": np.full(keep.sum(), dx),
            "level": np.full(keep.sum(), level, dtype=int),
        }
        for name in fields:
            values = np.asarray(grid[("boxlib", name)].d).reshape(-1)
            if values.size != ncell:
                raise ValueError(
                    f"Expected {ncell} 1D values for {name}, got {values.size}"
                )
            piece[name] = values[keep]
        pieces.append(piece)

    profile = {
        key: np.concatenate([piece[key] for piece in pieces])
        for key in pieces[0]
    }
    order = np.argsort(profile["x_m"])
    return {key: values[order] for key, values in profile.items()}


def front_metrics(profile):
    eta = profile.get("rigid_eta")
    if eta is None:
        return "", "", ""
    x = profile["x_m"]
    crossings = np.flatnonzero((eta[:-1] >= 0.5) & (eta[1:] < 0.5))
    if not len(crossings):
        return "", "", ""
    i = int(crossings[0])
    fraction = (0.5 - eta[i]) / (eta[i + 1] - eta[i])
    front_x = x[i] + fraction * (x[i + 1] - x[i])
    front_values = {}
    for name in ("temperature", "pressure"):
        values = profile.get(name)
        front_values[name] = (
            values[i] + fraction * (values[i + 1] - values[i])
            if values is not None
            else ""
        )
    return front_x, front_values["temperature"], front_values["pressure"]


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plot_dir", type=Path, help="1D AMReX output directory")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="destination for CSV files (default: <plot_dir>/csv_export)",
    )
    args = parser.parse_args()

    plot_dir = args.plot_dir.resolve()
    output_dir = args.output_dir.resolve() if args.output_dir else plot_dir / "csv_export"
    output_dir.mkdir(parents=True, exist_ok=True)

    time_rows = []
    final_profile = None
    available_fields = None
    for snapshot in snapshots(plot_dir):
        dataset = yt.load(str(snapshot))
        if int(dataset.domain_dimensions[1]) != 1 or int(dataset.domain_dimensions[2]) != 1:
            raise ValueError(f"{snapshot} is not a one-dimensional plotfile")
        available_fields = {
            name for namespace, name in dataset.field_list if namespace == "boxlib"
        }
        if "temperature" not in available_fields or "pressure" not in available_fields:
            raise ValueError(f"{snapshot} lacks temperature or pressure fields")
        fields = [name for name in BASE_FIELDS if name in available_fields]
        fields.extend(
            sorted(
                name
                for name in available_fields
                if name.startswith("component_density_") and name not in fields
            )
        )
        profile = leaf_profile(dataset, fields)
        for name, values in profile.items():
            if name not in ("x_m", "dx_m", "level") and not np.all(np.isfinite(values)):
                raise ValueError(f"{snapshot} contains nonfinite values in {name}")

        front_x, front_temperature, front_pressure = front_metrics(profile)
        row = {
            "time_s": float(dataset.current_time),
            "front_x_m": front_x,
            "front_temperature_K": front_temperature,
            "front_pressure_Pa": front_pressure,
            "temperature_max_K": float(np.max(profile["temperature"])),
            "pressure_min_Pa": float(np.min(profile["pressure"])),
            "pressure_max_Pa": float(np.max(profile["pressure"])),
            "velocity_min_m_s": (
                float(np.min(profile["velocity"])) if "velocity" in profile else ""
            ),
            "velocity_max_m_s": (
                float(np.max(profile["velocity"])) if "velocity" in profile else ""
            ),
        }
        ap_gas = profile.get("component_density_AP_gas")
        row["AP_gas_max_kg_m3"] = float(np.max(ap_gas)) if ap_gas is not None else ""
        component_fields = [
            name for name in fields if name.startswith("component_density_")
        ]
        if component_fields:
            row["mixture_mass_per_area_kg_m2"] = float(
                sum(np.sum(profile[name] * profile["dx_m"]) for name in component_fields)
            )
        else:
            row["mixture_mass_per_area_kg_m2"] = ""
        time_rows.append(row)
        final_profile = profile

    summary_fields = list(time_rows[0])
    write_csv(output_dir / "front_history.csv", time_rows, summary_fields)
    profile_fields = ["x_m", "dx_m", "level", *fields]
    profile_rows = [
        {name: final_profile[name][i] for name in profile_fields}
        for i in range(len(final_profile["x_m"]))
    ]
    write_csv(output_dir / "final_profile.csv", profile_rows, profile_fields)

    final = time_rows[-1]
    print(f"Snapshots exported: {len(time_rows)}")
    print(f"Final time: {final['time_s']:.9g} s")
    print(f"Final eta=0.5 front: {final['front_x_m']} m")
    print(f"Final maximum temperature: {final['temperature_max_K']:.9g} K")
    print(f"CSV output: {output_dir}")


if __name__ == "__main__":
    main()
