#!/usr/bin/env python3
"""Extract local regression data along the rigid_eta=0.5 burn surface."""

import argparse
import csv
import math
import re
from pathlib import Path

import contourpy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt

yt.set_log_level(50)


FIELD_NAMES = {
    "temperature": "temperature",
    "pressure": "pressure",
    "eta": "rigid_eta",
    "eta_ap": "rigid_species_eta_AP_solid",
    "eta_htpb": "rigid_species_eta_HTPB_solid",
}

CSV_COLUMNS = [
    "plotfile", "time_s", "contour_id", "point_id",
    "x_m", "y_m", "arc_length_m", "contour_length_m", "contour_closed",
    "normal_x", "normal_y", "grad_eta_per_m", "eta_dot_per_s",
    "AP_eta_dot_per_s", "HTPB_eta_dot_per_s",
    "mass_flux_kg_m2_s", "AP_mass_flux_kg_m2_s",
    "HTPB_mass_flux_kg_m2_s", "temperature_K", "pressure_Pa",
    "thermal_conductivity_W_m_K", "normal_temperature_gradient_K_m",
    "heat_flux_into_solid_W_m2", "heat_flux_outward_W_m2",
    "AP_solid_fraction", "HTPB_solid_fraction", "species",
    "reference_density_kg_m3", "distance_to_interface_m",
    "signed_distance_to_interface_m",
]


def metadata_values(path):
    values = {}
    for line in path.read_text().splitlines():
        if "=" not in line or line.lstrip().startswith("#"):
            continue
        key, value = line.split("=", 1)
        values[key.strip()] = value.strip().strip('"')
    return values


def leading_float(value):
    match = re.match(r"\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)", value)
    if not match:
        raise ValueError(f"Could not read a number from {value!r}")
    return float(match.group(1))


def plotfile_time(path):
    lines = (path / "Header").read_text().splitlines()
    field_count = int(lines[1])
    return float(lines[3 + field_count])


def smoother_step(eta):
    eta_lo = 0.11920292202211755
    x = np.clip((eta - eta_lo) / (1.0 - 2.0 * eta_lo), 0.0, 1.0)
    return x**3 * (x * (6.0 * x - 15.0) + 10.0)


def load_snapshot(path, level):
    dataset = yt.load(str(path))
    refinement = int(dataset.refine_by) ** level
    dimensions = np.asarray(dataset.domain_dimensions, dtype=int) * refinement
    dimensions[2] = 1
    grid = dataset.covering_grid(
        level=level, left_edge=dataset.domain_left_edge, dims=dimensions)

    fields = {}
    for name, plot_name in FIELD_NAMES.items():
        fields[name] = np.asarray(grid[("boxlib", plot_name)]).squeeze()
        if fields[name].ndim != 2:
            raise RuntimeError(f"{path}/{plot_name} is not two dimensional")
        if not np.all(np.isfinite(fields[name])):
            raise RuntimeError(f"{path}/{plot_name} contains non-finite values")

    lower = np.asarray(dataset.domain_left_edge, dtype=float)
    upper = np.asarray(dataset.domain_right_edge, dtype=float)
    nx, ny = fields["eta"].shape
    fields["x"] = lower[0] + (np.arange(nx) + 0.5) * (upper[0] - lower[0]) / nx
    fields["y"] = lower[1] + (np.arange(ny) + 0.5) * (upper[1] - lower[1]) / ny
    fields["time"] = float(dataset.current_time)
    fields["path"] = path
    return fields


def resample_contour(line, spacing):
    line = np.asarray(line, dtype=float)
    closed = len(line) > 2 and np.linalg.norm(line[0] - line[-1]) < 0.25 * spacing
    if closed:
        line = line[:-1]
        work = np.vstack((line, line[0]))
    else:
        work = line

    lengths = np.linalg.norm(np.diff(work, axis=0), axis=1)
    keep = np.concatenate(([True], lengths > 1.0e-14 * spacing))
    work = work[keep]
    if closed:
        if len(work) < 4:
            return None
        if not np.array_equal(work[-1], work[0]):
            work = np.vstack((work, work[0]))
    elif len(work) < 2:
        return None

    cumulative = np.concatenate(([0.0], np.cumsum(np.linalg.norm(
        np.diff(work, axis=0), axis=1))))
    total = cumulative[-1]
    if not total > spacing:
        return None
    arc = np.arange(0.0, total, spacing)
    if not closed:
        arc = np.append(arc, total)
    points = np.column_stack((
        np.interp(arc, cumulative, work[:, 0]),
        np.interp(arc, cumulative, work[:, 1]),
    ))
    return points, arc, total, closed


def interface_distances(arc, total, closed, ap_fraction):
    count = len(arc)
    pairs = [(i, i + 1) for i in range(count - 1)]
    if closed:
        pairs.append((count - 1, 0))

    crossings = []
    shifted = ap_fraction - 0.5
    for left, right in pairs:
        left_value = shifted[left]
        right_value = shifted[right]
        if left_value == 0.0:
            crossings.append(arc[left])
        if left_value * right_value < 0.0:
            segment = ((arc[right] - arc[left]) if right else
                       (total - arc[left]))
            fraction = abs(left_value) / (abs(left_value) + abs(right_value))
            crossings.append((arc[left] + fraction * segment) % total)

    if not crossings:
        return np.full(count, np.nan)
    crossings = np.unique(np.asarray(crossings))
    distance = np.min(np.abs(arc[:, None] - crossings[None, :]), axis=1)
    if closed:
        distance = np.minimum(distance, total - distance)
    return distance


def interpolator(x, y, values):
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    def sample(points):
        x_index = (points[:, 0] - x[0]) / dx
        y_index = (points[:, 1] - y[0]) / dy
        i = np.clip(np.floor(x_index).astype(int), 0, len(x) - 2)
        j = np.clip(np.floor(y_index).astype(int), 0, len(y) - 2)
        x_weight = np.clip(x_index - i, 0.0, 1.0)
        y_weight = np.clip(y_index - j, 0.0, 1.0)
        return (
            (1.0 - x_weight) * (1.0 - y_weight) * values[i, j]
            + x_weight * (1.0 - y_weight) * values[i + 1, j]
            + (1.0 - x_weight) * y_weight * values[i, j + 1]
            + x_weight * y_weight * values[i + 1, j + 1])

    return sample


def extract_snapshot(current, previous, following, properties, output_path):
    x = current["x"]
    y = current["y"]
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    if previous is None:
        dt = following["time"] - current["time"]
        eta_ap_dot = (following["eta_ap"] - current["eta_ap"]) / dt
        eta_htpb_dot = (following["eta_htpb"] - current["eta_htpb"]) / dt
        eta_dot = (following["eta"] - current["eta"]) / dt
    elif following is None:
        dt = current["time"] - previous["time"]
        eta_ap_dot = (current["eta_ap"] - previous["eta_ap"]) / dt
        eta_htpb_dot = (current["eta_htpb"] - previous["eta_htpb"]) / dt
        eta_dot = (current["eta"] - previous["eta"]) / dt
    else:
        dt = following["time"] - previous["time"]
        eta_ap_dot = (following["eta_ap"] - previous["eta_ap"]) / dt
        eta_htpb_dot = (following["eta_htpb"] - previous["eta_htpb"]) / dt
        eta_dot = (following["eta"] - previous["eta"]) / dt
    if not dt > 0.0:
        raise RuntimeError(f"Non-positive plotfile time interval around {current['path']}")

    grad_eta_x, grad_eta_y = np.gradient(current["eta"], dx, dy, edge_order=2)
    grad_t_x, grad_t_y = np.gradient(
        current["temperature"], dx, dy, edge_order=2)

    condensed_eta = np.maximum(current["eta_ap"], 0.0) + np.maximum(
        current["eta_htpb"], 0.0)
    solid_fraction = smoother_step(condensed_eta)
    solid_scale = np.divide(
        solid_fraction, condensed_eta,
        out=np.zeros_like(condensed_eta), where=condensed_eta > 0.0)
    gas_conductivity = properties["gas_lambda_a"] * current["temperature"] + \
        properties["gas_lambda_b"]
    conductivity = (
        (1.0 - solid_fraction) * gas_conductivity
        + solid_scale * (
            np.maximum(current["eta_ap"], 0.0) * properties["k_ap"]
            + np.maximum(current["eta_htpb"], 0.0) * properties["k_htpb"]))

    values = {
        "temperature": current["temperature"],
        "pressure": current["pressure"],
        "eta_ap": current["eta_ap"],
        "eta_htpb": current["eta_htpb"],
        "eta_dot": eta_dot,
        "eta_ap_dot": eta_ap_dot,
        "eta_htpb_dot": eta_htpb_dot,
        "grad_eta_x": grad_eta_x,
        "grad_eta_y": grad_eta_y,
        "grad_t_x": grad_t_x,
        "grad_t_y": grad_t_y,
        "conductivity": conductivity,
    }
    sample = {name: interpolator(x, y, value) for name, value in values.items()}

    contour_generator = contourpy.contour_generator(
        x=x, y=y, z=current["eta"].T, name="serial", line_type="Separate")
    lines = contour_generator.lines(0.5)
    rows = []
    contour_summary = []
    spacing = min(dx, dy)
    for contour_id, raw_line in enumerate(lines):
        resampled = resample_contour(raw_line, spacing)
        if resampled is None:
            continue
        points, arc, contour_length, closed = resampled
        # Contours ending on a cell-center boundary can differ from that
        # boundary by one floating-point ulp after arc-length resampling.
        points[:, 0] = np.clip(points[:, 0], x[0], x[-1])
        points[:, 1] = np.clip(points[:, 1], y[0], y[-1])
        local = {name: function(points) for name, function in sample.items()}
        if any(np.any(~np.isfinite(value)) for value in local.values()):
            raise RuntimeError(
                f"Non-finite interpolated data on contour {contour_id} in {current['path']}")

        grad_eta = np.hypot(local["grad_eta_x"], local["grad_eta_y"])
        if np.any(grad_eta <= 0.0):
            raise RuntimeError(
                f"Zero eta gradient on contour {contour_id} in {current['path']}")
        normal_x = -local["grad_eta_x"] / grad_eta
        normal_y = -local["grad_eta_y"] / grad_eta
        normal_temperature_gradient = (
            local["grad_t_x"] * normal_x + local["grad_t_y"] * normal_y)
        heat_flux_into_solid = local["conductivity"] * normal_temperature_gradient

        eta_sum = np.maximum(local["eta_ap"], 0.0) + np.maximum(
            local["eta_htpb"], 0.0)
        ap_fraction = np.divide(
            np.maximum(local["eta_ap"], 0.0), eta_sum,
            out=np.full_like(eta_sum, 0.5), where=eta_sum > 0.0)
        ap_fraction = np.clip(ap_fraction, 0.0, 1.0)
        htpb_fraction = 1.0 - ap_fraction
        species_ap = ap_fraction >= 0.5
        distance = interface_distances(
            arc, contour_length, closed, ap_fraction)

        ap_mass_flux = -properties["rho_ap"] * local["eta_ap_dot"] / grad_eta
        htpb_mass_flux = -properties["rho_htpb"] * local["eta_htpb_dot"] / grad_eta
        mass_flux = ap_mass_flux + htpb_mass_flux

        for point_id, point in enumerate(points):
            rows.append({
                "plotfile": current["path"].name,
                "time_s": current["time"],
                "contour_id": contour_id,
                "point_id": point_id,
                "x_m": point[0],
                "y_m": point[1],
                "arc_length_m": arc[point_id],
                "contour_length_m": contour_length,
                "contour_closed": int(closed),
                "normal_x": normal_x[point_id],
                "normal_y": normal_y[point_id],
                "grad_eta_per_m": grad_eta[point_id],
                "eta_dot_per_s": local["eta_dot"][point_id],
                "AP_eta_dot_per_s": local["eta_ap_dot"][point_id],
                "HTPB_eta_dot_per_s": local["eta_htpb_dot"][point_id],
                "mass_flux_kg_m2_s": mass_flux[point_id],
                "AP_mass_flux_kg_m2_s": ap_mass_flux[point_id],
                "HTPB_mass_flux_kg_m2_s": htpb_mass_flux[point_id],
                "temperature_K": local["temperature"][point_id],
                "pressure_Pa": local["pressure"][point_id],
                "thermal_conductivity_W_m_K": local["conductivity"][point_id],
                "normal_temperature_gradient_K_m": normal_temperature_gradient[point_id],
                "heat_flux_into_solid_W_m2": heat_flux_into_solid[point_id],
                "heat_flux_outward_W_m2": -heat_flux_into_solid[point_id],
                "AP_solid_fraction": ap_fraction[point_id],
                "HTPB_solid_fraction": htpb_fraction[point_id],
                "species": "AP" if species_ap[point_id] else "HTPB",
                "reference_density_kg_m3": (
                    properties["rho_ap"] if species_ap[point_id]
                    else properties["rho_htpb"]),
                "distance_to_interface_m": distance[point_id],
                "signed_distance_to_interface_m": (
                    distance[point_id] if species_ap[point_id]
                    else -distance[point_id]),
            })
        contour_summary.append((points, species_ap, distance))

    temporary_path = output_path.with_suffix(".csv.tmp")
    with temporary_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    temporary_path.replace(output_path)
    return rows, contour_summary


def write_validation_plot(current, contours, path):
    x = current["x"] * 1.0e3
    y = current["y"] * 1.0e3
    figure, axis = plt.subplots(figsize=(7.0, 8.0))
    image = axis.pcolormesh(
        x, y, current["temperature"].T, shading="nearest", cmap="inferno")
    for points, species_ap, distance in contours:
        axis.scatter(
            1.0e3 * points[:, 0], 1.0e3 * points[:, 1],
            c=np.where(species_ap, 1.0, 0.0), cmap="coolwarm",
            vmin=0.0, vmax=1.0, s=8.0, linewidths=0.0)
        interface = np.isfinite(distance) & (distance < 0.75 * min(
            current["x"][1] - current["x"][0],
            current["y"][1] - current["y"][0]))
        axis.scatter(
            1.0e3 * points[interface, 0], 1.0e3 * points[interface, 1],
            facecolors="none", edgecolors="white", s=35.0, linewidths=0.8)
    axis.set_aspect("equal")
    axis.set_xlabel("x [mm]")
    axis.set_ylabel("y [mm]")
    axis.set_title(f"eta=0.5 samples, t={1.0e3 * current['time']:.4f} ms")
    figure.colorbar(image, ax=axis, label="Temperature [K]")
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(
        description="Extract eta=0.5 regression data from Cyclone plotfiles")
    parser.add_argument(
        "output_directory", nargs="?", default=Path(__file__).resolve().parent,
        type=Path)
    parser.add_argument("--level", type=int, default=None)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--max-outputs", type=int, default=None)
    args = parser.parse_args()

    root = args.output_directory.resolve()
    metadata = metadata_values(root / "metadata")
    level = args.level if args.level is not None else int(metadata["amr.max_level"])
    calorie_per_centimeter = 418.4
    properties = {
        "rho_ap": leading_float(metadata["AP_solid.reference_density"]),
        "rho_htpb": leading_float(metadata["HTPB_solid.reference_density"]),
        "k_ap": leading_float(metadata["AP_solid.thermal_conductivity"]),
        "k_htpb": leading_float(metadata["HTPB_solid.thermal_conductivity"]),
        "gas_lambda_a": calorie_per_centimeter * leading_float(
            metadata["gas.transport.rocfire.lambda_a"]),
        "gas_lambda_b": calorie_per_centimeter * leading_float(
            metadata["gas.transport.rocfire.lambda_b"]),
    }

    plotfiles = sorted(root.glob("*cell"), key=plotfile_time)
    if len(plotfiles) < 2:
        raise RuntimeError("At least two plotfiles are required to compute eta_dot")
    selected = list(range(0, len(plotfiles), args.stride))
    if args.max_outputs is not None:
        selected = selected[:args.max_outputs]

    data_directory = root / "surface_data"
    validation_directory = root / "surface_validation"
    data_directory.mkdir(exist_ok=True)
    validation_directory.mkdir(exist_ok=True)
    validation_indices = {selected[0], selected[len(selected) // 2], selected[-1]}

    manifest_rows = []
    cache = {}
    for sequence, index in enumerate(selected, start=1):
        needed = {index}
        needed.add(index - 1 if index > 0 else index + 1)
        needed.add(index + 1 if index + 1 < len(plotfiles) else index - 1)
        for needed_index in needed:
            if needed_index not in cache:
                cache[needed_index] = load_snapshot(plotfiles[needed_index], level)
        current = cache[index]
        previous = cache.get(index - 1) if index > 0 else None
        following = cache.get(index + 1) if index + 1 < len(plotfiles) else None
        csv_path = data_directory / f"{current['path'].stem}_surface.csv"

        rows, contours = extract_snapshot(
            current, previous, following, properties, csv_path)
        if index in validation_indices:
            write_validation_plot(
                current, contours,
                validation_directory / f"{current['path'].stem}_surface.png")
        finite_distance = np.asarray([
            row["distance_to_interface_m"] for row in rows], dtype=float)
        manifest_rows.append({
            "plotfile": current["path"].name,
            "time_s": current["time"],
            "datafile": str(csv_path.relative_to(root)),
            "contour_count": len(contours),
            "point_count": len(rows),
            "points_with_interface_distance": int(np.count_nonzero(
                np.isfinite(finite_distance))),
            "mean_mass_flux_kg_m2_s": np.mean([
                row["mass_flux_kg_m2_s"] for row in rows]) if rows else math.nan,
            "mean_heat_flux_into_solid_W_m2": np.mean([
                row["heat_flux_into_solid_W_m2"] for row in rows]) if rows else math.nan,
        })
        keep = {index, index + 1}
        cache = {key: value for key, value in cache.items() if key in keep}
        print(
            f"[{sequence}/{len(selected)}] {current['path'].name}: "
            f"{len(contours)} contours, {len(rows)} points", flush=True)

    manifest_path = root / "surface_data_manifest.csv"
    with manifest_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=manifest_rows[0].keys())
        writer.writeheader()
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
