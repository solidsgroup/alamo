#!/usr/bin/env python3
"""Compute effective 2-D and plane-strain-equivalent 3-D RVE properties.

The input is an AMReX/BoxLib cell plotfile readable by yt, for example:

    python3 effective_rve_properties.py \
        results_caleb/Homogenisation/output.test_homogenisation/00001cell

The script volume-weights the domain averages of ``stress_xx`` and
``stress_yy``.  For a 2-D dataset yt's ``cell_volume`` is the cell area
times its unit out-of-plane thickness, so this is equivalently an
area-weighted average.

For a periodic RVE, the plotted displacement and ``strain`` fields describe
only the periodic fluctuation.  Their domain-average deformation gradient is
therefore the identity, even when a macroscopic deformation is applied
through ``F0``.  The macroscopic strain must instead be supplied from the
prescribed deformation gradient Fbar.  The defaults below correspond to the
periodic 1% x-stretch in ``input.circle_compression``:
Fbar = diag(1.01, 1.0).
"""

from __future__ import annotations

import argparse

import yt


FIELD_TYPE = "boxlib"
VOLUME_FIELD = ("index", "cell_volume")


def domain_average(data, field_name):
    """Return integral(field dV) / integral(dV) over all AMR cells."""
    field = (FIELD_TYPE, field_name)
    volume = data[VOLUME_FIELD]
    return (data[field] * volume).sum() / volume.sum()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Volume/area-average normal stress and strain from a yt-readable "
            "plotfile and calculate effective RVE elastic properties."
        )
    )
    parser.add_argument(
        "plotfile",
        help="cell plotfile to analyze (e.g. output.test_homogenisation/00001cell)",
    )
    parser.add_argument(
        "--fbar-xx",
        "--fbarxx",
        dest="fbar_xx",
        type=float,
        default=1.01,
        help=(
            "prescribed macroscopic Fbar_xx (default: 1.01, matching "
            "input.circle_compression)"
        ),
    )
    parser.add_argument(
        "--fbar-yy",
        "--fbaryy",
        dest="fbar_yy",
        type=float,
        default=1.0,
        help="prescribed macroscopic Fbar_yy (default: 1.0)",
    )
    args = parser.parse_args()

    if args.fbar_xx <= 0.0 or args.fbar_yy <= 0.0:
        parser.error("Fbar diagonal components must be positive.")

    ds = yt.load(args.plotfile)
    data = ds.all_data()

    # These are <sigma_xx> and <sigma_yy>, where <> denotes the volume (or,
    # in 2-D, area) average.  For finite.neohookeanpre, the reported stress
    # already includes the imposed Fbar through F0 = Fbar^{-1}.
    stress_x = domain_average(data, "stress_xx")
    stress_y = domain_average(data, "stress_yy")

    # Do not average the plotfile strain here: under periodic displacement
    # boundary conditions it is F_fluc = I + grad(u_fluc), not F_physical.
    # The engineering strain used for these secant effective properties is
    # obtained from the imposed macroscopic deformation gradient.
    strain_x = args.fbar_xx - 1.0
    strain_y = args.fbar_yy - 1.0

    denominator = stress_y * strain_y - stress_x * strain_x
    if denominator == 0:
        raise ZeroDivisionError(
            "stress_y*strain_y - stress_x*strain_x is zero; "
            "effective properties are undefined."
        )

    e_eff_2d = (stress_y**2 - stress_x**2) / denominator
    nu_eff_2d = (stress_x * strain_y - stress_y * strain_x) / denominator
    nu_eff_3d = nu_eff_2d / (1.0 + nu_eff_2d)
    e_eff_3d = e_eff_2d * (1.0 - nu_eff_3d**2)

    print(f"plotfile: {args.plotfile}")
    print("prescribed macroscopic deformation:")
    print(f"  Fbar_xx  = {args.fbar_xx}")
    print(f"  Fbar_yy  = {args.fbar_yy}")
    print("volume/area-weighted domain averages:")
    print(f"  stress_x  = {stress_x}")
    print(f"  stress_y  = {stress_y}")
    print(f"  strain_x  = {strain_x}")
    print(f"  strain_y  = {strain_y}")
    print("effective properties:")
    print(f"  E_eff_2D  = {e_eff_2d}")
    print(f"  nu_eff_2D = {nu_eff_2d}")
    print(f"  E_eff_3D  = {e_eff_3d}")
    print(f"  nu_eff_3D = {nu_eff_3d}")


if __name__ == "__main__":
    main()
