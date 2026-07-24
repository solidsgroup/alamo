#!/usr/bin/env python3
"""Render a concrete LowMach input from the AP/HTPB sandwich calibration template.

Substitutes the operating pressure and the two HTPB fullfeedback fit
parameters into ``input.lm.ap_htpb_fullfeedback.template`` and writes a
runnable input file. AP_decomposition's fullfeedback parameters are baked
into the template already (held fixed at the values calibrated against pure
AP monopropellant regression-rate data) and are not touched here. The
gas-side product density is derived from pressure the same way as the
original ``input.lm.ap_htpb_fullfeedback``: ``P / (GAS_CONSTANT * 700 K)``.

Example
-------
    python scripts/render_htpb_sandwich_input.py \\
        --template input.lm.ap_htpb_fullfeedback.template \\
        --pressure-mpa 1.5 --htpb-pre-exponential 0.001 \\
        --htpb-activation-temperature 4000 --out /tmp/input_P1.5
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Gas constant [J/kg/K] and inflow temperature [K] used to build the gas-side
# density from pressure -- matches the constants baked into the original
# input.lm.ap_htpb_fullfeedback (700 K hot-gas bootstrap).
GAS_CONSTANT = 319.787
INFLOW_TEMPERATURE = 700.0


def product_density(pressure_pa: float) -> float:
    return pressure_pa / (GAS_CONSTANT * INFLOW_TEMPERATURE)


def render(template_text: str, pressure_pa: float, htpb_pre_exponential: float,
           htpb_activation_temperature: float) -> str:
    replacements = {
        "@PRESSURE@": repr(float(pressure_pa)),
        "@PRODUCT_DENSITY@": repr(product_density(pressure_pa)),
        "@HTPB_PRE_EXPONENTIAL@": repr(float(htpb_pre_exponential)),
        "@HTPB_ACT_TEMP@": repr(float(htpb_activation_temperature)),
    }
    text = template_text
    for key, value in replacements.items():
        text = text.replace(key, value)
    leftover = [k for k in replacements if k in text]
    if leftover:
        raise RuntimeError(f"template still contains placeholders after render: {leftover}")
    return text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--template", type=Path, required=True,
                         help="path to the .template input file")
    parser.add_argument("--pressure-mpa", type=float, required=True,
                         help="operating pressure in MPa")
    parser.add_argument("--htpb-pre-exponential", type=float, required=True,
                         help="HTPB_pyrolysis fullfeedback pre_exponential in 1/Pa/s")
    parser.add_argument("--htpb-activation-temperature", type=float, required=True,
                         help="HTPB_pyrolysis fullfeedback activation_temperature in K")
    parser.add_argument("--out", type=Path, required=True,
                         help="path to write the rendered input file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    template_text = args.template.read_text()
    text = render(template_text, args.pressure_mpa * 1.0e6,
                  args.htpb_pre_exponential, args.htpb_activation_temperature)
    args.out.write_text(text)
    print(f"Wrote {args.out} (P = {args.pressure_mpa} MPa, "
          f"rho_product = {product_density(args.pressure_mpa * 1.0e6):.6f} kg/m^3, "
          f"htpb_pre_exponential = {args.htpb_pre_exponential}, "
          f"htpb_activation_temperature = {args.htpb_activation_temperature})")


if __name__ == "__main__":
    main()
