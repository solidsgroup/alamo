#!/usr/bin/env python3
"""Render a concrete LowMach input from the fullfeedback calibration template.

Substitutes the operating pressure and the two fit parameters into
``input.lm.ap_monopropellant_fullfeedback.template`` and writes a runnable
input file. The final-product density used by the gas-side initial condition
and the inflow boundary is derived from the pressure as
``P / (GAS_CONSTANT * INFLOW_TEMPERATURE)``, matching the original input's
inflow value (13.401777701139464 at 3 MPa, 700 K).

Example
-------
    python scripts/render_input.py --template input.lm...template \\
        --pressure-mpa 3.54 --pre-exponential 0.0027 \\
        --activation-temperature 3145 --out /tmp/input_P3.54
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Final-product gas constant [J/kg/K] and inflow temperature [K] used to build
# the gas-side density from pressure. These match the constants baked into the
# original input.lm.ap_monopropellant_fullfeedback.
GAS_CONSTANT = 319.787
INFLOW_TEMPERATURE = 700.0


def product_density(pressure_pa: float) -> float:
    return pressure_pa / (GAS_CONSTANT * INFLOW_TEMPERATURE)


def render(template_text: str, pressure_pa: float, pre_exponential: float,
           activation_temperature: float) -> str:
    replacements = {
        "@PRESSURE@": repr(float(pressure_pa)),
        "@PRODUCT_DENSITY@": repr(product_density(pressure_pa)),
        "@PRE_EXPONENTIAL@": repr(float(pre_exponential)),
        "@ACT_TEMP@": repr(float(activation_temperature)),
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
    parser.add_argument("--pre-exponential", type=float, required=True,
                         help="fullfeedback pre_exponential in 1/Pa/s")
    parser.add_argument("--activation-temperature", type=float, required=True,
                         help="fullfeedback activation_temperature in K")
    parser.add_argument("--out", type=Path, required=True,
                         help="path to write the rendered input file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    template_text = args.template.read_text()
    text = render(template_text, args.pressure_mpa * 1.0e6,
                  args.pre_exponential, args.activation_temperature)
    args.out.write_text(text)
    print(f"Wrote {args.out} (P = {args.pressure_mpa} MPa, "
          f"rho_product = {product_density(args.pressure_mpa * 1.0e6):.6f} kg/m^3, "
          f"pre_exponential = {args.pre_exponential}, "
          f"activation_temperature = {args.activation_temperature})")


if __name__ == "__main__":
    main()
