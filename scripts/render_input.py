#!/usr/bin/env python3
"""Render a concrete LowMach input from the AP regression calibration template.

Substitutes the operating pressure and the two fit parameters into
``input.lm.ap_monopropellant.template`` (or the _fine variant) and writes a
runnable input file. The gas-side product density is computed inside the
input itself (P / (319.787 * 300), matching tests/LMRFMonoAP/input's
cold-start convention) rather than being precomputed here.

Example
-------
    python scripts/render_input.py --template input.lm...template \\
        --pressure-mpa 3.54 --rate-multiplier 2450 \\
        --activation-temperature 3145 --out /tmp/input_P3.54
"""

from __future__ import annotations

import argparse
from pathlib import Path


def render(template_text: str, pressure_pa: float, rate_multiplier: float,
           activation_temperature: float) -> str:
    replacements = {
        "@PRESSURE@": repr(float(pressure_pa)),
        "@RATE_MULTIPLIER@": repr(float(rate_multiplier)),
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
    parser.add_argument("--rate-multiplier", type=float, required=True,
                         help="phase_change.rate_multiplier (dimensionless)")
    parser.add_argument("--activation-temperature", type=float, required=True,
                         help="phase_change.activation_temperature in K")
    parser.add_argument("--out", type=Path, required=True,
                         help="path to write the rendered input file")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    template_text = args.template.read_text()
    text = render(template_text, args.pressure_mpa * 1.0e6,
                  args.rate_multiplier, args.activation_temperature)
    args.out.write_text(text)
    print(f"Wrote {args.out} (P = {args.pressure_mpa} MPa, "
          f"rate_multiplier = {args.rate_multiplier}, "
          f"activation_temperature = {args.activation_temperature})")


if __name__ == "__main__":
    main()
