#!/usr/bin/env python3
"""
Generate the analytical reference data for the LowMachElasticPressure test.

This does NOT read simulation output - it evaluates the closed-form 1-D
confined (oedometric) compression solution directly (see Readme.rst for the
derivation), so `test` checks the solver against a known-correct answer
rather than against a previous simulation's own output.

Physics
-------
A rigid-solid slab occupies y < y0, loaded by a uniform, frozen fluid
pressure P above it. Roller BCs on the x-faces force u_x = 0 everywhere, so
elastic.apply_fluid_pressure's interfacial traction reduces to exact 1-D
confined compression:

    sigma_yy(y) = -P                          (uniform in the solid)
    disp_y(y)   = -P*y/M                      (linear, disp_y(0) = 0)

with M = kappa + 4*mu/3 the confined/oedometric modulus (plane-strain
reduction of Model::Solid::Finite::NeoHookeanPredeformed at F=I - see
Readme.rst). mu/kappa must match `input` (Solid.elastic.model.mu/kappa)
exactly - if you change one there, change it here too and re-run this
script.
"""
import numpy as np
import pandas as pd

# --- Must match input ---
mu = 140.0      # Solid.elastic.model.mu
kappa = 150.0   # Solid.elastic.model.kappa
lam = kappa - 2.0 * mu / 3.0
M = lam + 2.0 * mu

x_ray = 1.6e-3      # x coordinate of the comparison ray in `test`
y_lo, y_hi = 0.5e-3, 2.4e-3
n_pts = 200

pressures = {
    "pressure-off.csv": 0.0,
    "pressure-1Pa.csv": 1.0,
    "pressure-4Pa.csv": 4.0,
}


def analytic(y, P):
    disp_y = -P * y / M
    stress_yy = np.full_like(y, -P)
    return disp_y, stress_yy


if __name__ == "__main__":
    y = np.linspace(y_lo, y_hi, n_pts)
    for filename, P in pressures.items():
        disp_y, stress_yy = analytic(y, P)
        df = pd.DataFrame({
            "x": np.full_like(y, x_ray),
            "y": y,
            "z": np.zeros_like(y),
            "disp_y": disp_y,
            "stress_yy": stress_yy,
        })
        df.to_csv(f"reference/{filename}")
        print(f"Wrote {len(df)} rows to reference/{filename} (P={P})")
