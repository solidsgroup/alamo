#!/usr/bin/env python3
"""
Generate the analytical reference data for the LowMachElasticInclusion test.

This does NOT read simulation output -- it evaluates the closed-form
two-phase circular-inclusion elasticity solution directly (see Readme.rst
for the full Kolosov-Muskhelishvili derivation), so `test` checks the
solver against a known-correct answer rather than against a previous
simulation's own output.

Physics
-------
A circular AP inclusion of radius `a` sits at the origin, embedded in an
infinite matrix (the LowMach "soft void" ersatz-elastic Gas). The domain
boundary is traction-free (see `input`) -- the ONLY loading is
`elastic.apply_fluid_pressure`'s interfacial body force, which is itself
the fluid pressure acting as a traction on the AP/Gas free surface (see
LowMach::UpdateModel). This is the axisymmetric specialization of the
classical Kolosov-Muskhelishvili two-phase circular inhomogeneity
solution: for purely dilatational (equibiaxial) loading there is no
angular harmonic content, so the general complex potentials collapse to
the elementary monopole/dipole radial form
    u_r(r) = A*r + B/r
with B = 0 required in the (bounded) inclusion, and A = 0 required in the
matrix (traction-free as r -> infinity kills the non-decaying "A*r" term
there, leaving only the localized B/r dipole decay).

Plane-strain constitutive law (matching
Model::Solid::Finite::NeoHookeanPredeformed linearized about F=I -- see
tests/LowMachElasticPressure/Readme.rst): with lambda = kappa - 2*mu/3,
    sigma_rr = 2*(lambda+mu)*A - 2*mu*B/r^2
    sigma_tt = 2*(lambda+mu)*A + 2*mu*B/r^2
Define the "2-D areal modulus" k = 2*(lambda+mu) = 2*kappa + (2/3)*mu.

``elastic.apply_fluid_pressure`` does not enforce ordinary traction
continuity at the AP/Gas interface -- the RHS it adds is
``-P*grad(eta)``, i.e. (since P is spatially uniform)
``-grad(P*eta)`` exactly, so ``div(sigma) = -grad(P*eta)`` is
``div(sigma + P*eta*I) = 0``. Integrating across the diffuse interface
(``eta: 1 -> 0`` from inclusion to matrix) shows ``sigma + P*eta*I`` is
continuous, i.e. ``sigma`` itself has a genuine jump:

    sigma_rr(a+) = sigma_rr(a-) + P

this is the correct interfacial condition for a fluid at pressure P
pushing on the solid, not a continuity condition (verified directly
against raw simulation output: a naive continuity-only match
under-predicts the interior stress magnitude by ~90%).

Matching u_r (continuous, A1 = 0 in the matrix per above) and this jump in
sigma_rr at r=a between inclusion (2, B2=0) and matrix (1, A1=0) gives:

    A2 = -P / (k2+2*mu1)
    B1 = a^2*A2                 (= a^2*(A2-A1), A1=0)

Uniform stress inside the inclusion:
    sigma_in = k2*A2

Outside (matrix, r>=a) -- pure dipole decay, no remote offset (A1=0):
    sigma_rr(r) = -2*mu1*B1/r^2
    sigma_tt(r) = +2*mu1*B1/r^2
    u_r(r)      = B1/r

mu1/kappa1 (Gas/void) and mu2/kappa2 (AP) must match `input`
(elastic.void.model.mu/kappa, AP_solid.elastic.model.mu/kappa) exactly -- if
you change one there, change it here too and re-run this script.
"""
import numpy as np
import pandas as pd

# --- Must match input ---
mu1, kappa1 = 5.0e3, 4.0e5   # elastic.void.model.mu/kappa      (matrix / Gas)
mu2, kappa2 = 9.47e6, 1.00e8 # AP_solid.elastic.model.mu/kappa  (inclusion, AP)
a = 1.5e-4                   # inclusion radius R
w = 4.0e-5                   # diffuse interface width (must match `input`)
P = 1.0                      # fluid pressure

k1 = 2.0 * kappa1 + (2.0 / 3.0) * mu1
k2 = 2.0 * kappa2 + (2.0 / 3.0) * mu2

A1 = 0.0                      # traction-free domain boundary
A2 = -P / (k2 + 2.0 * mu1)
B1 = a * a * A2

sigma_in = k2 * A2

x_lo, x_hi = -1.3e-3, 1.3e-3
n_pts = 400


def analytic(x):
    """Evaluate along the y=0 ray. Returns disp_x, disp_y, stress_xx, stress_yy."""
    r = np.abs(x)
    sign = np.sign(x)
    sign[sign == 0.0] = 1.0

    inside = r < a

    u_r = np.where(inside, A2 * r, B1 / np.where(r > 0, r, 1.0))
    sigma_rr = np.where(inside, sigma_in, -2.0 * mu1 * B1 / np.where(r > 0, r * r, 1.0))
    sigma_tt = np.where(inside, sigma_in, 2.0 * mu1 * B1 / np.where(r > 0, r * r, 1.0))

    disp_x = u_r * sign
    disp_y = np.zeros_like(x)
    # Along y=0 (theta=0 or pi): sigma_xx = sigma_rr, sigma_yy = sigma_tt.
    stress_xx = sigma_rr
    stress_yy = sigma_tt
    return disp_x, disp_y, stress_xx, stress_yy


if __name__ == "__main__":
    x = np.linspace(x_lo, x_hi, n_pts)
    disp_x, disp_y, stress_xx, stress_yy = analytic(x)
    df = pd.DataFrame({
        "x": x,
        "y": np.zeros_like(x),
        "z": np.zeros_like(x),
        "disp_x": disp_x,
        "disp_y": disp_y,
        "stress_xx": stress_xx,
        "stress_yy": stress_yy,
    })
    df.to_csv("reference/inclusion.csv")
    print(f"Wrote {len(df)} rows to reference/inclusion.csv "
          f"(sigma_in={sigma_in:.6g}, A1={A1:.6g}, A2={A2:.6g}, B1={B1:.6g})")

    zero = pd.DataFrame({
        "x": x,
        "y": np.zeros_like(x),
        "z": np.zeros_like(x),
        "disp_x": np.zeros_like(x),
        "disp_y": np.zeros_like(x),
        "stress_xx": np.zeros_like(x),
        "stress_yy": np.zeros_like(x),
    })
    zero.to_csv("reference/pressure-off.csv")
    print(f"Wrote {len(zero)} rows to reference/pressure-off.csv (P=0 control)")
