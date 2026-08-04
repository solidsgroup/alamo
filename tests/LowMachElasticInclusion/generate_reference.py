#!/usr/bin/env python3
"""
Generate the analytical reference data for the LowMachElasticInclusion test.

This does NOT read simulation output -- it evaluates the closed-form
solution directly (see Readme.rst for the full derivation), so `test`
checks the solver against a known-correct answer rather than against a
previous simulation's own output.

Physics
-------
A circular AP inclusion of radius `a` sits at the origin, embedded in a
square HTPB block that fills the domain width and extends from y=-1.5e-3 to
y=ytop=1.5e-3; a Gas layer above y=ytop delivers `elastic.apply_fluid_
pressure`'s interfacial body force onto the top face of the block, pushing
straight down. Side and bottom faces of the domain (which coincide with the
HTPB block's own side/bottom faces) are rollers; the top (Gas's free
surface) is traction-free.

Step 1: far-field state in the HTPB block. Ignoring the AP inclusion, the
block/roller/top-load system is exactly 1-D (uniform in x, by the symmetry
of a full-width block with rollers on both sides): u_x = 0 everywhere,
eps_xx = eps_zz = 0 (plane strain), so this is confined ("oedometer")
compression, not free uniaxial stress. The Gas layer above is, by the same
1-D argument, itself statically determinate (traction-free top + roller
sides + uniform interfacial load below) with sigma_yy = 0 throughout,
independent of Gas's own modulus (see `input`) -- so the HTPB/Gas
interfacial pressure jump (same mechanism as the original circular-
inclusion-in-void version of this test: sigma_nn(fluid side) =
sigma_nn(solid side) + P, going from solid to fluid in the direction of the
outward normal) gives sigma_yy(HTPB top) = -P directly. Plane-strain
constitutive law (lambda = kappa - 2*mu/3) then gives, uniformly through
the HTPB block:

    sigma_yy_inf = -P
    sigma_xx_inf = -P * lambda1 / (lambda1 + 2*mu1)

-- NOT equibiaxial (that only happens in the incompressible limit
lambda1 -> infinity; HTPB's nu ~ 0.499 puts this test close to, but not
exactly at, that limit).

Step 2: superpose the classical two-phase circular-inhomogeneity-under-
remote-stress solution (AP embedded in "infinite" HTPB, valid since the
block is 10x the inclusion radius on every side -- Saint-Venant/image
corrections are O((a/L)^2), negligible next to the ~5-8% error this test
tolerates). AP/HTPB is an ordinary bonded interface here (ordinary
continuity of traction and displacement -- no fluid-pressure jump; that
only applies at the Gas/HTPB interface). Decompose the remote state into:

  - an isotropic part p0 = (sigma_xx_inf+sigma_yy_inf)/2, handled by the
    SAME axisymmetric Lame solution the original version of this test used
    (u_r = A*r+B/r), but now sourced by a REMOTE stress rather than an
    interfacial jump, so A1 != 0 in the matrix:

        A1 = p0/k1
        A2 = p0*(k1+2*mu1) / (k1*(k2+2*mu1))
        B1 = a^2*(A2 - p0/k1)
        sigma_in_iso = k2*A2                                    (r<a)
        sigma_rr(r) = p0 - 2*mu1*B1/r^2, sigma_tt(r) = p0 + 2*mu1*B1/r^2  (r>=a)

    where k = 2*kappa + (2/3)*mu (1=HTPB matrix, 2=AP inclusion), same
    definition as before.

  - a deviatoric part s = (sigma_xx_inf-sigma_yy_inf)/2 (equivalent to a
    remote pure-shear state at 45 degrees), handled by the classical
    circular-inhomogeneity-under-remote-shear solution via Kolosov-
    Muskhelishvili complex potentials phi(z), psi(z). Solving the bonded-
    interface matching problem (see derivation notes in this repo's commit
    history / Readme.rst) gives, with kM1 = 3-4*nu1 (Muskhelishvili's plane-
    strain material constant for the MATRIX only -- remarkably, for a
    CIRCULAR inhomogeneity under remote shear, the inclusion's own kM2
    drops out of the solution entirely):

        B  = s*a^2*(mu1-mu2) / (mu1 + kM1*mu2)
        D3 = a^2*B
        gamma2p = -s*mu2*(1+kM1) / (mu1 + kM1*mu2)

    Interior (r<a, uniform Cartesian stress -- the classic 2-D "Eshelby"
    result that a circular inhomogeneity's interior field under remote
    uniform stress is itself uniform):
        sigma_xx_dev_in = -gamma2p,  sigma_yy_dev_in = +gamma2p

    Along the y=0 ray (x=r, theta=0/pi), exterior (r=|x|>=a):
        sigma_xx_dev(x) = -4*B/x^2 + s + 3*D3/x^4
        sigma_yy_dev(x) =        -s - 3*D3/x^4
    (even in x, so this holds for x<0 too without extra sign handling)

    Displacement (2*mu*(ux+i*uy) = kM*phi(z) - z*conj(phi'(z)) - conj(psi(z))),
    evaluated on the real axis (uy_dev = 0 there by symmetry, matching the
    isotropic part):
        ux_dev(x) = [(kM1+1)*B/x + s*x - D3/x^3] / (2*mu1)      (r>=a, odd in x)
        ux_dev(x) = -gamma2p*x / (2*mu2)                        (r<a, odd in x)

Total: sigma_xx = sigma_xx_iso + sigma_xx_dev, etc.; ux = ux_iso + ux_dev,
with ux_iso(x) = A1*x + B1/x (r>=a), A2*x (r<a) -- same odd-in-x form as the
original test (u_r = A*r+B/r with disp_x = u_r*sign(x) collapses to this
single formula in x). uy = 0 exactly along y=0 for both parts.

mu1/kappa1 (HTPB) and mu2/kappa2 (AP) must match `input` exactly
(HTPB_solid.elastic.model.mu/kappa, AP_solid.elastic.model.mu/kappa) -- if
you change one there, change it here too and re-run this script.
"""
import numpy as np
import pandas as pd

# --- Must match input ---
mu1, kappa1 = 1.67e6, 8.33e8  # HTPB_solid.elastic.model.mu/kappa (matrix)
mu2, kappa2 = 9.47e6, 1.00e8  # AP_solid.elastic.model.mu/kappa  (inclusion)
a = 1.5e-4                    # inclusion radius R
w = 4.0e-5                    # diffuse interface width (must match `input`)
P = 1.0                       # fluid pressure

lambda1 = kappa1 - (2.0 / 3.0) * mu1
k1 = 2.0 * kappa1 + (2.0 / 3.0) * mu1
k2 = 2.0 * kappa2 + (2.0 / 3.0) * mu2

nu1 = (3.0 * kappa1 - 2.0 * mu1) / (2.0 * (3.0 * kappa1 + mu1))
kM1 = 3.0 - 4.0 * nu1

# --- Step 1: far-field confined-compression state in the HTPB block ---
sigma_yy_inf = -P
sigma_xx_inf = -P * lambda1 / (lambda1 + 2.0 * mu1)
p0 = 0.5 * (sigma_xx_inf + sigma_yy_inf)
s = 0.5 * (sigma_xx_inf - sigma_yy_inf)

# --- Step 2a: isotropic (mean-stress) inhomogeneity correction ---
A1 = p0 / k1
A2 = p0 * (k1 + 2.0 * mu1) / (k1 * (k2 + 2.0 * mu1))
B1 = a * a * (A2 - p0 / k1)
sigma_in_iso = k2 * A2

# --- Step 2b: deviatoric (remote-shear) inhomogeneity correction ---
Bdev = s * a * a * (mu1 - mu2) / (mu1 + kM1 * mu2)
D3dev = a * a * Bdev
gamma2p = -s * mu2 * (1.0 + kM1) / (mu1 + kM1 * mu2)

sigma_xx_in = sigma_in_iso - gamma2p
sigma_yy_in = sigma_in_iso + gamma2p

x_lo, x_hi = -1.3e-3, 1.3e-3
n_pts = 400


def analytic(x):
    """Evaluate along the y=0 ray. Returns disp_x, disp_y, stress_xx, stress_yy."""
    r = np.abs(x)
    inside = r < a
    r_safe = np.where(r > 0.0, r, 1.0)
    x_safe = np.where(x != 0.0, x, 1.0)

    # Isotropic part (odd-in-x u_r*sign(x) form collapses to a single
    # formula in x -- see module docstring).
    ux_iso = np.where(inside, A2 * x, A1 * x + B1 / x_safe)
    sigma_rr_iso = np.where(inside, sigma_in_iso, p0 - 2.0 * mu1 * B1 / r_safe**2)
    sigma_tt_iso = np.where(inside, sigma_in_iso, p0 + 2.0 * mu1 * B1 / r_safe**2)

    # Deviatoric part.
    ux_dev = np.where(
        inside,
        -gamma2p * x / (2.0 * mu2),
        ((kM1 + 1.0) * Bdev / x_safe + s * x - D3dev / x_safe**3) / (2.0 * mu1),
    )
    sigma_xx_dev = np.where(inside, -gamma2p, -4.0 * Bdev / r_safe**2 + s + 3.0 * D3dev / r_safe**4)
    sigma_yy_dev = np.where(inside, gamma2p, -s - 3.0 * D3dev / r_safe**4)

    # Along y=0: sigma_xx=sigma_rr, sigma_yy=sigma_tt for the isotropic part.
    disp_x = ux_iso + ux_dev
    disp_y = np.zeros_like(x)
    stress_xx = sigma_rr_iso + sigma_xx_dev
    stress_yy = sigma_tt_iso + sigma_yy_dev
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
          f"(sigma_xx_in={sigma_xx_in:.6g}, sigma_yy_in={sigma_yy_in:.6g}, "
          f"p0={p0:.6g}, s={s:.6g})")

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
