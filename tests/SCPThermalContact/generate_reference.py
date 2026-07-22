#!/usr/bin/env python3
"""
Generate the analytical reference data for the SCPThermalContact test.

This does NOT read simulation output - it evaluates a closed-form solution
of the heat equation directly, so the regression test in `test` is checked
against a known-correct answer rather than against a previous simulation's
own output (which would only catch changes in behavior, not existing bugs).

Physics
-------
Hydro::AdvanceSolidEnergy conducts heat across the diffuse solid/gas
interface with two terms (see src/Integrator/Hydro.cpp and Readme.rst):

  1. Bulk conduction inside the solid:      d/dt(T) = alpha * d^2T/dx^2
  2. Interfacial exchange with the gas:     q = h_interface*|grad(phi_s)|*(Tg-Ts)

Away from the interface only (1) applies (grad(phi_s)=0 in the bulk); at the
interface (2) acts like a surface heat-transfer coefficient h_interface
(hydro.solid.h_interface [W/m^2/K]) - a physical contact conductance,
independent of eps. |grad(phi_s)| (first power, not squared) is a proper
surface delta function - its integral across the interface is exactly 1 for
any eps - so as eps->0 this term converges to a genuine sharp-interface Robin
condition set by h_interface alone. Together this is exactly the classical
problem of a semi-infinite solid, initially at a uniform temperature T_i,
whose surface (x=0) exchanges heat by convection with an ambient fluid at
T_g via coefficient h - solved in closed form in Carslaw & Jaeger,
"Conduction of Heat in Solids", 2nd ed., section 2.8:

  T(x,t) = T_i + (T_g-T_i) * [ erfc(x/(2*sqrt(a*t)))
             - exp(h*x/k + h^2*a*t/k^2) * erfc(x/(2*sqrt(a*t)) + h*sqrt(a*t)/k) ]

where x is depth into the solid measured from the interface and
a = alpha_solid = k_solid/(rho_phys*cp_solid). Unlike the earlier version of
this test, h here does NOT depend on eps at all - it is set directly by
hydro.solid.h_interface, so the comparison is exact (up to the diffuse-vs-
sharp-interface gap noted in Readme.rst) for any eps, and should genuinely
converge as eps->0, not just at whatever eps happened to be tuned.

All the physical constants below (k_ap, rho_ap, cp_ap, h_interface, T_i,
T_g, stop_time) must match `input` exactly - if you change one there,
change it here too and re-run this script.
"""
import numpy as np
import pandas as pd
from scipy.special import erfc, erfcx

# --- Must match input ---
k_solid      = 0.4186          # propellant.fullfeedback.k_ap [W/m/K]
rho_phys     = 1950.0          # propellant.fullfeedback.rho_ap [kg/m^3]
cp_solid     = 1297.90         # propellant.fullfeedback.cp_ap [J/kg/K]
h            = 4723395.6       # hydro.solid.h_interface [W/m^2/K] - NOT eps-dependent
Ti           = 300.0           # temp.ic.expression.constant.Ti [K]
Tg           = 1000.0          # temp.ic.expression.constant.Tg [K]
stop_time    = 1.6e-8          # stop_time [s]
x_ray        = 2.0e-8          # x coordinate of the comparison ray in `test`
n_cell_y     = 100
y_lo, y_hi   = -5.0e-7, 5.0e-7

alpha = k_solid / (rho_phys * cp_solid)


def T_robin(depth, t):
    """Semi-infinite solid with surface convection (Carslaw & Jaeger 2.8)."""
    depth = np.asarray(depth, dtype=float)
    T = np.full_like(depth, Ti)
    pos = depth[depth > 0]
    if t > 0 and pos.size:
        s = pos / (2.0 * np.sqrt(alpha * t))
        arg2 = s + h * np.sqrt(alpha * t) / k_solid
        expo = h * pos / k_solid + (h**2) * alpha * t / (k_solid**2)
        # erfcx(z) = exp(z^2)*erfc(z), used to avoid overflow in exp(expo)*erfc(arg2)
        val = erfc(s) - np.exp(expo - arg2**2) * erfcx(arg2)
        T[depth > 0] = Ti + (Tg - Ti) * val
    return T


if __name__ == "__main__":
    dy = (y_hi - y_lo) / n_cell_y
    y = y_lo + (np.arange(n_cell_y) + 0.5) * dy
    y_solid = y[y < 0.0]
    depth = -y_solid

    T = T_robin(depth, stop_time)

    df = pd.DataFrame({
        "x": np.full_like(y_solid, x_ray),
        "y": y_solid,
        "z": np.zeros_like(y_solid),
        "temperature": T,
    })
    df.to_csv("reference/reference-2d.csv")
    print(f"Wrote {len(df)} rows to reference/reference-2d.csv")
    print(f"alpha={alpha:.6e} m^2/s  h={h:.6e} W/m^2/K")
    print(f"T at interface (y->0-): {T[-1]:.3f} K")
