NIST solid-AP enthalpy (opt-in)
========================================

LowMach provides an experimental solid ammonium-perchlorate caloric model::

   AP_solid.thermal_model = nist_ap_shomate
   AP_solid.nist_ap.transition_width = 2_K

The flame diagnostics use these stricter nonlinear-conduction controls::

   diffusion.enthalpy.max_iterations = 150
   diffusion.enthalpy.relative_tolerance = 1e-9
   diffusion.enthalpy.absolute_tolerance = 1e-8_K

The default is ``constant_cp``; existing input decks keep their original
behavior. ``AP_solid.specific_heat`` remains required by the material parser,
but NIST cp(T) replaces that constant in the active thermal calculations.
Specify the solid Arrhenius prefactor and activation energy explicitly.
The constant-cp Chen reference helpers are not supported with this model.

Thermodynamics
--------------

The two Shomate fits are from the `NIST Chemistry WebBook AP table
<https://webbook.nist.gov/cgi/cbook.cgi?ID=C7790989&Type=JANAFS&Table=on#JANAFS>`_.
Molecular weight is 117.489 g/mol. Molar heat capacity and enthalpy are
converted to J/(kg K) and J/kg, with h(300 K)=0. The solid-solid transition
at 513.15 K is included through the enthalpy difference between the fits.

A cubic Hermite enthalpy bridge across ``transition_width`` preserves both
endpoint enthalpies and endpoint heat capacities. The active cp is the
analytic derivative of that same enthalpy. This is a regularized equilibrium
transition, not a kinetic or hysteretic polymorph model. Check sensitivity
to the width (allowed range 0.1--20 K), timestep, mesh, and diffuse-interface
thickness before interpreting a calibrated parameter.

The source fits span 298--1500 K. Constant-cp linear enthalpy tails guard
numerical trial states outside this range; they do not validate material
properties outside it. Audit temperatures and solid mass outside the range.
The high-temperature solid fit does not establish that AP remains
undecomposed there.

Coupling and conventions
------------------------

Mixture heat capacity uses solid partial density times NIST cp(T).
Implicit conduction iterates in integrated enthalpy; gas chemistry integrates
an enthalpy coordinate; kinetic phase transfer uses the same nonlinear
heat-to-temperature inverse. Gas species retain their common constant
Gross-model cp. Explicit transport retains the existing temperature-based
RK formulation and therefore requires timestep-convergence checks.

Solid Q retains the existing local phase-heat convention at frozen
pre-transfer composition. It is not redefined as a standard formation
enthalpy. When auditing a global sensible-enthalpy budget, include the
material enthalpy difference associated with phase conversion. For the
NIST solid, use integrated h(T), not cp(T)*(T-300). A Q fitted to an outlet
temperature is an effective model parameter, not a NIST-derived
decomposition enthalpy.

The pyrolysis surface-temperature closure is unchanged: the reconstructed
eta=0.5 temperature receives the local implicit heat correction. The legacy
``energy`` output is not a NIST total-enthalpy diagnostic; reconstruct the
mixture enthalpy from partial densities and the caloric functions.

Supported diagnostic scope
--------------------------

The current opt-in implementation requires CPU, stationary prescribed rigid
AP solid, Gross constant-cp gas thermodynamics, and Gross or frozen chemistry.
It rejects GPU, moving or deformable solids, liquids, other equilibrium
phase-change mechanisms, AP temperature overrides, homogeneous-binder
substitution, and constant-cp Chen reference helpers.

Constitutive inversion, closed conduction, heating through the transition,
phase transfer, reacting-mixture enthalpy, and disabled-model compatibility
are tested separately from flame calibration. The calibration inputs,
test drivers, and raw outputs are maintained outside the repository.
