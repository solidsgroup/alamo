Experimental surface-temperature reconstructions
================================================

For an irreversible solid-to-gas ``arrhenius_surface_flux`` mechanism,
``reconstruct_surface_temperature = true`` retains its original default:
temperature interpolated at the total condensed volume fraction contour
``eta = 0.5``. No production input is changed by the additional options.

The following keys are relative to ``<mechanism>.phase_change.``::

    reconstruct_surface_temperature = true
    surface_temperature_reconstruction = conversion_weighted
    surface_temperature_band_cutoff = 0.01
    surface_temperature_contour_fraction = 0.5

``surface_temperature_contour_fraction`` is the total **condensed** volume
fraction, not the gas fraction used by some diffuse-boundary formulations.
It must lie strictly in ``(1e-10,1-1e-10)`` and, for a weighted reconstruction,
inside its averaging band. A nondefault contour requires reconstruction to
be enabled. Changing this value changes the kinetic temperature closure;
it does not merely change a reported surface location. The default remains
0.5. A contour fitted on one grid is not an independently calibrated physical
surface: assess interface-width and grid sensitivity.

Available reconstructions:

* ``contour``: original contour interpolation, unchanged default.
* ``solid_mass_weighted``: average of the local caloric temperature weighted
  by the reacting solid's partial mass density along the interface normal.
  For pure AP with constant intrinsic density, this is also volume-fraction
  weighting; it is not an independent measurement of solid-phase temperature.
* ``conversion_weighted``: average weighted by predicted converted mass per
  unit volume. The predictor executes the existing implicit kinetic solve
  using the contour temperature on scratch density and temperature. All its
  physical changes are discarded. The frozen weights are communicated before
  the single actual mass/heat/volume/recoil update.
* ``arrhenius_equivalent``: uses the same predicted-conversion weights but
  averages ``exp(-Ta/T)`` before inverting for temperature. Stable log-sum-exp
  evaluation avoids exponential overflow. At zero activation temperature the
  arithmetic conversion-weighted average is used.

All weighted averages use the connected normal-line band containing the
contour, bounded by the first exits from ``cutoff < eta < 1-cutoff``.
The cutoff must be in ``(1e-10,0.5)``. Boundary locations use bisection and
multilinear interpolation; trapezoidal quadrature spacing is at most half
the smallest cell width. The input interface thickness sets the search
reach (ten thicknesses from the contour), not an assumed eta profile.
Communication allows the additional distance from a cell to the contour.
MPI, periodic and coarse/fine ghosts are filled for the frozen fields and
the predictor weights. Incomplete bands abort instead of silently using
gas temperature. If predicted conversion is identically zero, its average
is undefined and the contour temperature is retained.

These are experimental kinetic closures, not thermodynamic separation of
solid and gas temperatures, nor a solution of Chen's surface energy balance.
In particular, the conversion weights predict the current step from the
control closure; they are not a converged nonlocal fixed point of the final
conversion. No predictor history is needed on restart.

The chosen reconstructed temperature is the base temperature for kinetics.
The existing local implicit enthalpy-induced increment is still added to
that base. The actual local caloric temperature remains separate and is
updated only by the physical energy transfer. Heat capacities, conductivity,
gas chemistry, geometry, and Arrhenius A/Ea are not modified by selecting a
reconstruction. Timestep, mesh and band-cutoff sensitivity must be assessed
before using a closure as a calibrated physical model.

In the initial pure-AP comparison, solid-mass weighting nearly extinguished
regression at unchanged A/Ea/Q. Arrhenius-equivalent weighting produced
extreme mixed-cell cooling (below 100 K) while those cells' caloric temperatures
still entered gas EOS and chemistry. That diagnostic was excluded from
calibration even though its regression rate increased. These observations
are reasons to retain the contour default; a rate match alone does not
validate a reconstructed temperature or its resulting local thermal state.

Comparison tests and simulation outputs for this experiment are kept outside
the repository under ``/tmp/pure_ap_ea_a_calibration/``.

Constant-cp reference-enthalpy diagnostic
---------------------------------------

An independent opt-in thermal update is available for irreversible
solid-to-gas Arrhenius conversion::

    thermal_closure = reference_enthalpy
    enthalpy_reference_temperature = 300_K

The default ``thermal_closure = local_heat`` retains the prior update,
including its frozen-composition heat capacity. The opt-in model interprets
the specified signed phase heat Qs as the change of reference-state enthalpy.
For positive converted mass density dm and constant phase heat capacities::

    C_new = C_old + (cp_products - cp_solid) * dm
    C_new * (T_new-Tref) = C_old * (T_old-Tref) + Qs * dm

Heat capacities are obtained from the configured materials, not separate
adjustable mechanism parameters. Gas thermodynamics must be ``gross_model``
or ``cpconstant``; NIST and other temperature-dependent caloric models are
rejected. The same temperature increment is included in the implicit kinetic
solve. Final physical T still determines product EOS volume and the thermal
volume change of pre-existing gas. ``phase_change_heat`` remains the explicit
Qs*dm increment; unlike the default, there is no additional incidental
composition-change sensible-enthalpy source.

At a sharp interface this convention gives a phase enthalpy jump
``h_g(Ts)-h_s(Ts) = -Qs + (cp_g-cp_s)*(Ts-Tref)``. Do not simultaneously add
that sensible correction to Qs. Changing the reference temperature while
keeping Qs fixed changes the physical model when phase heat capacities
differ. Use a consistent formation-enthalpy reference in interpreting fits.

This is a conservative single-temperature diagnostic, not independent solid
and gas temperatures and not a nonlocal interface-flux solver. It does not
by itself remove the finite-width overlap of gas chemistry and solid heat
storage. Defaults, material A/Ea and gas kinetics are unchanged.
