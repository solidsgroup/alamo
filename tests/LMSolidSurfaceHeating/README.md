# Prescribed heat flux on a solid surface

An `interfacial_heat_source` may select `condensed = <rigid solid or liquid>`.
The existing `liquid = <liquid species>` input remains supported with its
original liquid-only validation. Specify exactly one of these keys.

The source deposits `heat_flux * activation(T)` per unit diffuse exposed
area. Its volume density is proportional to
`|eta_g grad(eta_c) - eta_c grad(eta_g)|`; buried condensed interfaces receive
no heat. It changes no species. All inputs and geometric calculations use
the existing device-compatible mechanism.

For a constant prescribed flux, choose an activation temperature and width
that fully activate the source over the simulated temperature range. The
test uses 1 K and 0.01 K, respectively. For reactive flame predictions omit
this source once external heating is no longer part of the physical setup.

This short annulus test checks localization, unchanged mass, and deposited
heat against the independent reference `heat_flux * 2*pi*radius * timestep`.
It provides a prescribed-flux check without a heated gas slab whose thermal
storage can change the flux reaching a moving solid surface.
