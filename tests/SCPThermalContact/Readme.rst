This test checks the conduction of heat across the diffuse interface between
a solid propellant and the surrounding gas, independent of everything else
in the Flame/Hydro coupling (no combustion, no interface regression, no
chemistry). It exists to validate :code:`Hydro::AdvanceSolidEnergy` against
a known closed-form solution of the heat equation.

.. figure:: ../../../tests/SCPThermalContact/reference/comparison.png
   :align: center

   Final temperature profile (:code:`t=1.6e-8s`) along a ray through the
   solid and gas, compared against the analytical solution below. Regenerate
   with :code:`python3 plot_comparison.py output` after running the input
   file. Points within a few interface widths of :code:`y=0` depart from the
   analytical curve for the reason given in the Notes section below - this
   is expected, not a bug.

**Setup**

A flat, static interface sits at :code:`y=0`: solid (:code:`y<0`) initially
at :code:`T_i=300K`, gas (:code:`y>0`) held near :code:`T_g=1000K`. The phase
field's own kinetics are frozen (:code:`propellant.fullfeedback.bound` set
far above any temperature reached in the run), so the interface never moves
and no mass/energy is injected by combustion - the only thing changing the
temperature field over time is conduction.

**How heat conduction from solid to gas works**

:code:`Hydro::AdvanceSolidEnergy` (in :code:`src/Integrator/Hydro.cpp`)
advances the solid's caloric energy field (:code:`solid.energy_mf`, from
which a solid-only temperature :code:`T_solid = E_solid/(rho_phys*cp) + T_ref`
is recovered) via two physically distinct mechanisms:

1. **Bulk conduction inside the solid**:
   :code:`d/dt(T_solid) = div(phi_s * alpha_solid * grad(T_solid))`,
   ordinary Fourier diffusion with :code:`alpha_solid = k_solid/(rho_phys*cp)`,
   applied only to the solid-only temperature field - never the shared,
   reactive gas/solid blended field. (An earlier version of this operator
   differentiated the blended field directly; at the interface its gradient
   inherits however steep the reacting gas temperature happens to be, and a
   tiny physical diffusivity multiplying an unbounded gradient produced a
   runaway. See the two preceding commits on this file for the fix.)

2. **Interfacial exchange with the gas**: a bounded Robin/convective-type
   flux,

   .. code-block::

      q_interface = h_interface * |grad(phi_s)| * (T_gas - T_solid)

   where :code:`h_interface` (:code:`hydro.solid.h_interface`, [W/m^2/K]) is
   a *physical* interfacial contact conductance, :code:`phi_s` is the local
   solid fraction (from the phase field :code:`eta`), and
   :code:`T_gas`/:code:`T_solid` are the *un-blended* physical temperatures
   on each side, recovered algebraically from the previous step's shared
   temperature field. :code:`|grad(phi_s)|` is largest right at the
   interface and ~zero away from it, so this term is localized to the
   diffuse interface layer.

   :code:`|grad(phi_s)|` (first power, not squared) is what makes this
   converge to a genuine sharp-interface condition as the phase field's
   diffuse width :code:`eps` -> 0: for any monotonic profile running from 1
   to 0, :code:`integral(|grad(phi_s)|) = 1` exactly, for *any* eps (this is
   just the fundamental theorem of calculus) - so :code:`|grad(phi_s)|` is a
   properly normalized surface delta function, and
   :code:`integral(q_interface) -> h_interface*(T_gas(0)-T_solid(0))` as
   eps->0, independent of eps. An earlier version of this term used
   :code:`k_solid * |grad(phi_s)|^2` instead (i.e. an implicit
   :code:`h = k_solid/(eps*sqrt(pi))`, tied to eps) - its integral diverges
   as :code:`1/eps`, so refining eps silently strengthened the coupling
   every time rather than converging to anything fixed. A mesh/eps
   refinement study with that version showed the error growing (not
   shrinking) with refinement - the tell that something was wrong, since
   the discretization itself was fine (see the bulk-only sine-decay
   companion check, which agrees with its analytical solution to <0.02K
   throughout). With :code:`h_interface` fixed and eps/mesh refined
   together, the same study now converges properly - see the git history on
   this file for the before/after numbers.

Away from the interface, only term (1) is active (bulk Fourier diffusion).
Right at the interface, term (2) dominates and acts like a surface heat
transfer coefficient equal to :code:`h_interface` directly - not derived
from :code:`k_solid` or :code:`eps` at all.

**Analytical solution**

Together, these two mechanisms are exactly the classical problem of a
semi-infinite solid, initially at a uniform temperature :code:`T_i`, whose
surface exchanges heat by convection with an ambient fluid at :code:`T_g`
via a heat transfer coefficient :code:`h` (= :code:`h_interface`)
(Carslaw & Jaeger, *Conduction of Heat in Solids*, 2nd ed., section 2.8):

.. code-block::

   T(x,t) = T_i + (T_g-T_i) * [ erfc(x/(2*sqrt(a*t)))
              - exp(h*x/k + h^2*a*t/k^2) * erfc(x/(2*sqrt(a*t)) + h*sqrt(a*t)/k) ]

where :code:`x` is depth into the solid from the interface and
:code:`a = alpha_solid`. :code:`generate_reference.py` evaluates this
formula directly (it does not read any simulation output) and writes
:code:`reference/reference-2d.csv`; :code:`test` compares the simulation's
final :code:`temperature` field along a ray through the solid against it.

**Notes**

- The initial temperature condition (:code:`temp.ic`) uses a *sharp* step at
  :code:`y=0`, not the diffuse :code:`eps` profile used for :code:`eta` -
  the analytical solution assumes a uniformly cold solid at :code:`t=0`, with
  the interface profile developing entirely through the conduction physics
  being tested, not baked into the initial condition.
- The gas pressure (:code:`hydro.density.ic`, :code:`hydro.pressure.ic`) is
  set far above a realistic chamber pressure. This does not change the
  ideal-gas sound speed (:code:`c=sqrt(gamma*R*T)`, independent of pressure),
  so it doesn't affect the timestep - it exists purely to inflate the gas's
  volumetric heat capacity so it behaves as a near-constant-temperature
  reservoir over the (very short, nanosecond-scale) test duration. Without
  it, the gas region would need to be many times deeper (and the grid
  correspondingly larger) to avoid measurably cooling as it feeds the solid.
- Close to the interface (within a few multiples of :code:`eps`), the
  simulation's *diffuse* interface inherently departs somewhat from the
  analytical solution's idealized *infinitely sharp* boundary - this is an
  expected modeling gap, not a solver error, and is why the comparison
  tolerance in :code:`test` is looser than a typical bit-for-bit regression
  check (observed relative error is ~1.8% as of this writing). Beyond
  roughly :code:`4*eps` into the solid the two agree to well under 1K.
- This test is also what caught a real bug in the temperature reconstruction
  used throughout :code:`Hydro.cpp` (not specific to :code:`AdvanceSolidEnergy`):
  for any cell at or below :code:`hydro.cutoff`, :code:`ApplyCutoffToConserved`
  overwrites the conserved energy to the pure-solid value, but a separate
  "C1-continuous" temperature blend still tried to reconstruct a *gas* state
  from that same now-solid-only energy - dividing energy built on the
  *physical* solid density/cp scale by a reconstructed density on the *tame*
  Riemann-blend scale (~195x smaller), producing spurious temperatures of
  tens of thousands of Kelvin that were still visible after being weighted
  by the small (but nonzero) true gas fraction in the final blend. Confirmed
  via a mesh/eta refinement study to get *worse*, not better, with
  refinement - a useful signal that it wasn't just discretization error.
  Fixed by skipping the gas-state reconstruction entirely below
  :code:`hydro.cutoff` (there is no real gas state left to reconstruct
  there) and reporting the solid caloric temperature directly instead.
- The domain and timestep are both deliberately tiny (nanometers, picoseconds)
  purely so the test runs in seconds; see the header comment in
  :code:`generate_reference.py` for how those relate to the physical
  constants being tested. The physics being checked (the discretized
  conduction operator) is scale-invariant, so there is nothing special about
  these particular numbers - they were chosen only to keep the run fast on a
  laptop.
