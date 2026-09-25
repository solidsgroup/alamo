Spherical nanoindentation proof of concept
=========================================

This example uses Integrator::Mechanics and Model::Solid::Finite::CrystalPlastic,
not the periodic, linear MechanicsFFT solver. The crystal has its cubic axes
aligned with x/y/z, and the indenter loads the [001] surface.

Geometry: a 4 x 4 x 2 um block, 4.0 um rigid spherical radius, and 0.05 um
maximum penetration. The sides and bottom are clamped, with a free top outside
contact. Load to maximum depth by 0.1 s, hold to 0.15 s, then retract until the
tip is 0.01 um above the original surface at 0.25 s.

This is SMALL-SLOPE contact: vertical nonpenetration with zero lateral traction,
evaluated at reference x/y coordinates. It is not full finite-sliding contact
normal to the sphere. Separation and the contact footprint are solved, not
prescribed from the undeformed sphere intersection. CPU, nonperiodic geometry,
and the upper face of the last coordinate direction are supported. AMR contact
requires ``amr.nsubsteps=1``. The expression BC tags every upper-surface cell
touching ``contact.region`` and optionally tags a subsurface volume with
``contact.refinement_region``. This ensures that the finest level covers every
node inspected by the contact update. Contact cannot include clamped domain
edges.

Expression boundary conditions
-----------------------------

The existing type keys accept expressions, in addition to literal disp/trac/
neumann/periodic names. For example::

    bc.expression.type.zhi = trac trac "if(x*x+y*y<0.04,disp,trac)"

An expression must return disp, trac, or neumann, and may depend on x,y,z,t.
Spatial type expressions currently require nonperiodic geometry. Edge/corner
keys remain separate, as for existing Expression boundaries.

With contact enabled, both type and value expressions additionally receive
``contact`` (the frozen active-set flag) and ``obstacle`` (the upper admissible
vertical displacement). The demo uses::

    bc.expression.type.zhi = trac trac "if(contact,disp,trac)"
    bc.expression.val.zhi = "0" "0" "if(contact,obstacle,0)"

``contact.obstacle``, ``contact.region``, and ``contact.refinement_region`` are
expressions in x,y,z,t. The refinement expression is evaluated at cell centers;
it defaults to zero because contact-region coverage at the upper surface is
always enforced separately.
The contact update uses p + k*(u_vertical-obstacle), with p=-P_zz (or -P_yy
in 2D), and freezes the resulting mask during each Newton/multigrid solve.
``contact.stiffness`` is k in stress/length units: an active-set selection
parameter, NOT a penalty compliance. ``contact.tolerance`` is a stress-valued
hysteresis band. The contact solve aborts if the mask does not settle within
``contact.max_iter`` equilibrium solves. Plasticity advances only after this
outer iteration completes.

``traction_scale`` multiplies both sides of traction equations, including the
tangent, to improve numerical conditioning. Inputs/outputs retain physical
stress units. Its default is one, preserving existing cases.

Material and interpretation
---------------------------

Internal units are um, MPa, s (system.mass=kg). Elastic constants are
C11=168000, C12=121000, C44=75000 MPa; slip resistance is 1000 MPa,
reference slip rate 0.1/s, and rate exponent 1. The linear rate law is chosen
to produce visible hysteresis while keeping the explicit plastic update bounded.
These are illustrative FCC parameters, NOT a calibrated material data set. The
current model has fixed slip resistance, explicit time integration, and no
intrinsic length scale. Do not claim quantitative hardness or indentation-size
effects from this demo.

The input uses a 16 x 16 x 8 base grid with two contact-driven refinement levels:
250 nm in the far field, 125 nm at level one, and 62.5 nm beneath the tip. With
the current refinement expression, startup produces about 20,480 cells across
levels. The expected peak geometric contact diameter spans about ten finest
cells. This is a fast demonstration case, not a demonstrated mesh-converged
calculation. Check timestep, grid resolution, and specimen size before using
quantitative results.

Output
------

Plotfiles include displacement, stress (first Piola for this model), deformation
gradient in the strain fields, and the twelve signed slip variables. Signed
slip is not accumulated absolute slip. Contact log records include actual solve
time, mask stability, active-node count, compressive reaction, maximum
penetration and maximum tensile normal traction. In 3D the reaction unit is
MPa*um^2 = microNewtons; in 2D it is force per out-of-plane length.

Use the contact log for load-depth curves, not thermo.dat: the latter has
known parallel limitations and does not integrate the z face in 3D.
The integrator writes plots after advancing the timestep; the logged contact
time is the time at which the mechanical fields were equilibrated.

The tests/ElasticTypeExpression and tests/ElasticContact cases provide small
manufactured mixed-boundary and elastic load/release checks, respectively.
The 2D contact regression represents a cylinder, not a sphere.

Running
-------

For the default adaptive 3D mesh, use a fresh output directory::

    mpiexec -n 8 bin/mechanics-3d-clang++ tests/Nanoindentation/input \
      solver.verbose=1 \
      plot_file=tests/Nanoindentation/output_amr

Keep the existing Newton and linear tolerances. The default mesh has been checked
through its first loaded step; a complete load/unload cycle on the revised
geometry remains unverified.

From the repository root, a coarse 3D plastic smoke run is::

    mpiexec -n 2 bin/mechanics-3d-clang++ tests/Nanoindentation/input \
      'amr.n_cell=8 8 4' \
      solver.max_iter=300 solver.verbose=1 \
      timestep=0.005 stop_time=0.251 \
      plot_file=tests/Nanoindentation/output_amr_smoke

Set ``model1.gammadot0=0`` for the otherwise identical elastic control. The
default input uses a finer mesh and costs more than this smoke test. Summarize
a completed run with::

    python3 tests/Nanoindentation/summarize.py tests/Nanoindentation/output_amr_smoke

This produces ``load_depth.csv``, ``load_depth.png``, ``slip_section.png``, and ``summary.json`` in
the run directory. The summarizer's tip-depth formula matches this input's
loading schedule; change it if the obstacle motion changes.

The example explicitly enables the conservative stress-divergence operator
and refreshes its tangent coefficients at every Newton iteration. Traction rows
are scaled by 1e-5. These settings matter for consistent residual/tangent
discretization and mixed displacement/traction conditioning. Newton updates
are undamped; convergence failure aborts instead of silently accepting a step.

Verified smoke results
----------------------

The current radius-4.0 viscoplastic settings completed a full coarse
8 x 8 x 4 base-grid cycle with dt=0.005 s. The reaction peaked at 4897.87 uN,
relaxed to 4181.39 uN during the fixed-depth hold, and lost contact between
14 nm and 11 nm on withdrawal. The load-depth loop dissipated 41.5 percent of
the loading work. This verifies the intended hysteresis qualitatively, not mesh
convergence or material calibration.

The results below predate the radius-4.0 contact-driven AMR default and are
retained as historical solver checks for the original radius-0.5 setup.

``output_cp_verified`` completed 51 equilibrium solves on 32 x 32 x 16 cells
using dt=0.005 s, the same rate prefactor (0.001/s), and smoother
relaxation (0.5). The matching-mesh elastic control is ``output_elastic32``.
Both completely separated after withdrawal. The plastic run's maximum
reported penetration was 3.47e-18 um, maximum tensile normal traction was
7.90e-10 MPa, final reaction was 2.85e-13 uN, and final maximum absolute
cell-averaged signed slip was 0.0012274. The contact footprint reached only
nine nodes: these are machinery checks, not a resolved indentation benchmark.

The coarse plastic run shows a small reaction increase during the hold;
constitutive/discretization and timestep validation remain necessary before
interpreting that response physically. Do not present the curve as calibrated
hardness or a validated relaxation prediction.

The earlier ``output_cp_smoke`` run used a ten-times larger rate prefactor
and failed during the hold; it is incomplete and is not a verified result.
The revised adaptive default has not yet been run through a complete cycle.

The serial/MPI mixed-boundary/contact regression cases include scaled nonzero
tractions, a different active-set stiffness, fixed AMR, and contact-driven AMR.
The 2D and 3D C++ unit suites pass, including new FCC hydrostatic/Schmid-factor
checks (the hydrostatic check is 3D only). The FCC slip-vector constructor,
shared slip-system temporary, and zero-model activation-time initialization
were corrected while preparing this example.

AMR solver repair
-----------------

The conservative residual evaluates a first coarse/fine ghost row. Its incoming
face tangent must also be initialized, and that row must be included in the
multigrid smoother and diagonal. Leaving out these contributions caused the
linearization mismatch and poor AMR convergence. Boundary equations continue
to be evaluated by the elastic BC operator, with physical-domain derivative
stencils; they are not replaced by ghost-cell boundary prescriptions.

The fixed-footprint 2D linear reproducer takes 32 multigrid cycles and a
zero second Newton correction. The legacy fixed-patch 3D mesh's first loaded
solve took 76 cycles (previously 904), with Newton corrections 5e-4, 5.43e-7,
2.88e-11 um.
The two existing soft-void regression cases also pass.

``output_amr_elastic_cycle_coarse`` completed the full load/release schedule
with an 8 x 8 x 4 base grid, two nested levels (32 x 32 x 16 finest-equivalent),
and dt=0.025 s. Peak reaction was 1903.13 uN and final contact count was zero.
This is a machinery check, not a resolved material prediction. A separate
single 50 nm jump on a finer mesh failed nonlinear convergence; do not use
that loading shortcut. ``output_amr_partial`` is an intentionally stopped
plastic run, not a completed load/unload result.

``output_amr_cp_cycle_coarse`` completed all 51 crystal-plasticity solves on
the same coarse nested mesh with dt=0.005 s. Peak reaction was 1924.33 uN,
final contact count was zero, final reaction was -1.82e-9 uN, and maximum
absolute cell-averaged signed slip was 0.00123058. Maximum reported penetration
was 1.55e-13 um and maximum tensile normal traction was 3.08e-7 MPa.
