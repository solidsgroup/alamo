# Free circular casing: t=1 s validation

## Outcome

- Added an optional static nodal `casing_support` field to Flame.  With no IC,
  the support is exactly one and the legacy material blend is unchanged.
- The ideal rod-and-tube deck now uses a diffuse circular support of radius
  `0.0877 m`, leaving the square-domain corners as the existing soft void
  material.
- The external rectangle uses zero traction; `xhi` and `ylo` retain the two
  quarter-domain symmetry constraints.
- Face-consistent plot stress now reconstructs across those symmetry planes by
  even material / odd-normal-displacement parity.  Phase-field interfaces are
  not treated as geometric boundaries.
- The chamberutils analytical solver supports `sigma_r(a4)=0`, and the new
  `free_circular_rod_and_tube.json` supplies the actual circular radius rather
  than the old square-domain area-equivalent radius.

## Solver note

The `26 GPa -> 0.5 MPa` casing/exterior contrast made the rediscretized
intra-level coarse hierarchy diverge.  The stable deck recipe keeps all three
AMR levels but sets `elastic.max_coarsening_level=0` and uses the BiCGStab
bottom solve.  The t=1 s static solve converged in two Newton iterations.

## Verification

- `make -j4 bin/test-2d-g++ bin/alamo-2d-g++`: pass.
- `./bin/test-2d-g++`: pass, including conservative and symmetry-parity nodal
  stress reconstruction.
- Eight-rank geometry/elastic smoke runs: pass with the stable solver recipe.
- Eight-rank t=1 s run: pass.  For this validation run only,
  `elastic.interval=4999` solved the uncoupled static elasticity immediately
  before the 1 s plot.
- Corner modulus: `0.5 MPa`; supported casing approaches `26 GPa`.
- Finest-level symmetry-plane phase-band adjacent-node stress jumps are about
  `8-9 kPa`; symmetry-plane shear is a few kPa, against MPa-scale normal stress.
- Correct-radius combined analytical comparison:
  - radial-stress RMS difference: `0.195 MPa`;
  - mean absolute polar shear: `0.0207 MPa`;
  - rod traction transmission: `92.55%`;
  - casing radial stress approaches zero at `a4`;
  - the diffuse Alamo casing hoop curve approaches the sharp-interface
    analytical curve through the support transition.

## Artifacts

- Alamo output: `output_ideal_rod_and_tube_free_circular_t1/05000node`
- Analytical comparison: `results/stress_validate_t1/compare_alamo.png`
- Direct spatial stress/material plot: `results/spatial_stress_t1.png`
- Run log: `run_t1_np8.log`
- Comparison log: `stress_validate_t1.log`

No full 6.5 s run was started, and no commit was created.
