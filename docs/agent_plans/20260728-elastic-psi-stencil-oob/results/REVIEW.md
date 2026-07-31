# Adversarial review — N11 elastic psi stencil OOB

Reviewed 2026-07-31 against fix commit `c465264cf` and call-site adjudication
commit `3bb84ff7d`. The core `Elastic::Diagonal` fix has no finding: growing only
the cell-centered `psi` coefficient to three ghosts matches the two-nodal-ghost
launch plus one backing cell, preserves boundary interpolation, and was not a
sanitizer-only suppression.

The follow-up adjudication has two findings.

## High — PhaseFieldMicrostructure dynamic eta halo remains undersized

`mechanics.type=dynamic` is accepted in
`src/Integrator/PhaseFieldMicrostructure.H`, but only static mechanics raises
eta storage to three ghosts. `UpdateModel` launches over a two-ghost nodal model
box and calls the default lower-corner `CellToNodeAverage`; dynamic eta therefore
has one ghost, or two with anisotropy, where the launch contract requires three.
`UpdateEigenstrain` carries the same two-ghost launch for every non-disabled
mechanics type.

This is a supported CPU configuration even though the GPU dynamic path aborts
earlier. Calling it dormant and declining a follow-on in `RESULT.md` §5 is not
adequate. A separate task must either require three eta ghosts whenever
mechanics is enabled or constrain the launch paths to the available source
extent.

Evidence: `PhaseFieldMicrostructure.H:96,114`;
`PhaseFieldMicrostructure.cpp:251,260,270,336,344,353`.

## Medium — `amr.print_ghost_nodes` exceeds the admitted source contract

`amr.print_ghost_nodes` directly grows the nodal plot BoxArray, while the
cell-source admission check requires only one ghost.
`AverageCellcenterToNode` then launches over that enlarged nodal box. At
`print_ghost_nodes=1`, a one-ghost cell field can be read at its second lower
ghost by the default cell-to-node average.

The option is user-facing and therefore not non-live merely because its default
is zero. A separate task must validate/reject this setting for cell-to-node
plotting or require and fill `print_ghost_nodes + 1` source ghosts.

Evidence: `Integrator.cpp:1032,1062`; `Util.cpp:421`;
`Numeric/Stencil.H:1570`.

## Verdict

**REVISE.** The core N11 source fix may remain. Do not mark the task done until
the two follow-on task folders are opened and a human confirms that disposition.
