# TASK: free-circular-casing
# Folder: docs/agent_plans/20260722-free-circular-casing/

---

## Header

| Field        | Value |
|--------------|-------|
| Risk tier    | 2 (new optional material mask and changed rod/tube mechanics) |
| Model        | Codex root session |
| Verification | partial-oracle |
| Est. scope   | Flame field/blend, mechanics output boundary, ideal deck, chamberutils reference |
| Parallel-safe| no: shared Flame/Mechanics source and live dirty worktree |

## Operating rules

1. Preserve existing behavior unless the new casing-support IC is configured.
2. Do not reuse burn field `eta` for exterior space or apply chamber pressure
   at the casing/air interface.
3. Keep right/bottom quarter-symmetry constraints; make top/left zero traction
   only in the ideal rod-and-tube deck.
4. Do not commit until the user accepts the new t=1 s plots.

## Objective

Represent the rod-and-tube casing as a circular annulus surrounded by soft
void inside the rectangular AMR domain, with a free exterior rather than a
clamped exterior. Remove the output-stress scheme discontinuity on the right
and bottom symmetry planes using symmetry parity. Update the analytical
validation to use a traction-free casing exterior.

## Design

- Add an optional static nodal casing-support field, defaulting to one so all
  existing Flame decks retain their current material blend.
- With support `c`, use the partition
  `prop=c*phi*eta`, `gas/exterior=c*phi*(1-eta)+(1-c)`,
  `casing=c*(1-phi)`.
- Configure `c` in the ideal rod-and-tube deck as a diffuse circle of radius
  `a4=0.0877 m`, centered at `(0.0877, 0.0877)`, with a resolved transition.
- Tag the support transition for refinement; do not couple it to pressure RHS.
- Use zero traction on `xlo` and `yhi`; retain `xhi=disp trac` and
  `ylo=trac disp` symmetry.
- Extend displacement across only `xhi` and `ylo` with the appropriate vector
  parity before applying the existing centered face reconstruction. Do not
  treat `phi`, `eta`, or the casing-support interface as geometric boundaries.
- Change the chamberutils analytical outer condition from `u(a4)=0` to
  `sigma_r(a4)=0` for this validation case.

## Oracle

1. Focused 2-D build and unit test.
2. Initialization/short run confirms void-valued corner modulus, circular
   casing support, finite stress, and unchanged pressure-loading surface.
3. Eight-rank ideal rod-and-tube run to t=1 s.
4. Generate stress and stress_validate plots; compare symmetry-plane and
   interior phi-interface stress without running the full 6.5 s case.

## Checkpoints

- [x] Human confirms the optional support-mask/parity/free-outer design before source edits
- [ ] Diff and t=1 s plots reviewed before completion or commit

## Closeout

- [x] Results and limitations recorded in `results/RESULT.md`
- [ ] Accepted outcome appended to `docs/llm/SESSION_LOG.tsv`
- [ ] No commit without explicit user approval
