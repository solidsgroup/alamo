# RESULT: traction-multiplier

Added `elastic.traction_multiplier` (Set::Scalar, default 1.0) to
`elastic_struct` in `src/Integrator/Flame.H`, parsed in `Flame.cpp` via
`pp_query_default("elastic.traction_multiplier", ..., 1.0)`.

Applied at `Flame.cpp:485`:
```cpp
const Set::Scalar traction = (elastic.traction_from_chamber ? chamber.pressure : elastic.traction) * elastic.traction_multiplier;
```

This multiplies the traction actually used regardless of source (constant
`elastic.traction` or dynamic `chamber.pressure`), unlike `elastic.traction`
which is dead when `traction_from_chamber=1`.

## Evidence
- `benchmark/status.sh`: device-lint PASS, golden-compare PASS, a100-sanitizer
  PASS. Golden compare confirms default multiplier=1.0 is bit-identical to
  pre-change behavior (no-op multiply).

## Deviations from plan
None.

## Next step (not part of this task)
User will set `elastic.traction_multiplier = 1.16` in the 5 confirm decks
(anchor, star, centre_bore, rod_and_tube, multifin) for the next batch —
previous t116/t126 batches were no-ops since all decks have
`traction_from_chamber=1`.
