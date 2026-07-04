# MLMG high-contrast evidence bundle

This directory holds concise notes for the 2026-07-02 to 2026-07-04
high-contrast elastic MLMG campaign. Raw logs, plots, and probe scripts are
local evidence and are ignored by this directory's `.gitignore`.

Load these first:

- `../MLMG_HIGH_CONTRAST_FINDINGS.md`: short index and current decision.
- `SETTINGS.md`: current best measured input recipe.
- `RESYNC_COEFFS.md`: stale-hierarchy source fix.
- `NR_PSI_CONVERGENCE.md`: psi-weighted Newton convergence gate.
- `NEGATIVE_RESULTS.md`: refuted or removed experiment paths.

Raw outputs are intentionally not summarized here in full. The most important
artifact names are captured in the focused notes above.

## Key artifact directories

| path | contents |
|---|---|
| `validation_20260703/` | Early 3D CPU/GPU validation and diagonal-inflation refutation. |
| `validation_20260704/` | Resync, psi-update convergence, rod-and-tube seam, and lower-floor robustness screens. |

## Current best evidence

| artifact | headline |
|---|---|
| `validation_20260704/codex_pf0_resync1_20260704_123913.out` | `resync_coeffs=1` makes the floor-0 anchor Newton linear solves converge. |
| `validation_20260704/codex_pf0_resync0_20260704_123929.out` | Control with `resync_coeffs=0` reproduces stale-hierarchy divergence. |
| `validation_20260704/codex_perf_new_pf010_resync_20260704_124819.out` | `psi_floor=0.01`, resync, line search, and zeroed displacement complete the first elastic window in about 8.0 s for a fixed 3-Newton screen. |
| `validation_20260704/codex_nrpsi_pf010_20260704_134719.out` | `nr_convergence=psi_update` stops at 13 Newton solves, about 27.4 s. |
| `validation_20260704/codex_nrdiag_pf010_update_nrtol_20260704_134300.out` | Raw-update reference takes 178 Newton solves, about 215 s. |
