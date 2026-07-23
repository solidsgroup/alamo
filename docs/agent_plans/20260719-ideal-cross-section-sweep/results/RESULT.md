# RESULT — ideal-cross-section-sweep (2026-07-19)

## Outcome: 5/5 PASS

All five geometries ran to t=6.4998s (stop_time 6.5s) at the IDEAL recipe:
aluminum casing 70/26 GPa, void 0.5/0.5 MPa, psi OFF, n_cell 128 full /
64 quarter (rod_and_tube), dt 2.0e-4, np=8 sequential, bin/alamo-2d-clang++.
Modulus contrast 1.4e5x (70 GPa vs 0.5 MPa) held on every geometry with
ZERO solver failures / NaN / aborts in all five run logs.

## Per-geometry

| geometry     | exit | wall     | output size | final t |
|--------------|------|----------|-------------|---------|
| rod_and_tube | 0    | 1631s  (27m)  | 5.7G  | 6.4998 |
| multifin     | 0    | 4518s  (75m)  | 22G   | 6.4998 |
| anchor       | 0    | 10405s (2h53m)| 18G   | 6.4998 |
| star         | 0    | 5448s  (91m)  | 17G   | 6.4998 |
| centre_bore  | 0    | 4096s  (68m)  | 12G   | 6.4998 |

Total wall: 7h15m (13:55:07 -> 21:10:35). Total disk: ~75G (130G free after).

## Artifacts

- Decks: input_ideal_{rod_and_tube,multifin,anchor,star,centre_bore}
  (sed-derived from input_sweep_*; diffs verified = only casing/void/n_cell/
  plot_file lines)
- Driver: run_ideal_sweep.sh; log ideal_sweep_driver.log
- Outputs: output_ideal_<g>/ each with thermo.dat + thermo_plots.png
  (chamberutils thermo.py --show-regression, exit 0 x5)
- Run logs: run_ideal_<g>.log

## Deviations from plan

None. Oracle met (exit=0 x5, t >= 6.49 x5).

## Notes

- multifin (historical cliff geometry) clean at full contrast — psi-OFF +
  0.5/0.5 void recipe needs no per-geometry tuning at n_cell 128.
- Known Mechanics.H:424 parallel trac/disp thermo.dat warning present as
  always; use boxlib output for stress-derived quantities.
- Next: cross-geometry comparison figures via chamberutils.
