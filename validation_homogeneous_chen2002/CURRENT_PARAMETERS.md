# Selected binder heat capacity and parameter set

Retain **cp=2130 J/(kg K)** with its previously calibrated binder A and E/R.
The user selected this model after comparing the cp=2418.29 alternative.
`parameters_current.json` is a byte-for-byte copy of
`ap_endpoint_calibration/parameters_frozen.json`. All 18 associated standard
simulations were checked again for successful execution, complete duration,
matching input/executable hashes and matching parameters. Reusing these
completed results requires no new calibration or simulations.

## Physical basis

Hanson-Parr and Parr report an HTPB heat capacity of approximately
0.51 cal/(g K) at 21°C, equivalent to **2133.84 J/(kg K)**. Our rounded 2130
value therefore has a direct measurement basis near the cold-solid temperature.
The value is visible on printed page 26 of the indexed original-paper copy.
The complete specimen/cure description was not verified, so this does not
identify the heat capacity of a particular user formulation or its hot state.
[Original paper (1999)](https://doi.org/10.1080/07370659908216094),
[indexed original PDF](https://electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20Energetic%20Materials/1999%20Volume%2017/1/1-48.pdf).

The **2418.29 J/(kg K)** value also appears as an HTPB thermal property in
Table 2 of *Diffuse interface method for solid composite propellant ignition
and regression*. This establishes its use in a combustion model; it does
not establish that it is the more accurate value for this binder.
[Model paper, Table 2](https://engrxiv.org/preprint/download/3085/5643/4453).

The two constants differ by 13.53%. Neither is established as a universal
heat capacity across the heated solid layer. For a physical material model,
the relevant sensible enthalpy is the integral of formulation-specific cp(T)
from T0 to Ts; constant cp approximates that integral. The available evidence
supports 2130 as a plausible modeling choice, without proving that either
constant is the uniquely correct value at the simulated surface temperatures.
The smaller Chen error is the user's selection criterion among these model
choices, not independent evidence that the material property is exact.

## Retained values and results

| Parameter | Binder | AP |
|---|---:|---:|
| Density [kg/m³] | 920 | 1950 |
| cp [J/(kg K)] | 2130 | 1297.9 |
| Conductivity [W/(m K)] | 0.213 | 0.4186 |
| Q [cal/g] | −66 | −100 |
| A [m/s] | 24.6083329584 | 1179.58476192 |
| E/R [K] | 5568.85064240 | 10739.1532211 |
| E [kJ/mol] | 46.30200049 | 89.29028801 |
| E [kcal/mol] | 11.06644371 | 21.34089102 |

T0=300 K. Volume-additive density, mass-weighted cp and Q, volume-weighted
ln(A) and E/R, the Chen conductivity rule, and frozen gas chemistry are retained.
The cp=2418.29 kinetic pair is not combined with cp=2130.

| Study | Mean absolute error [%] | Maximum absolute error [%] |
|---|---:|---:|
| Selected cp=2130, matching A/E | 0.1314 | 0.2790 |
| Archived cp=2418.29, refitted A/E | 0.6164 | 1.2143 |

These errors compare all 18 LowMach points with Chen's Figure 4 solid model
curves. See the [selected plot](ap_endpoint_calibration/figure4_comparison.png)
and [report](ap_endpoint_calibration/REPORT.md). The
[cp=2418.29 study](binder_cp2418_calibration/REPORT.md) and its raw outputs
remain intact as the comparison study.

For future runs, explicitly select the retained configuration with
`--parameters validation_homogeneous_chen2002/parameters_current.json --dt-scale .2`
when invoking `run_study.py`, using a fresh run tag. The generator's older
default and all historical study inputs remain intact. `current_selection.json`
records the selected source, parameter hash, figure/report paths and errors.
