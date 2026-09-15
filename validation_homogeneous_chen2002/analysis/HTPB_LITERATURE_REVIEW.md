# AP-free HTPB: literature audit and calibration disposition

Reviewed 2026-09-09. Literature search completed; no new calibration accepted,
no simulations restarted, and no material input changed. The last two
Arrhenius runs have completion receipts, but their parameter set remains
rejected. Earlier outputs are preserved.

## Material identity and physical properties

Neat liquid HTPB prepolymer is not the same material as an AP-free cured
elastomer. Cured binder contains reacted curative and can contain plasticizer
without containing AP. A numerical-model constant is not automatically a
measurement of either material.

| Property | Published evidence, converted to SI | Limitation |
|---|---|---|
| Density | Approximately 900 kg/m³ at 23°C [1] | Liquid R45 HTLO resin, not a cured-binder measurement |
| Heat capacity | Approximately 2130 J/(kg K) at 21°C [2] | Original paper's HTPB sample; cure formulation not verified in accessible excerpt |
| Conductivity | Approximately 0.201 W/(m K) at 21°C [2] | Same sample qualification as cp |
| Thermal diffusivity | 1.08 × 10⁻⁷ m²/s at 21°C [2] | Not independently adjustable once density, cp and k are specified |
| Heat-capacity temperature dependence | Approximately 1880 to 1710 J/(kg K), from 340 to 410 K [3] | Separately measured HTPB constituent; not the paper's cured composite |
| Glass transition | −75°C [1] | R45 HTLO specification, not a melting or pyrolysis temperature |

[1] Cray Valley, *Poly bd R45 HTLO Technical Data Sheet*, revision April
2020, [manufacturer PDF](https://crayvalley.com/download/6/low-vinyl-homopolymers-tech-datasheets/456/poly-bd-r-45).
Specific gravity is 0.90 at 23°C. The sheet explicitly identifies a liquid
resin; its density does not establish that 920 kg/m³ is wrong for a cured
binder.

[2] D. M. Hanson-Parr and T. P. Parr, *Thermal properties measurements of
solid rocket propellant oxidizers and binder materials as a function of
temperature*, Journal of Energetic Materials 17 (1999), 1–48,
[DOI](https://doi.org/10.1080/07370659908216094).
The indexed original-paper excerpt on printed page 26 reports HTPB
cp = 0.51 cal/(g K), k = 0.00048 cal/(cm s K), and diffusivity
0.00108 cm²/s at 21°C. Conversions use 1 cal = 4.184 J.
[Indexed paper copy](https://electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20Energetic%20Materials/1999%20Volume%2017/1/1-48.pdf).
Direct download failed; the full specimen description and temperature
curves were not inspected. These rounded values are screening evidence,
not a fully verified hot-state property set.

[3] *Research on the Specific Heat Capacity of PBX Formulations Based on
RDX*, Journal of Aerospace Technology and Management 8 (2016),
[full article and Table 1](https://www.scielo.br/j/jatm/a/dsKjKT3cZWxC4jSkt6bDs8S/?lang=en),
DOI 10.5028/jatm.v8i3.655. The HTPB-only constituent data must not be confused
with the composite results. The tabulated temperature polynomial gives the
endpoints above. It should not be extrapolated to pyrolysis temperatures.
Disagreement among HTPB measurements is not grounds to choose whichever
value best matches a regression curve.

[4] Veals et al., *Property Estimates for Hydroxyl-Terminated Polybutadiene
(HTPB) Type R45M Derived from Atomistic Molecular Dynamics Simulations*,
ARL-TR-9714, June 2023,
[original report reproduced online](https://www.scribd.com/document/918681583/Ad-1204916).
Section 2 distinguishes liquid-resin density measurements from thermal
measurements on IPDI-cured R45M over 263–413 K. Figure 7 compares measured
and predicted cp, increasing with temperature. Reliable numerical figure
ordinates were not available, so no curve was fitted. The DTIC PDF could
not be retrieved. Its modeled oligomer vaporization enthalpy is not the
heat of chemical pyrolysis.

## Decomposition energetics

[5] Y.-C. Lu and K. K. Kuo, *Thermal decomposition study of
hydroxyl-terminated polybutadiene (HTPB) solid fuel*, Thermochimica Acta
275 (1996), 181–191,
[original-paper abstract](https://www.sciencedirect.com/science/article/pii/0040603195027262).
Their resin/curative/cured-polymer comparison finds that curing affects
decomposition energetics and heating rate changes the thermal response.
It does not establish a universal decomposition heat for all AP-free HTPB.

[6] E. Hagen et al., *Effective enthalpy of pyrolysis of HTPB under
diffusion flame conditions*, Chemical Engineering Journal, article 180235,
online 3 August 2026,
[publisher abstract and introduction](https://www.sciencedirect.com/science/article/abs/pii/S1385894726076965),
DOI 10.1016/j.cej.2026.180235. Experiments support condition-dependent
effective pyrolysis enthalpy as product distributions change. The
introduction discusses contemporary estimates around 1–2 MJ/kg; these are
not universal bounds or measurements of our unspecified sample. Their
material is a cured, plasticized formulation. The accessible text does not
justify choosing a replacement constant Q here; the full paper was not
retrieved.

[7] H. Arisawa and T. B. Brill, *Flash pyrolysis of hydroxyl-terminated
polybutadiene (HTPB) I: Analysis and implications of the gaseous products*,
Combustion and Flame 106 (1996), 131–143,
[original-paper abstract](https://www.sciencedirect.com/science/article/pii/0010218096002532).
Measurements cover approximately 723–882 K and identify multiple products,
including oligomers. This is an experimental range, not a universal onset
temperature; it does not support a below-room-temperature hot-pyrolysis fit.

Heat of formation, physical vaporization heat, chemical decomposition heat,
and total gasification heat including sensible heating are different
quantities. Matching units do not make them interchangeable. Cp must not
silently absorb a reaction heat already included elsewhere in the energy
equation.

## Why an Arrhenius-only calibration is still unresolved

Chen Figure 4 is a model/DNS comparison, not a pure-HTPB experimental
dataset. The supplied article refers constituent properties to Table 2 of
its reference [3]. See `~/Downloads/homogeneous_model_paper.pdf`, printed
page 2926, equations 12–14. The prescribed flux and target rate determine
the required surface energy balance before Arrhenius parameters enter.

The existing [feasibility assessment](arrhenius_htpb/FEASIBILITY.md)
identifies a conflict with the retained Gross endothermic Q. Rechecking the
energy requirement with the approximate literature resin density does not
resolve it: at Chen's lowest-flux target the available energy is about
1.16 MJ/kg, below the retained decomposition cost of 1.2552 MJ/kg before
any positive sensible heating. At 920 kg/m³ the available energy is about
1.13 MJ/kg. Positive cp cannot change the sign of that deficit. Conductivity
changes the thermal profile but cannot remove this steady prescribed-flux
energy-balance conflict. No new A/E values were fitted.

This is conditional on the retained Q and flux definition, not a proof that
every possible HTPB model conflicts with Chen. The literature search has not
established a formulation- and condition-matched replacement Q. Choosing
one solely to recover Chen's curve would repeat thermal-parameter fitting.

## Disposition

The room-temperature evidence suggests revisiting k = 0.13 W/(m K), but
does not establish one replacement valid throughout decomposition. Existing
cp = 2418.29 J/(kg K) and density = 920 kg/m³ need formulation/temperature
provenance; the search does not establish that either already includes AP.

These literature values are not installed as a runnable parameter set.
No source, production input, AP parameter, Q, kinetic parameter or numerical
cutoff is changed. A physically supported calibration needs a matched
material/thermal basis and compatible observations, not another
unconstrained fit to the same model curve.
