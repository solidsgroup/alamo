# Standalone homogeneous LowMach comparison with Chen et al. (2002)

Completed study: nine baseline runs and two key-case refinements.

The completed planar LowMach runs reproduce the **chosen homogeneous
closure** to within 3.23% on the baseline setup.
Baseline speed errors range from -3.23% to -2.22%, with
late-window speed drift at most 0.19%.
The most refined key-case result is 1.03421 cm/s,
-0.89% from that closure.
This is a numerical comparison using a Figure-4-calibrated material
realization; it is **not independent validation of uniquely recovered material
properties or a reproduction of the article's heterogeneous DNS**.

The mixed-composition baseline predictions differ from the article's Fig. 4
solid curves by -12.64% to -12.33%. Most of
that discrepancy already exists in the selected analytic closure before
running LowMach. It is conditional on the assumed constituent densities,
the arithmetic mass-weighted density rule, and the digitized endpoint fits.
It is not evidence that the authors used an incorrect rule: the supplied
article does not uniquely identify the separate material parameters, and
graphical extraction has uncertainty.

## Scope and reproducibility

The study uses the existing `bin/lowmach-2d-clang++` and input files only.
No production source changes are introduced for this study. Case inputs,
material assumptions, binary/input hashes, execution receipts, and raw AMReX
plots are retained under `../runs/`. Failed setup pilots remain as diagnostic
history and are excluded from the completed-run tables.

The nine baseline combinations are q=200, 500, 1000 cal/(cm² s) and AP
volume fraction t=0, 0.5, 1. The strip is periodic transversely, with 32
transverse cells and spatially uniform initial conditions in that direction.
There is no resolved particle packing or gas chemistry. The homogeneous
`phase_change` implementation evolves the interface and temperature.

An input-defined fixed volumetric heat-source slab in the gas supplies the
requested net heat input. The upper temperature boundary is insulated and
the lower solid boundary is held at 300 K. The source does not prescribe
surface motion. The deliberately artificial conducting gas has k=100 W/(m K),
cp=1000 J/(kg K), R=319.787 J/(kg K), and P=100 MPa; temperature advection
is disabled. Projection is retained for rigid phase-change coupling. This
gas layer is a numerical surrogate for the paper's imposed surface flux,
not a physical combustion-gas model.

The reference thermal length is delta=k_s/(rho_s*c_s*r_ref). Baseline
interface width parameter ell=delta/4 and dy=ell/8. Each case runs for
six delta/r_ref, with roughly forty output intervals. The analytic
temperature profile initializes the run and sizes the domain; it does not
force the evolving rate. Surface positions are independently interpolated
at eta=0.5 in every transverse column, then fitted over the last 40% of
physical time. An independent solid-volume-loss rate checks that fit.

## Article-only reference and parameter identifiability

The scientific reference is only the supplied Chen, Buckmaster, Jackson,
and Massa article, *Homogenization issues and the combustion of heterogeneous
solid propellants*, Proceedings of the Combustion Institute 29 (2002),
2923–2929. Fig. 4 is printed p. 2927 / PDF page 5.

Its actual PDF vector paths were extracted, including solid Eq. 15 curves,
dashed curves, 12 Eq. 19 squares, 24 2D DNS circles, and 24 3D DNS asterisks.
`figure4_provenance.json` records the source hash and coordinate calibration.
The reconstructed reference is `figure4_extracted.pdf`; CSV decimals preserve
vector coordinates, not scientific precision. A full printed stroke spans
approximately 0.005 cm/s and 0.0011 in packing fraction. This is graphical
resolution, not a confidence interval; source DNS uncertainties are unavailable.

The article specifies A_B=1036 cm/s, E_B/R=7500 K, A_AP=94800 cm/s,
E_AP/R=11000 K. At pure-component endpoints, Eqs. 12 and 14 imply

    Ts = (E/R)/ln(A/r),
    q/r = alpha*(Ts-beta),  alpha=rho*c,  beta=T0+Q/c.

Fitting the three endpoints gives alpha_B≈1.155 MJ/(m³ K), beta_B≈143 K,
alpha_AP≈2.435 MJ/(m³ K), beta_AP≈−38 K. Separate rho, c, Q, and T0 are
not identified. The numerical realization assumes rho_B=920 kg/m³,
rho_AP=1950 kg/m³, T0=300 K, then derives c=alpha/rho and Q=c*(beta−T0).
The densities are assumed repository values, not article-recovered constants.
Conductivities are estimated from the article's approximate thermal lengths
47 µm binder and 33 µm AP at r=0.5 cm/s. These transport estimates are not
independently validated by steady burning-rate agreement.

For density ratio d=rho_AP/rho_B, arithmetic mass-weighted density gives

    (rho*c)_blend = [(1−t)+d²*t]/[(1−t)+d*t]²
                   * [(1−t)*alpha_B+t*alpha_AP].

Conditional on our endpoint fits, the Fig. 4 Eq. 19 squares are consistent
to graphical resolution with a multiplier near one on the last bracket.
The assumed density ratio gives a larger multiplier and lower mixed rates.
This is a density-convention/parameter sensitivity, not proof of an error
in the publication. The article's prose, unidentified parameters, and
graphical sensitivity prevent a unique resolution. Full derivations and
endpoint-perturbation sensitivity are retained in `REFERENCE_NOTES.md`,
`article_inferred_combinations.json`, and `figure4_blend_capacity_audit.csv`.

Pure-endpoint agreement is calibration to the same figure. Eq. 15 combines
pure-component speeds at the same q but separate surface temperatures;
Eq. 19 uses geometric Arrhenius parameters at the blend temperature from
Eq. 14. These definitions must be distinguished even where their published
curves nearly coincide.

## Nine-case baseline comparison

Rates are in cm/s. The chosen closure is Eqs. 14+19 with the disclosed
material realization. The Eq. 15 column uses the same fitted pure endpoints;
the last error uses the independently extracted Fig. 4 solid curve directly.

| q [cal/(cm² s)] | t | LowMach | Chosen closure | Error [%] | Eq. 15 | Error vs Fig. 4 [%] |
|---:|---:|---:|---:|---:|---:|---:|
| 200 | 0.0 | 0.78089 | 0.80197 | -2.63 | 0.80197 | -2.62 |
| 200 | 0.5 | 0.44513 | 0.45599 | -2.38 | 0.50907 | -12.64 |
| 200 | 1.0 | 0.36440 | 0.37289 | -2.28 | 0.37289 | -2.25 |
| 500 | 0.0 | 1.70534 | 1.75525 | -2.84 | 1.75525 | -2.85 |
| 500 | 0.5 | 1.01884 | 1.04352 | -2.37 | 1.16436 | -12.53 |
| 500 | 1.0 | 0.85133 | 0.87111 | -2.27 | 0.87111 | -2.32 |
| 1000 | 0.0 | 3.04723 | 3.14888 | -3.23 | 3.14888 | -3.22 |
| 1000 | 0.5 | 1.89790 | 1.94270 | -2.31 | 2.16516 | -12.33 |
| 1000 | 1.0 | 1.61308 | 1.64976 | -2.22 | 1.64976 | -2.20 |

`comparison.pdf` distinguishes baseline crosses, smaller-time-step triangles,
and thinner-interface diamonds. Its lower panel reports error against the
chosen closure, not against the paper's Eq. 19 squares. The paper DNS markers
are resolved heterogeneous calculations, not measurements from this study.

Baseline eta=0.5 surface temperatures exceed their chosen-closure values
by 2.26–18.53 K. At finite interface width,
temperature varies across the phase-change region; applying the Arrhenius
law to one interpolated surface temperature is not identical to integrating
the distributed phase-field kinetics. The discrepancy requires interface
and time-step assessment rather than an unqualified surface-law claim.

## Key-case numerical refinement: q=500, t=0.5

| Setup | ell/delta | dy [µm] | Step-ceiling factor | Rate [cm/s] | Closure error [%] | Late speed drift [%] |
|---|---:|---:|---:|---:|---:|---:|
| Baseline | 0.250 | 0.4876 | 1.00 | 1.01884 | -2.37 | -0.011 |
| Smaller time step | 0.250 | 0.4876 | 0.10 | 1.03066 | -1.23 | -0.024 |
| Thinner interface + smaller step | 0.125 | 0.2438 | 0.05 | 1.03421 | -0.89 | -0.035 |

The smaller-step case reduces the allowed maximum step; dynamic constraints
also limit the baseline, so the actual step ratio is not uniformly ten.
The thinner-interface case halves both ell and dy and also halves the
absolute time-step ceiling relative to the smaller-step case. It is a
coupled interface/grid/time refinement, not an observed spatial convergence
order or a complete uncertainty estimate. This selected refinement does not
establish convergence for all nine cases.

## Incident-flux and thermal-equation checks

The discrete source integral equals its prescribed q. Conductive flux is
reconstructed with harmonic face conductivity from stored temperature and
`thermal_conductivity_coeff`. The incident sample uses eta=10⁻⁶, outside
transient cooling ripples in the diffuse tail. Fluxes inside that tail can
be large because high conductivity magnifies small split-step temperature
ripples; they are not interpreted as the imposed incident flux.

Gas and solid thermal storage are reconstructed using the actual smoothed
solid fraction in the thermal-capacity model. Between outputs, the gas
temperature primitive is logarithmic because C_g is proportional to 1/T;
solid storage uses rho*c*H(eta). The residual compares integrated source,
latent heat from measured volume loss, bottom loss, and thermal storage.
This snapshot estimate has output quadrature and splitting error; it is not
an exact solver-native conservation audit. In particular, output spacing
corresponds to approximately 0.15 delta of front motion, larger than the
thin case's ell=0.125 delta. Its moving-interface storage can therefore be
under-resolved in time even when the solver's own timestep is small.

A second late-time check replaces only the solid storage with
r_measured*integral(C_s*dT/dy), averaged over the late profiles. This assumes
the condensed thermal profile translates steadily, supported by the measured
late-speed/temperature stability. Gas storage is still measured separately;
the fixed gas heating layer is not assumed to translate. Both residuals
are reported so that output-quadrature error is not mistaken for solver
error or a failure of interface refinement.

| Key-case setup | Incident flux error [% q] | Gas storage [% q] | Snapshot residual [% q] | Steady-solid residual [% q] |
|---|---:|---:|---:|---:|
| Baseline | -1.136 | +0.581 | +1.814 | +1.353 |
| Smaller time step | -0.537 | +0.535 | +0.839 | +0.320 |
| Thinner interface + smaller step | -0.487 | +0.637 | +1.287 | +0.094 |

Across the baseline cases the incident-flux error is -1.26% to
-1.13%, gas storage is 0.53–0.88% of q,
and the reconstructed residual is 1.58–1.94% of q.
These numbers bound the demonstrated fidelity of the conducting-layer
surrogate; independence from gas conductivity or pressure was not established.
Near q=500, t=0.5, the chosen closure has d(log r)/d(log q)≈0.90, so a
1% incident-flux bias alone can contribute roughly 0.9% rate bias. It must
not be hidden inside a claim of fully converged agreement.

## Independent planar kinetic calibration

The isothermal calibration has exact speed 250 mm/s with activation set to
zero. Independent raw-output extraction gives:

| T [K] | dy [µm] | Measured rate [mm/s] | Error [%] |
|---:|---:|---:|---:|
| 300 | 5.0 | 240.687 | −3.73 |
| 300 | 2.5 | 246.367 | −1.45 |
| 600 | 5.0 | 241.848 | −3.26 |
| 600 | 2.5 | 247.423 | −1.03 |

Fine-grid late drift is below 0.004%, and integrated-volume rates agree
with interface rates within 0.01%. Two mesh levels show improvement but do
not determine an observed order; Richardson estimates in the calibration
report assume second order. This verifies the kinetic-to-speed mapping,
not the coupled Figure 4 closure. Raw independent results are retained in
`calibration_independent_summary.csv`.

## Reproducing the analysis

From the repository root, using the configured Python environment:

```bash
MPLCONFIGDIR=/tmp/chen-mpl /home/esandall/Software/anaconda3/bin/python \
  validation_homogeneous_chen2002/analysis/digitize_figure4.py \
  /home/esandall/Downloads/homogeneous_model_paper.pdf
/home/esandall/Software/anaconda3/bin/python \
  validation_homogeneous_chen2002/analysis/infer_article_combinations.py
MPLCONFIGDIR=/tmp/chen-mpl /home/esandall/Software/anaconda3/bin/python \
  validation_homogeneous_chen2002/analysis/compare_study.py
/home/esandall/Software/anaconda3/bin/python \
  validation_homogeneous_chen2002/analysis/write_report.py
```

`simulation_summary.csv` and `.json` contain all case-level measurements;
`timeseries_*.csv` retain raw-output-derived histories;
`comparison.pdf` and `numerical_diagnostics.pdf` are vector figures with
matching PNG previews. Exact simulation commands and their input/binary
hashes are in each case's `run.json`; retained input decks are the executable
specification of each run. `../run_study.py` provides generation and execution.
