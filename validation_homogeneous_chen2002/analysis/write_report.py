#!/usr/bin/env python3
"""Write the reviewable study report from independently measured case results.

Updates STUDY_REPORT.md throughout the study, explicitly marking the missing
interface refinement until all nine baseline cases and both refinements exist.
"""
import csv
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def main():
    rows=json.loads((HERE/"simulation_summary.json").read_text())
    baseline=sorted([r for r in rows if r["width_ratio"]==.25 and r.get("dt_scale",1)==1 and r["name"].endswith("_nx32")],key=lambda r:(r["q_cal_cm2_s"],r["ap_volume_fraction"]))
    key=sorted([r for r in rows if r["q_cal_cm2_s"]==500 and r["ap_volume_fraction"]==.5 and (r in baseline or r.get("dt_scale",1)<1)],key=lambda r:(r["width_ratio"]<.25,r.get("dt_scale",1)<1))
    complete=len(baseline)==9 and any(r["width_ratio"]<.25 for r in key) and any(r["width_ratio"]==.25 and r.get("dt_scale",1)<1 for r in key)
    if not baseline:raise SystemExit("No completed baseline cases")
    def limits(items,field,absolute=False):
        vals=[abs(r[field]) if absolute else r[field] for r in items]
        return min(vals),max(vals)
    errors=limits(baseline,"error_eq14_19_percent")
    drift=max(abs(r["late_speed_drift_percent"]) for r in baseline)
    flux=limits(baseline,"incident_flux_error_percent")
    storage=limits(baseline,"gas_storage_percent_input")
    balance=limits(baseline,"thermal_balance_residual_percent_input")
    Ts_error=limits(baseline,"error_surface_temperature_K")
    paper_middle=[r for r in baseline if r["ap_volume_fraction"]==.5]
    middlerange=limits(paper_middle,"error_figure4_solid_percent")
    refined=min(key,key=lambda r:(r["width_ratio"],r["dx_m"],r.get("dt_scale",1)))
    text=f"""# Standalone homogeneous LowMach comparison with Chen et al. (2002)

{'Completed study: nine baseline runs and two key-case refinements.' if complete else f'DRAFT: {len(baseline)}/9 baseline cases and {len(key)-1}/2 refinements are currently included; the final interface assessment remains pending.'}

The completed planar LowMach runs reproduce the **chosen homogeneous
closure** to within {max(abs(x) for x in errors):.2f}% on the baseline setup.
Baseline speed errors range from {errors[0]:+.2f}% to {errors[1]:+.2f}%, with
late-window speed drift at most {drift:.2f}%.
The most refined key-case result is {refined['r_cm_s']:.5f} cm/s,
{refined['error_eq14_19_percent']:+.2f}% from that closure.
This is a numerical comparison using a Figure-4-calibrated material
realization; it is **not independent validation of uniquely recovered material
properties or a reproduction of the article's heterogeneous DNS**.

The mixed-composition baseline predictions differ from the article's Fig. 4
solid curves by {middlerange[0]:+.2f}% to {middlerange[1]:+.2f}%. Most of
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
"""
    for r in baseline:
        text+=f"| {r['q_cal_cm2_s']:.0f} | {r['ap_volume_fraction']:.1f} | {r['r_cm_s']:.5f} | {r['expected_r_eq14_19_cm_s']:.5f} | {r['error_eq14_19_percent']:+.2f} | {r['expected_r_eq15_cm_s']:.5f} | {r['error_figure4_solid_percent']:+.2f} |\n"
    text+=f"""
`comparison.pdf` distinguishes baseline crosses, smaller-time-step triangles,
and thinner-interface diamonds. Its lower panel reports error against the
chosen closure, not against the paper's Eq. 19 squares. The paper DNS markers
are resolved heterogeneous calculations, not measurements from this study.

Baseline eta=0.5 surface temperatures exceed their chosen-closure values
by {Ts_error[0]:.2f}–{Ts_error[1]:.2f} K. At finite interface width,
temperature varies across the phase-change region; applying the Arrhenius
law to one interpolated surface temperature is not identical to integrating
the distributed phase-field kinetics. The discrepancy requires interface
and time-step assessment rather than an unqualified surface-law claim.

## Key-case numerical refinement: q=500, t=0.5

| Setup | ell/delta | dy [µm] | Step-ceiling factor | Rate [cm/s] | Closure error [%] | Late speed drift [%] |
|---|---:|---:|---:|---:|---:|---:|
"""
    for r in key:
        name="Baseline" if r in baseline else "Smaller time step" if r["width_ratio"]==.25 else "Thinner interface + smaller step"
        factor=r.get("dt_scale",1)*r["width_ratio"]/.25
        text+=f"| {name} | {r['width_ratio']:.3f} | {r['dx_m']*1e6:.4f} | {factor:.2f} | {r['r_cm_s']:.5f} | {r['error_eq14_19_percent']:+.2f} | {r['late_speed_drift_percent']:+.3f} |\n"
    text+="""
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
"""
    for r in key:
        name="Baseline" if r in baseline else "Smaller time step" if r["width_ratio"]==.25 else "Thinner interface + smaller step"
        text+=f"| {name} | {r['incident_flux_error_percent']:+.3f} | {r['gas_storage_percent_input']:+.3f} | {r['thermal_balance_residual_percent_input']:+.3f} | {r.get('steady_profile_balance_residual_percent_input',float('nan')):+.3f} |\n"
    text+=f"""
Across the baseline cases the incident-flux error is {flux[0]:+.2f}% to
{flux[1]:+.2f}%, gas storage is {storage[0]:.2f}–{storage[1]:.2f}% of q,
and the reconstructed residual is {balance[0]:.2f}–{balance[1]:.2f}% of q.
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
MPLCONFIGDIR=/tmp/chen-mpl /home/esandall/Software/anaconda3/bin/python \\
  validation_homogeneous_chen2002/analysis/digitize_figure4.py \\
  /home/esandall/Downloads/homogeneous_model_paper.pdf
/home/esandall/Software/anaconda3/bin/python \\
  validation_homogeneous_chen2002/analysis/infer_article_combinations.py
MPLCONFIGDIR=/tmp/chen-mpl /home/esandall/Software/anaconda3/bin/python \\
  validation_homogeneous_chen2002/analysis/compare_study.py
/home/esandall/Software/anaconda3/bin/python \\
  validation_homogeneous_chen2002/analysis/write_report.py
```

`simulation_summary.csv` and `.json` contain all case-level measurements;
`timeseries_*.csv` retain raw-output-derived histories;
`comparison.pdf` and `numerical_diagnostics.pdf` are vector figures with
matching PNG previews. Exact simulation commands and their input/binary
hashes are in each case's `run.json`; retained input decks are the executable
specification of each run. `../run_study.py` provides generation and execution.
"""
    target=HERE/"STUDY_REPORT.md"
    target.write_text(text)
    if not complete:
        (HERE/"STUDY_REPORT_DRAFT.md").write_text(text)
    elif (HERE/"STUDY_REPORT_DRAFT.md").exists():
        (HERE/"STUDY_REPORT_DRAFT.md").write_text("This draft has been superseded by the completed [study report](STUDY_REPORT.md).\n")
    print(target)


if __name__=="__main__":main()
