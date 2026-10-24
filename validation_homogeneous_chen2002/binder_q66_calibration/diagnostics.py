#!/usr/bin/env python3
"""Diagnose density mixing and fit-window/runtime sensitivity; never refit kinetics."""
import argparse
import csv
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
sys.path.insert(0, str(STUDY / "analysis"))
from compare_pure_htpb import closure, digest, table
from report_calibration import audit
from extract_runs import summarize, write_csv
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--long-summary", type=Path)
    args = parser.parse_args()
    out = HERE / "sweep"
    parameters_path = HERE / "parameters_frozen.json"
    parameters = json.loads(parameters_path.read_text())
    rows = json.loads((out / "simulation_summary.json").read_text())
    assert len(rows) == 18
    b, ap = parameters["binder"], parameters["AP"]
    mixing = []
    for r in rows:
        t, q = r["ap_volume_fraction"], r["q_cal_cm2_s"]
        density = (1-t)*b["density_kg_m3"]+t*ap["density_kg_m3"]
        w = t*ap["density_kg_m3"]/density
        mass_averaged_density = (1-w)*b["density_kg_m3"]+w*ap["density_kg_m3"]
        cp = (1-w)*b["cp_J_kg_K"]+w*ap["cp_J_kg_K"]
        Q = (1-w)*b["heat_release_J_kg"]+w*ap["heat_release_J_kg"]
        rhoQ = (1-t)*b["density_kg_m3"]*b["heat_release_J_kg"]+t*ap["density_kg_m3"]*ap["heat_release_J_kg"]
        assert np.isclose(density*Q, rhoQ, rtol=1.e-12)
        material = json.loads((STUDY / "runs" / r["name"] / "case.json").read_text())["reference"]
        assert np.isclose(material["heat_release_J_kg"], Q, rtol=1.e-12)
        old_r, old_T, _ = closure(parameters, t, q)
        new_r, new_T, _ = closure(parameters, t, q, volume_density=True)
        target = r["figure4_solid_r_cm_s"]
        mixing.append(dict(q_cal_cm2_s=q, ap_volume_fraction=t,
            current_density_kg_m3=mass_averaged_density, volume_additive_density_kg_m3=density,
            current_volumetric_cp_J_m3_K=mass_averaged_density*cp,
            volume_additive_volumetric_cp_J_m3_K=density*cp,
            mass_weighted_Q_J_kg=Q, volume_additive_rhoQ_J_m3=rhoQ,
            current_analytic_r_cm_s=old_r, volume_density_analytic_r_cm_s=new_r,
            chen_r_cm_s=target, current_analytic_error_percent=100*(old_r/target-1),
            volume_density_analytic_error_percent=100*(new_r/target-1),
            current_analytic_Ts_K=old_T, volume_density_analytic_Ts_K=new_T))
    write_csv(out / "density_mixing_diagnostic.csv", mixing)
    example = next(r for r in mixing if r["q_cal_cm2_s"] == 500 and r["ap_volume_fraction"] == .4)
    windows = []
    paths = []
    for r in rows:
        matches = list((STUDY / "analysis/binder_q66").glob(f"**/timeseries_{r['name']}.csv"))
        assert len(matches) == 1
        paths.extend(matches)
        raw = [{k: v if k == "plotfile" else float(v) if v else float("nan")
                for k, v in item.items()} for item in csv.DictReader(matches[0].open())]
        for fraction in (.2, .3, .4, .5):
            s = summarize(raw, fraction)
            windows.append(dict(name=r["name"], q_cal_cm2_s=r["q_cal_cm2_s"],
                ap_volume_fraction=r["ap_volume_fraction"], late_fraction=fraction,
                n_late=s["n_late"], rate_cm_s=s["r_cm_s"],
                change_from_default_percent=100*(s["r_cm_s"]/r["r_cm_s"]-1),
                fit_standard_error_cm_s=s["fit_standard_error_cm_s"],
                late_speed_drift_percent=s["late_speed_drift_percent"]))
    write_csv(out / "fit_window_sensitivity.csv", windows)
    maximum = max(windows, key=lambda r: abs(r["change_from_default_percent"]))
    duration_text = "The longer-run comparison is pending; fit-window agreement alone does not prove duration independence."
    duration = None
    if args.long_summary:
        long_rows = json.loads(args.long_summary.read_text())
        assert len(long_rows) == 1
        long = long_rows[0]
        assert long["q_cal_cm2_s"] == 500 and long["ap_volume_fraction"] == .4
        evidence = audit(long_rows, parameters, arrhenius=True,
                         fixed_parameters=HERE / "fixed_properties.json", relaxations=12)
        (out / "duration_execution_audit.json").write_text(json.dumps(evidence, indent=2) + "\n")
        short = next(r for r in rows if r["q_cal_cm2_s"] == 500 and r["ap_volume_fraction"] == .4)
        duration = dict(q_cal_cm2_s=500, ap_volume_fraction=.4,
            short_name=short["name"], long_name=long["name"], short_relaxations=6, long_relaxations=12,
            short_duration_s=short["end_time_s"], long_duration_s=long["end_time_s"],
            short_r_cm_s=short["r_cm_s"], long_r_cm_s=long["r_cm_s"],
            relative_rate_change_percent=100*(long["r_cm_s"]/short["r_cm_s"]-1),
            short_late_drift_percent=short["late_speed_drift_percent"],
            long_late_drift_percent=long["late_speed_drift_percent"],
            short_incident_flux_error_percent=short["incident_flux_error_percent"],
            long_incident_flux_error_percent=long["incident_flux_error_percent"])
        (out / "duration_comparison.json").write_text(json.dumps(duration, indent=2) + "\n")
        paths.append(args.long_summary)
        duration_text = f"""The representative worst-error composition (t=0.4, q=500) was rerun for twelve relaxation times, with identical constituent properties, interface width, spatial resolution, timestep ceiling, heating prescription and other solver options. Duration increased from {duration['short_duration_s']:.8g} to {duration['long_duration_s']:.8g} s. The domain extends six thermal lengths deeper so the final surface retains the same eight-thermal-length separation from the cold boundary; thus this is a combined duration/deep-boundary check, not a same-domain restart.

The measured rate changes from **{duration['short_r_cm_s']:.8f} to {duration['long_r_cm_s']:.8f} cm/s**, a **{duration['relative_rate_change_percent']:+.5f}%** change. Late-window drift is {duration['short_late_drift_percent']:+.5f}% in the six-time run and {duration['long_late_drift_percent']:+.5f}% in the twelve-time run. Incident-flux error also changes from {duration['short_incident_flux_error_percent']:+.5f}% to {duration['long_incident_flux_error_percent']:+.5f}%. Both rate measurements use the final 40% of their respective durations. This supports duration sufficiency for diagnosing a roughly 13.5% blend discrepancy, but does not establish 0.01% accuracy in the asymptotic rate; only this one case has been tested with doubled duration. `duration_comparison.json` and `duration_execution_audit.json` retain the comparison and execution checks."""
    compact = table(["AP volume fraction", "Current analytic error [%]", "Volume-additive density analytic error [%]"],
        [[f"{r['ap_volume_fraction']:.1f}", f"{r['current_analytic_error_percent']:+.3f}",
          f"{r['volume_density_analytic_error_percent']:+.3f}"] for r in mixing if r["q_cal_cm2_s"] == 500])
    report = f'''# Mixing and runtime diagnostics

## Main recommendation: use a volume-additive mixture density

With AP volume fraction t, the volume-additive density is `rho=(1-t)*rho_binder+t*rho_AP`. Equivalently, for AP mass fraction w, `1/rho=(1-w)/rho_binder+w/rho_AP`. Keep cp and Q mass-weighted. This makes `rho*cp=(1-t)*rho_binder*cp_binder+t*rho_AP*cp_AP`, with the corresponding volume-weighted energy density for Q.

The current implementation in `src/Model/Mechanism/PhaseChange.H` instead uses an arithmetic mass-weighted density. At t=0.4 it gives {example['current_density_kg_m3']:.3f} kg/m³ instead of {example['volume_additive_density_kg_m3']:.3f} kg/m³. This raises both rho*cp and rho*Q by {100*(example['current_density_kg_m3']/example['volume_additive_density_kg_m3']-1):.3f}%, altering the energy available per unit regressed volume and lowering the coupled surface temperature.

The calculation below changes **only density in the analytic closure**, retaining the frozen A/E, cp, Q, T0 and AP parameters. At q=500, t=0.4, the predicted surface temperature increases from {example['current_analytic_Ts_K']:.3f} to {example['volume_density_analytic_Ts_K']:.3f} K. The analytic rate error against Chen falls from {example['current_analytic_error_percent']:+.3f}% to {example['volume_density_analytic_error_percent']:+.3f}%.

{compact}

These are diagnostic analytic predictions, not new LowMach runs with a changed density rule. No production code, mixing rule or constituent parameter was changed. The correction would need a consistent implementation in phase change, thermal capacity, initial density and reference/deck generation, followed by a new sweep. The earlier article audit records a density-weighting ambiguity in Chen's text versus its plotted blend results; this diagnostic does not establish which exact properties the authors used. See `../../analysis/REFERENCE_NOTES.md`.

Pure AP is also an unfitted endpoint and is underpredicted by about 3.7–3.9% in the present solver sweep. After resolving density mixing, I recommend calibrating AP A/E separately against the pure-AP endpoints using fixed, supported AP thermal properties and the same held-out strategy. Density mixing alone cannot change a pure endpoint. Only after those checks should residual blend errors motivate a different effective kinetic interpolation; Chen's Eq. (15) supplies a target from the two pure rates, with Eq. (14) supplying its coupled temperature.

## Surface temperature is an output of the prescribed-flux problem

Chen equations (12)–(14) determine r and Ts together:

`q = rho*r*[cp*(Ts-T0)-Q]`, `r = A*exp[-(E/R)/Ts]`.

Prescribing q while solving these equations does not require prescribing Ts independently. Changing density changes that energy balance, then changes Ts and the Arrhenius rate. Thus the observed mixing error is mediated by temperature, but is not evidence that the problem needs an independently prescribed Ts. The measured LowMach/analytic discrepancy in the completed sweep is only {min(r['error_eq14_19_percent'] for r in rows):+.3f}% to {max(r['error_eq14_19_percent'] for r in rows):+.3f}%, much smaller than the largest error against Chen. Incident-flux and interface-width effects remain a separate numerical check. Conductivity controls the thermal profile and relaxation time, but it cancels from this ideal steady surface balance, so changing k to fit steady rates would obscure the issue.

Source: [Chen et al. (2002)](https://doi.org/10.1016/S1540-7489(02)80357-1), user-supplied PDF, printed pp. 2926–2927. The density alternative is our conservation-based diagnostic, not an asserted quotation of their implemented rule.

## Q units, averaging and gas chemistry

The supplied binder Q is −66 cal/g = −276144 J/kg; the retained AP value is −100 cal/g = −418400 J/kg. `PhaseChange.H` parses both with `Unit::Energy()/Unit::Mass()`, averages them as `Q=(1-w)*Q_binder+w*Q_AP`, then applies heat as `-mass_change*Q`. During regression mass_change is negative, so these negative Q values remove heat. The reference generator uses the same mass-specific values and averaging. The gas molecular weight and `system.amount=kmol` do not turn this Q into a molar quantity.

At 40% AP by volume, w={.4*ap['density_kg_m3']/((1-.4)*b['density_kg_m3']+.4*ap['density_kg_m3']):.9f}, giving Q={example['mass_weighted_Q_J_kg']:.6f} J/kg = {example['mass_weighted_Q_J_kg']/4184:.6f} cal/g. This Q averaging is appropriate for additive mass-specific constituent heats and should be retained. With volume-additive density, `rho*Q=(1-t)*rho_binder*Q_binder+t*rho_AP*Q_AP`; the diagnostic verifies this equality for all 18 cases. The present inconsistency is the density used to convert Q into heat per regressed volume, rather than the mass weighting of Q itself.

For a constituent Q given in J/mol, use mole/amount fractions: `Q_molar=sum(x_i*Q_molar_i)`, where `x_i=(t_i*rho_i/M_i)/sum(t_j*rho_j/M_j)`. Volume weighting is equivalent only when the constituents have equal molar density rho_i/M_i. An alternative is to convert each molar Q to J/kg using its constituent molar mass, then apply mass weighting. [IUPAC amount-fraction definition](https://goldbook.iupac.org/terms/view/A00296) and [volume-fraction definition](https://goldbook.iupac.org/terms/view/V06643) distinguish these weights. These formulas assume additive constituent heats without an extra heat of mixing.

Every executed sweep and duration-check deck uses `chemistry.model.type = frozen` and only the binder-regression phase-change mechanism. Gas-phase reactions therefore add no heat; the prescribed gas-slab heat source and the endothermic solid-to-gas Q are the configured thermal sources. The execution audits verify these chemistry settings.

## Rate fitting and runtime sufficiency

The 18 sweep cases use six thermal relaxation times, `tau=delta/r=k/(rho*cp*r²)`. Rates are deterministic planar interface-position slopes, fit over the final 40% (last 2.4 tau), and checked against independent integrated solid-volume loss. The diagnostic is convergence to a steady rate, rather than statistical sampling of randomly packed particles or experimental noise.

Across the completed sweep, maximum absolute first-half/second-half late-window speed drift is {max(abs(r['late_speed_drift_percent']) for r in rows):.5f}%. Recomputing every case with the last 20%, 30%, 40% and 50% of its data changes the measured rate by at most **{abs(maximum['change_from_default_percent']):.5f}%**, at `{maximum['name']}` with a {100*maximum['late_fraction']:g}% fitting window. The windows use {min(r['n_late'] for r in windows)}–{max(r['n_late'] for r in windows)} saved interface positions. `fit_window_sensitivity.csv` records all 72 fits. Their slope standard errors measure numerical fit scatter; they are not confidence intervals for the physical model or Chen's data.

{duration_text}

Small window/runtime sensitivity supports using the reported steady rates. It does not establish spatial convergence, remove the approximately 1.5–2% numerical bias, or imply comparable precision in a figure with about 0.005 cm/s stroke resolution.
'''
    (out / "DIAGNOSTICS.md").write_text(report)
    source_paths = [parameters_path, out / "simulation_summary.json", Path(__file__),
                    STUDY / "analysis/compare_pure_htpb.py", STUDY / "analysis/extract_runs.py",
                    STUDY.parent / "src/Model/Mechanism/PhaseChange.H", *paths]
    (out / "diagnostic_provenance.json").write_text(json.dumps({str(p): digest(p) for p in source_paths}, indent=2) + "\n")
    print(compact)
    print("Maximum fit-window sensitivity [%]:", abs(maximum["change_from_default_percent"]))
    if duration:
        print(json.dumps(duration, indent=2))


if __name__ == "__main__":
    main()
