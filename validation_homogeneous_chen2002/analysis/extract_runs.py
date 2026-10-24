#!/usr/bin/env python3
"""Read raw AMReX plotfiles and measure planar solid regression independently.

Usage: extract_runs.py MANIFEST.json [--out DIRECTORY] [--late-fraction .4]
Manifest is a list or {"cases": [...]} with each case containing "name", "output"
(relative to manifest directory or absolute), and arbitrary scalar metadata.
Numerical quantities are interpreted in the solver's SI code units. No analytic
velocity is used to determine a measured surface location or regression speed.
"""
import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/chen-analysis-mpl")
import numpy as np
import yt

yt.set_log_level(40)


def crossing(y, eta, temperature, level):
    indices = np.flatnonzero((eta[:-1]-level)*(eta[1:]-level) <= 0)
    indices = [i for i in indices if eta[i] != eta[i+1]]
    if len(indices) != 1:
        raise ValueError(f"Expected one eta={level} crossing; found {len(indices)}")
    i = indices[0]
    w = (level-eta[i])/(eta[i+1]-eta[i])
    return (1-w)*y[i]+w*y[i+1], (1-w)*temperature[i]+w*temperature[i+1]


def read_plot(path, case=None):
    ds = yt.load(str(path))
    lev = ds.index.max_level
    dims = ds.domain_dimensions * (ds.refine_by**lev)
    if ds.dimensionality < 3:
        dims[ds.dimensionality:] = 1
    grid = ds.covering_grid(lev, ds.domain_left_edge, dims)
    eta = np.asarray(grid["boxlib", "rigid_eta"]).squeeze(axis=2)
    temp = np.asarray(grid["boxlib", "temperature"]).squeeze(axis=2)
    left, right = np.asarray(ds.domain_left_edge), np.asarray(ds.domain_right_edge)
    dy = (right[1]-left[1])/dims[1]
    y = left[1]+(np.arange(dims[1])+.5)*dy
    locations, temperatures = [], []
    for level in (.1,.5,.9):
        values = np.array([crossing(y,e,t,level) for e,t in zip(eta,temp)])
        locations.append(values[:,0])
        temperatures.append(values[:,1])
    grad = np.abs(np.gradient(eta,dy,axis=1))
    row = dict(time_s=float(ds.current_time), plotfile=str(path), dy_m=dy,
               surface_y_m=np.mean(locations[1]),
               surface_roughness_std_m=np.std(locations[1]),
               solid_column_m=np.mean(np.sum(eta,axis=1)*dy),
               surface_temperature_K=np.mean(temperatures[1]),
               surface_temperature_x_std_K=np.std(temperatures[1]),
               temperature_eta01_K=np.mean(temperatures[0]),
               temperature_eta09_K=np.mean(temperatures[2]),
               interface_width_01_09_m=np.mean(np.abs(locations[2]-locations[0])),
               integral_abs_grad_eta=np.mean(np.sum(grad,axis=1)*dy),
               cold_boundary_temperature_K=np.mean(temp[:,0]),
               gas_boundary_temperature_K=np.mean(temp[:,-1]),
               temperature_min_K=np.min(temp),temperature_max_K=np.max(temp),
               eta_min=np.min(eta),eta_max=np.max(eta))
    if case is not None:
        conductivity=np.asarray(grid["boxlib","thermal_conductivity_coeff"]).squeeze(axis=2)
        # Diffusion.cpp requests harmonic cell-to-face averaging. Sampling far
        # outside the constitutive transition avoids interface-flux ambiguity.
        kl,kr=conductivity[:,:-1],conductivity[:,1:]
        kface=np.divide(2*kl*kr,kl+kr,out=np.zeros_like(kl),where=kl+kr>0)
        downward_flux=kface*np.diff(temp,axis=1)/dy
        yface=.5*(y[:-1]+y[1:])
        for level,label in [(1.e-6,"gas"),(.001,"gas_tail"),(.1,"eta01"),(.5,"eta05"),(.9,"eta09"),(.999,"solid")]:
            values=[]
            for e,t,f in zip(eta,temp,downward_flux):
                position,_=crossing(y,e,t,level)
                values.append(np.interp(position,yface,f))
            row[f"conductive_flux_down_{label}_W_m2"]=np.mean(values)
        source_lo,source_hi=case["source_y_m"]
        source_integral=case["reference"]["heat_flux_W_m2"]*np.count_nonzero((y>source_lo)&(y<source_hi))*dy/(source_hi-source_lo)
        row["input_source_integral_W_m2"]=source_integral
        row["bottom_heat_loss_W_m2"]=np.mean(conductivity[:,0]*(temp[:,0]-case["T0_K"])/(.5*dy))
        # The scalar 'energy' field is GAS total energy, not the mixed thermal
        # equation's primitive; reconstruct its actual capacity from eta/T.
        z=np.clip((eta-.11920292202211755)/(.8807970779778824-.11920292202211755),0,1)
        solid_fraction=z*z*z*(z*(6*z-15)+10)
        row["solid_capacity_temperature_gradient_J_m3"]=case["reference"]["density_kg_m3"]*case["reference"]["cp_J_kg_K"]*np.mean(np.sum(solid_fraction*np.gradient(temp,dy,axis=1),axis=1))*dy
        row["_thermal_profile"]={"T":temp,"H":solid_fraction,"dy":dy}
        if case.get("advect_temperature",False):
            u=np.asarray(grid["boxlib","velocityx"]).squeeze(axis=2)
            v=np.asarray(grid["boxlib","velocityy"]).squeeze(axis=2)
            dx=(right[0]-left[0])/dims[0]
            gradx=(np.roll(temp,-1,axis=0)-np.roll(temp,1,axis=0))/(2*dx)
            grady=np.gradient(temp,dy,axis=1)
            gas_factor=case["pressure_Pa"]*case["gas_cp_J_kg_K"]/case.get("gas_constant_J_kg_K",319.7870237751246)
            row["gas_thermal_advection_W_m2"]=gas_factor*np.mean(np.sum((1-solid_fraction)/temp*(u*gradx+v*grady),axis=1))*dy
        else:
            row["gas_thermal_advection_W_m2"]=0.0
    return {k:float(v) if isinstance(v,np.number) else v for k,v in row.items()}


def read_case_plots(plots,case):
    """Add interval-integrated gas/solid thermal storage and energy residual.

    For zero temperature advection, integrate C*dT/dt = div(k grad T)+S+Q*mdot.
    Gas C=p*cp_g/(R*T)*(1-H); the exact logarithmic T increment is combined
    with a trapezoid in H. Solid C=rho*cp_s*H. Time sampling error is assessed
    separately from solver error. These are thermal-equation storage terms,
    not a conservative total-energy integral or a statistical uncertainty.
    """
    rows=[];previous=None
    ref=case["reference"]
    rho,cp,Q=(ref[k] for k in ("density_kg_m3","cp_J_kg_K","heat_release_J_kg"))
    gas_factor=case["pressure_Pa"]*case["gas_cp_J_kg_K"]/case.get("gas_constant_J_kg_K",319.7870237751246)
    for path in plots:
        row=read_plot(path,case)
        profile=row.pop("_thermal_profile")
        if previous is not None:
            old,old_profile=previous
            dt=row["time_s"]-old["time_s"]
            if dt <= 0:continue
            H=.5*(profile["H"]+old_profile["H"])
            dT=profile["T"]-old_profile["T"]
            integrate=lambda x:np.mean(np.sum(x,axis=1))*profile["dy"]/dt
            gas=gas_factor*integrate((1-H)*np.log(profile["T"]/old_profile["T"]))
            solid=rho*cp*integrate(H*dT)
            latent=Q*rho*(old["solid_column_m"]-row["solid_column_m"])/dt
            source=.5*(row["input_source_integral_W_m2"]+old["input_source_integral_W_m2"])
            bottom=.5*(row["bottom_heat_loss_W_m2"]+old["bottom_heat_loss_W_m2"])
            advection=.5*(row["gas_thermal_advection_W_m2"]+old["gas_thermal_advection_W_m2"])
            row.update(storage_interval_start_s=old["time_s"],
                gas_thermal_storage_W_m2=gas,solid_thermal_storage_W_m2=solid,
                latent_heat_source_W_m2=latent,
                gas_storage_percent_input=100*gas/source,
                thermal_balance_residual_W_m2=source+latent-bottom-gas-solid-advection,
                thermal_balance_residual_percent_input=100*(source+latent-bottom-gas-solid-advection)/source)
        rows.append(row);previous=row,profile
    return rows


def fit_speed(rows, quantity):
    t = np.array([r["time_s"] for r in rows])
    s = np.array([r[quantity] for r in rows])
    if len(t) < 3 or np.ptp(t) == 0:
        raise ValueError("Need at least three distinct late-time plotfiles")
    coeff = np.polyfit(t-t.mean(),s,1)
    residual = s-np.polyval(coeff,t-t.mean())
    sse = np.sum(residual**2)
    stderr = np.sqrt(sse/(len(t)-2)/np.sum((t-t.mean())**2))
    return -coeff[0], stderr, np.sqrt(np.mean(residual**2))


def summarize(rows, late_fraction):
    rows = sorted(rows,key=lambda r:r["time_s"])
    cutoff = rows[-1]["time_s"]-late_fraction*(rows[-1]["time_s"]-rows[0]["time_s"])
    late = [r for r in rows if r["time_s"] >= cutoff]
    r,stderr,rms = fit_speed(late,"surface_y_m")
    volume,_,_ = fit_speed(late,"solid_column_m")
    # Compare independent halves of the late interval. These measure drift,
    # not statistical uncertainty; the regression stderr is also descriptive.
    mid=(late[0]["time_s"]+late[-1]["time_s"])/2
    halves=[[p for p in late if p["time_s"] <= mid],
            [p for p in late if p["time_s"] >= mid]]
    speeds=[fit_speed(h,"surface_y_m")[0] if len(h)>=3 else float("nan") for h in halves]
    result=dict(n_plotfiles=len(rows),n_late=len(late),end_time_s=rows[-1]["time_s"],
                fit_start_s=late[0]["time_s"],fit_end_s=late[-1]["time_s"],
                r_m_s=r,r_cm_s=100*r,fit_standard_error_cm_s=100*stderr,
                fit_position_residual_rms_m=rms,volume_r_cm_s=100*volume,
                late_first_half_r_cm_s=100*speeds[0],late_second_half_r_cm_s=100*speeds[1],
                late_speed_drift_percent=100*(speeds[1]-speeds[0])/r,
                surface_temperature_K=np.mean([p["surface_temperature_K"] for p in late]),
                surface_temperature_range_K=np.ptp([p["surface_temperature_K"] for p in late]),
                temperature_eta01_K=np.mean([p["temperature_eta01_K"] for p in late]),
                temperature_eta09_K=np.mean([p["temperature_eta09_K"] for p in late]),
                integral_abs_grad_eta=np.mean([p["integral_abs_grad_eta"] for p in late]),
                final_solid_depth_m=rows[-1]["solid_column_m"],
                max_surface_roughness_std_m=max(p["surface_roughness_std_m"] for p in rows))
    energy_keys=["conductive_flux_down_gas_W_m2","conductive_flux_down_solid_W_m2",
        "conductive_flux_down_gas_tail_W_m2",
        "conductive_flux_down_eta01_W_m2","conductive_flux_down_eta05_W_m2","conductive_flux_down_eta09_W_m2",
        "input_source_integral_W_m2","bottom_heat_loss_W_m2","gas_thermal_storage_W_m2",
        "solid_thermal_storage_W_m2","latent_heat_source_W_m2","gas_storage_percent_input",
        "gas_thermal_advection_W_m2",
        "solid_capacity_temperature_gradient_J_m3",
        "thermal_balance_residual_W_m2","thermal_balance_residual_percent_input"]
    for key in energy_keys:
        values=[p[key] for p in late if key in p and np.isfinite(p[key])]
        if values:result[key]=np.mean(values)
    if "input_source_integral_W_m2" in result:
        result["incident_flux_error_percent"]=100*(result["conductive_flux_down_gas_W_m2"]/result["input_source_integral_W_m2"]-1)
        if "solid_capacity_temperature_gradient_J_m3" in result:
            steady_solid=result["r_m_s"]*result["solid_capacity_temperature_gradient_J_m3"]
            residual=result["input_source_integral_W_m2"]+result["latent_heat_source_W_m2"]-result["bottom_heat_loss_W_m2"]-result["gas_thermal_storage_W_m2"]-steady_solid-result["gas_thermal_advection_W_m2"]
            result["steady_profile_solid_storage_W_m2"]=steady_solid
            result["steady_profile_balance_residual_percent_input"]=100*residual/result["input_source_integral_W_m2"]
    return result


def write_csv(path, rows):
    if not rows:
        return
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader();w.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest",type=Path)
    parser.add_argument("--out",type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument("--late-fraction",type=float,default=.4)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads(args.manifest.read_text())
    cases=manifest["cases"] if isinstance(manifest,dict) else manifest
    summaries, errors=[],[]
    for case in cases:
        output=Path(case["output"])
        if not output.is_absolute(): output=args.manifest.resolve().parent/output
        plots=sorted(p.parent for p in output.glob("*cell/Header"))
        try:
            rows=[read_plot(p) for p in plots]
            if not rows: raise ValueError("No plotfiles")
            # Restart/final writers may emit duplicate physical times.
            rows=sorted({r["time_s"]:r for r in rows}.values(),key=lambda r:r["time_s"])
            write_csv(args.out/f"timeseries_{case['name']}.csv",rows)
            metadata={k:v for k,v in case.items() if isinstance(v,(str,int,float,bool))}
            summaries.append(metadata|summarize(rows,args.late_fraction))
            print(case["name"],summaries[-1]["r_cm_s"],flush=True)
        except Exception as exc:
            errors.append(dict(name=case["name"],error=str(exc)))
            print(case["name"],"ERROR",str(exc),flush=True)
    write_csv(args.out/"simulation_summary.csv",summaries)
    (args.out/"simulation_summary.json").write_text(json.dumps(summaries,indent=2)+"\n")
    (args.out/"extraction_errors.json").write_text(json.dumps(errors,indent=2)+"\n")
    if errors: raise SystemExit(1)


if __name__ == "__main__":
    main()
