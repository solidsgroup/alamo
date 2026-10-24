"""Reconstruct this study's instantaneous gas heat from saved physical fields.

Mirrors Rocfire::compute_progress_rates/source_from_progress and
LowMach::UpdateDerivedDiagnostics. This permits omitting redundant diagnostic
arrays during time stepping. Validate against solver qdot before use.
"""
import numpy as np

SPECIES=('AP_gas','HTPB_gas','Mono','Premixed','Primary','Final')
RG=8.31446261815


def heat_release(data, meta):
    model=meta['gas_model']
    temp=data['boxlib','temperature'].d
    rho=np.stack([np.maximum(data['boxlib','component_density_'+s].d,0) for s in SPECIES])
    total=rho.sum(axis=0)
    fractions=rho/np.maximum(total,1e-300)
    P=meta['pressure_Pa']
    mw=model['molecular_weight_g_mol']*.001
    gas_fraction=np.clip(total*RG/mw*temp/P,0,1)
    condensed=(np.maximum(data['boxlib','component_density_AP_solid'].d,0)/1950+
               np.maximum(data['boxlib','component_density_HTPB_solid'].d,0)/meta['material_properties']['matrix_density_kg_m3'])
    eta_lo=.11920292202211755
    z=np.clip((condensed-eta_lo)/(1-2*eta_lo),0,1)
    phase_H=z**3*(z*(6*z-15)+10)
    occupancy=np.maximum(0,np.minimum(gas_fraction,1-phase_H))
    A=np.array(model['A_g_cm3_s'])*1000
    n=np.array(model['pressure_exponent'])
    energy=np.array(model['activation_energy_kcal_mol'])*4184
    heat=np.array(model['qgas_cal_g'])*4184
    heat[:2]+=np.array(model['qsolid_cal_g'])*4184
    powers=np.stack((fractions[0],fractions[1],fractions[0]*fractions[1],fractions[2]*fractions[3]))
    shape=(4,)+(1,)*temp.ndim
    progress=(A*(P/model['pressure_reference_Pa'])**n).reshape(shape)*np.exp(-energy.reshape(shape)/(RG*np.maximum(temp,1e-300)))*powers
    return occupancy*np.sum(heat.reshape(shape)*progress,axis=0)


def stored_or_reconstructed_heat(data, meta):
    if ('boxlib','qdot') in data.ds.field_list:
        return data['boxlib','qdot'].d
    return heat_release(data,meta)
