#include <vector>
#include <cmath>
#include <memory>
#include "Util/Util.H"
#include "Set/Base.H"
#include "Set/Set.H"
#include "Model/Gas/Gas.H"
#include "Model/Gas/Thermo/Thermo.H"
#include "Model/Gas/Transport/Transport.H"
#include "Model/Gas/EOS/EOS.H"

namespace Model {
namespace Gas {

// Methods that need to be defined by inherited class

// Thermodynamic quantities
double Gas::cp_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Specific heat (constant pressure), J/(kmol-K)
    if (!thermo) Util::Abort(INFO, "[Gas::cp_mol] No Thermo model attached.");
    return thermo->cp_mol(T, X, i , j, k);
}
double Gas::enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Specific enthalpy, J/kmol)
    if (!thermo) Util::Abort(INFO, "[Gas::enthalpy_mol] No Thermo model attached.");
    return thermo->enthalpy_mol(T, X, i , j, k);
}
double Gas::entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // specific entropy, J/(kmol-K)
    if (!thermo) Util::Abort(INFO, "[Gas::entropy_mol] No Thermo model attached.");
    return thermo->entropy_mol(T, X, i , j, k);
}
double Gas::cp_mol_species(double T, int n) const {
    // Specific heat (constant pressure) for species n, J/(kmol-K)
    if (!thermo) Util::Abort(INFO, "[Gas::cp_mol_species] No Thermo model attached.");
    return thermo->cp_mol_species(T, n);
}

// Transport quantities
double Gas::dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Dynamic viscosity, Pa-s
    if (!transport) Util::Abort(INFO, "[Gas::dynamic_viscosity] No Transport model attached.");
    return transport->dynamic_viscosity(T, X, i , j, k);
}
double Gas::thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Thermal conductivity coefficient, W/(m-K)
    if (!transport) Util::Abort(INFO, "[Gas::thermal_conductivity] No Transport model attached.");
    return transport->thermal_conductivity(T, X, i , j, k);
}
void Gas::diffusion_coeffs(Set::Patch<Set::Scalar>& DKM, double T, double P, Set::Patch<const Set::Scalar>& X, int i, int j, int k) {
    // Species diffusion coefficients, m^2/s
    if (!transport) Util::Abort(INFO, "[Gas::diffusion_coeffs] No Transport model attached.");
    return transport->diffusion_coeffs(DKM, T, P, X, i , j, k);
}

// EOS
double Gas::ComputeT(
        double density, double momentumx, double momentumy, double E, double Tguess,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const {
    // Temperature, K
    if (!eos) Util::Abort(INFO, "[Gas::ComputeT] No EOS model attached.");
    return eos->ComputeT(density, momentumx, momentumy, E, Tguess, X, i , j, k, rtol);
}
double Gas::ComputeT(
        double pressure, double density,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Temperature, K
    if (!eos) Util::Abort(INFO, "[Gas::ComputeT] No EOS model attached.");
    return eos->ComputeT(pressure, density, X, i , j, k);
}
double Gas::ComputeP(double density, double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Pressure, Pa
    if (!eos) Util::Abort(INFO, "[Gas::ComputeP] No EOS model attached.");
    return eos->ComputeP(density, T, X, i , j, k);
}
double Gas::ComputeE(    
        double density, double momentumx, double momentumy, double T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Energy, J/m^3
    if (!eos) Util::Abort(INFO, "[Gas::ComputeE] No EOS model attached.");
    return eos->ComputeE(density, momentumx, momentumy, T, X, i , j, k);
}

} // namespace Gas
} // namespace Model
