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
    return thermo.cp_mol(T, X, i , j, k);
}
double Gas::enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Specific enthalpy, J/kmol)
    return thermo.enthalpy_mol(T, X, i , j, k);
}
double Gas::entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // specific entropy, J/(kmol-K)
    return thermo.entropy_mol(T, X, i , j, k);
}
double Gas::cp_mol_species(double T, int n) const {
    // Specific heat (constant pressure) for species n, J/(kmol-K)
    return thermo.cp_mol_species(T, n);
}
double Gas::enthalpy_mol_species(double T, int n) const {
    // Specific enthalpy for species n, J/kmol
    return thermo.enthalpy_mol_species(T, n);
}
double Gas::entropy_mol_species(double T, int n) const {
    // specific entropy for species n, J/(kmol-K)
    return thermo.entropy_mol_species(T, n);
}

// Transport quantities
double Gas::dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Dynamic viscosity, Pa-s
    return transport.dynamic_viscosity(T, X, i , j, k);
}
double Gas::thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Thermal conductivity coefficient, W/(m-K)
    return transport.thermal_conductivity(T, X, i , j, k);
}
void Gas::diffusion_coeffs(Set::Patch<Set::Scalar>& DKM, double T, double P, Set::Patch<const Set::Scalar>& X, int i, int j, int k) {
    // Species diffusion coefficients, m^2/s
    return transport.diffusion_coeffs(DKM, T, P, X, i , j, k);
}

// EOS
#if AMREX_SPACEDIM == 2
double Gas::ComputeT(
        double density, double momentumx, double momentumy, double E, double Tguess,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const 
{
    // Temperature, K
    return eos.ComputeT(*this, density, momentumx, momentumy, E, Tguess, X, i, j, k, rtol);
}
double Gas::ComputeE(    
        double density, double momentumx, double momentumy, double T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const 
{
    // Energy, J/m^3
    return eos.ComputeE(density, momentumx, momentumy, T, R(X,i,j,k), gamma(T,X,i,j,k));
}
#elif AMREX_SPACEDIM == 3
double Gas::ComputeT(
        double density, double momentumx, double momentumy, double momentumz, double E, double Tguess,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const 
{
    // Temperature, K
    return eos.ComputeT(*this, density, momentumx, momentumy, momentumz, E, Tguess, X, i, j, k, rtol);
}
double Gas::ComputeE(    
        double density, double momentumx, double momentumy, double momentumz, double T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const 
{
    // Energy, J/m^3
    return eos.ComputeE(density, momentumx, momentumy, momentumz, T, R(X,i,j,k), gamma(T,X,i,j,k));
}
#endif
double Gas::ComputeT_from_primitives(
        double pressure, double density,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const 
{
    // Temperature, K
    return eos.ComputeT_from_primitives(pressure, density, R(X,i,j,k));
}
double Gas::ComputeP(double density, double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const 
{
    // Pressure, Pa
    return eos.ComputeP(density, T, R(X,i,j,k));
}

} // namespace Gas
} // namespace Model
