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
Set::Scalar Gas::cp_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Specific heat (constant pressure), J/(kmol-K)
    return thermo.cp_mol(T, X, i , j, k);
}
Set::Scalar Gas::enthalpy_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Specific enthalpy, J/kmol)
    return thermo.enthalpy_mol(T, X, i , j, k);
}
Set::Scalar Gas::entropy_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // specific entropy, J/(kmol-K)
    return thermo.entropy_mol(T, X, i , j, k);
}
Set::Scalar Gas::cp_mol_species(Set::Scalar T, int n) const {
    // Specific heat (constant pressure) for species n, J/(kmol-K)
    return thermo.cp_mol_species(T, n);
}
Set::Scalar Gas::enthalpy_mol_species(Set::Scalar T, int n) const {
    // Specific enthalpy for species n, J/kmol
    return thermo.enthalpy_mol_species(T, n);
}
Set::Scalar Gas::entropy_mol_species(Set::Scalar T, int n) const {
    // specific entropy for species n, J/(kmol-K)
    return thermo.entropy_mol_species(T, n);
}

// Transport quantities
Set::Scalar Gas::dynamic_viscosity(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Dynamic viscosity, Pa-s
    return transport.dynamic_viscosity(T, X, i , j, k);
}
Set::Scalar Gas::thermal_conductivity(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Thermal conductivity coefficient, W/(m-K)
    return transport.thermal_conductivity(T, X, i , j, k);
}
void Gas::diffusion_coeffs(Set::Patch<Set::Scalar>& DKM, Set::Scalar T, Set::Scalar P, Set::Patch<const Set::Scalar>& X, int i, int j, int k) {
    // Species diffusion coefficients, m^2/s
    return transport.diffusion_coeffs(DKM, T, P, X, i , j, k);
}

// EOS
#if AMREX_SPACEDIM == 2
Set::Scalar Gas::ComputeT(
        Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar E, Set::Scalar Tguess,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k, Set::Scalar rtol) const
{
    // Temperature, K
    return eos.ComputeT(*this, density, momentumx, momentumy, E, Tguess, X, i, j, k, rtol);
}
Set::Scalar Gas::ComputeE(
        Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
{
    // Energy, J/m^3
    return eos.ComputeE(*this, density, momentumx, momentumy, T, X, i, j, k);
}
#elif AMREX_SPACEDIM == 3
Set::Scalar Gas::ComputeT(
        Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar E, Set::Scalar Tguess,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k, Set::Scalar rtol) const
{
    // Temperature, K
    return eos.ComputeT(*this, density, momentumx, momentumy, momentumz, E, Tguess, X, i, j, k, rtol);
}
Set::Scalar Gas::ComputeE(
        Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
{
    // Energy, J/m^3
    return eos.ComputeE(*this, density, momentumx, momentumy, momentumz, T, X, i, j, k);
}
#endif
Set::Scalar Gas::ComputeT_from_primitives(
        Set::Scalar pressure, Set::Scalar density,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
{
    // Temperature, K
    return eos.ComputeT_from_primitives(pressure, density, R(X,i,j,k));
}
Set::Scalar Gas::ComputeP(Set::Scalar density, Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
{
    // Pressure, Pa
    return eos.ComputeP(density, T, R(X,i,j,k));
}

} // namespace Gas
} // namespace Model
