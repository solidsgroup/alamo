//
// :code:Model::Gas::EOS::CPG:code: implements a calorically perfect gas equation of
// state for mixture-averaged gas properties.
//
// Pressure, temperature, and conserved energy are related through the ideal-gas law and
// the calorically perfect internal-energy relation using the mixture gas constant and
// heat-capacity ratio supplied by the gas model.
//

#ifndef MODEL_GAS_EOS_CPG_H_
#define MODEL_GAS_EOS_CPG_H_

#include "IO/ParmParse.H"
#include "Model/Gas/EOS/EOS.H"
#include <array>
#include <memory>

namespace Model {
namespace Gas {

namespace EOS {

class CPG
{
public:
    static constexpr const char* name = "cpg";

private:
    std::array<Set::Scalar, NSPECIES> cv_tref_mol = {};

public:
    template<typename THERMO>
    static void Parse(CPG & value, IO::ParmParse & /*pp*/, const THERMO & thermo)
    {
        for (int n=0; n<NSPECIES; ++n)
        {
            const Set::Scalar Tref = thermo.reference_temperature_species(n);
            const Set::Scalar cv_mol =
                thermo.cp_mol_species(Tref, n) - Set::Constant::Rg;
            Util::Assert(INFO, TEST(std::isfinite(Tref) && Tref >= 0.0),
                "EOS::CPG received an invalid reference temperature for species ",
                n, ": ", Tref);
            Util::Assert(INFO, TEST(std::isfinite(cv_mol) && cv_mol > 0.0),
                "EOS::CPG requires positive finite cv for species ", n, ": ", cv_mol);
            value.cv_tref_mol[n] = cv_mol * Tref;
        }
    }

private:
    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar reference_internal_energy_mass(const GAS &gas,
            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        Set::Scalar reference_energy_mol = 0.0;
        for (int n=0; n<NSPECIES; ++n)
            reference_energy_mol += X(i,j,k,n) * cv_tref_mol[n];
        return reference_energy_mol / gas.GetMW(X, i, j, k);
    }

public:
    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar E, Set::Scalar Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, Set::Scalar /*rtol=1e-12*/) const
    {
        const Set::Scalar kinetic =
            0.5*(momentumx*momentumx + momentumy*momentumy)/density;
        const Set::Scalar internal_energy = (E - kinetic) / density
            + reference_internal_energy_mass(gas, X, i, j, k);
        return internal_energy * (gas.gamma(Tguess, X, i, j, k) - 1.0)
            / gas.R(X, i, j, k);
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar E, Set::Scalar Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, Set::Scalar /*rtol=1e-12*/) const
    {
        const Set::Scalar kinetic =
            0.5*(momentumx*momentumx + momentumy*momentumy + momentumz*momentumz)/density;
        const Set::Scalar internal_energy = (E - kinetic) / density
            + reference_internal_energy_mass(gas, X, i, j, k);
        return internal_energy * (gas.gamma(Tguess, X, i, j, k) - 1.0)
            / gas.R(X, i, j, k);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT_from_primitives(Set::Scalar pressure, Set::Scalar density,
                    Set::Scalar R) const
    {
        return pressure / density / R;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeP(Set::Scalar density, Set::Scalar T,
                    Set::Scalar R) const
    {
        return density * R * T;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar T,
                    Set::Scalar R, Set::Scalar gamma) const
    {
        Set::Scalar P = density * R * T;
        Set::Scalar rhoE = P / (gamma - 1.0);
        Set::Scalar E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy)/density;
        return E;
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        return ComputeE(density, momentumx, momentumy, T, gas.R(X, i, j, k),
            gas.gamma(T, X, i, j, k))
            - density * reference_internal_energy_mass(gas, X, i, j, k);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar T,
                    Set::Scalar R, Set::Scalar gamma) const
    {
        Set::Scalar P = density * R * T;
        Set::Scalar rhoE = P / (gamma - 1.0);
        Set::Scalar E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy + momentumz*momentumz)/density;
        return E;
    }

    template<typename GAS>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        return ComputeE(density, momentumx, momentumy, momentumz, T, gas.R(X, i, j, k),
            gas.gamma(T, X, i, j, k))
            - density * reference_internal_energy_mass(gas, X, i, j, k);
    }

}; // class CPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif
