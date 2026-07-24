//
// :code:Model::Gas::Thermo::Rocfire:code: provides a reduced thermodynamic model used
// with Rocfire-style propellant gas closures.
//
// The model assumes a constant species molecular weight and constant mass heat capacity.
// Species enthalpy is linear in temperature relative to :code:Tref:code:, and entropy is
// treated as zero for this reduced closure.
//
// Reference: Gross, M. L., Hedman, T. D., Son, S. F., Jackson, T. L., &
// Beckstead, M. W. (2013). Coupling micro and meso-scale combustion models of
// AP/HTPB propellants. Combustion and Flame, 160(5), 982-992.
// https://doi.org/10.1016/j.combustflame.2013.01.016
//

#ifndef MODEL_GAS_THERMO_ROCFIRE_H_
#define MODEL_GAS_THERMO_ROCFIRE_H_

#include <vector>
#include <memory>
#include <cmath>
#include "IO/ParmParse.H"
#include "Model/Gas/Thermo/Thermo.H"

namespace Model {
namespace Gas {
namespace Thermo {

class Rocfire
{
public:
    static constexpr const char* name = "rocfire";

private:
    Set::Scalar mw;
    Set::Scalar cp_mass_value;
    Set::Scalar Tref;

public:
    Rocfire() {} ;

    static void Parse(Rocfire & value, IO::ParmParse & pp,
            const std::vector<double>& molecular_weight)
    {
        Util::Assert(INFO, TEST(!molecular_weight.empty()),
            "gas.mw must contain at least one molecular weight.");
        value.mw = molecular_weight[0];
        // gravimetric specific heat
        pp.query_default("cp_mass", value.cp_mass_value, "0.3_cal/g/K", Unit::SpecificHeatCapacity());
        // reference temperature
        pp.query_default("Tref", value.Tref, "300.0_K", Unit::Temperature());

        const Set::Scalar mw_tol = 1.0e-12 *
            ((std::fabs(value.mw) > 1.0) ? std::fabs(value.mw) : 1.0);
        for (std::size_t n = 0; n < molecular_weight.size(); ++n)
        {
            Util::Assert(INFO,
                TEST(std::isfinite(molecular_weight[n]) && molecular_weight[n] > 0.0),
                "Invalid molecular weight for species ", n, ": ", molecular_weight[n]);
            Util::Assert(INFO,
                TEST(std::fabs(molecular_weight[n] - value.mw) <= mw_tol),
                "gas.thermo.rocfire requires all gas molecular weights to equal ",
                value.mw, "; species ", n, " has molecular weight ",
                molecular_weight[n], ".");
        }

        Util::Assert(INFO, TEST(value.cp_mass_value > 0.0), "gas.thermo.rocfire.cp_mass must be positive.");
        Util::Assert(INFO, TEST(value.Tref >= 0.0), "gas.thermo.rocfire.Tref must be greater than or equal to zero.");
    }

public:
    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar cp_mol(Set::Scalar /*T*/, const Composition& /*X*/,
            int /*i*/, int /*j*/, int /*k*/) const
    {
        return cp_mol_species(0.0, 0);
    }
    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar enthalpy_mol(Set::Scalar T, const Composition& /*X*/,
            int /*i*/, int /*j*/, int /*k*/) const
    {
        return enthalpy_mol_species(T, 0);
    }
    template<class Composition>
    Set::Scalar entropy_mol(Set::Scalar /*T*/, const Composition& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return 0.0;
    }
    Set::Scalar cp_mol_species(Set::Scalar /*T*/, int /*n*/) const
    {
        return cp_mass_value * mw;
    }
    Set::Scalar dcp_mol_species_dT(Set::Scalar /*T*/, int /*n*/) const
    {
        return 0.0;
    }
    Set::Scalar enthalpy_mol_species(Set::Scalar T, int /*n*/) const
    {
        return cp_mass_value * mw * (T - Tref);
    }
    Set::Scalar entropy_mol_species(Set::Scalar /*T*/, int /*n*/) const
    {
        return 0.0;
    }

    template<class Fraction>
    AMREX_FORCE_INLINE
    Fraction MoleFractions(
            const amrex::Array4<const Set::Scalar>& partial_density,
            int a_nspecies, int i, int j, int k) const
    {
        return Fraction(partial_density, a_nspecies, i, j, k);
    }

    template<class Composition>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar GetMW(const Composition& /*X*/,
            int /*i*/, int /*j*/, int /*k*/) const
    {
        return mw;
    }

}; // class Rocfire

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif
