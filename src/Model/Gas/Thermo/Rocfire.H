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

#include <array>
#include <vector>
#include <memory>
#include <cmath>
#include <string>
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
    Set::Scalar cp_mass;
    Set::Scalar Tref;

public:
    Rocfire() {} ;

    static void Parse(Rocfire & value, IO::ParmParse & pp)
    {
        // gravimetric specific heat
        pp.query_default("cp_mass", value.cp_mass, "0.3_cal/g/K", Unit::SpecificHeatCapacity());
        // reference temperature
        pp.query_default("Tref", value.Tref, "300.0_K", Unit::Temperature());

        const std::string suffix = ".thermo.rocfire";
        const std::string prefix = pp.getPrefix();
        Util::Assert(INFO,
            TEST(prefix.size() > suffix.size() &&
                 prefix.compare(prefix.size() - suffix.size(), suffix.size(), suffix) == 0),
            "Thermo::Rocfire expected parser prefix ending in ", suffix,
            ", got ", prefix);
        IO::ParmParse gas_pp(prefix.substr(0, prefix.size() - suffix.size()));
        std::array<Set::Scalar, NSPECIES> MW;
        gas_pp.queryarr_required("mw", MW, Unit::Mass()/Unit::Amount());

        value.mw = MW[0];
        const Set::Scalar mw_tol = 1.0e-12 * ((std::fabs(value.mw) > 1.0) ? std::fabs(value.mw) : 1.0);
        for (int n=0; n<NSPECIES; ++n)
        {
            Util::Assert(INFO, TEST(std::isfinite(MW[n]) && MW[n] > 0.0),
                "Invalid molecular weight for species ", n, ": ", MW[n]);
            if (std::fabs(MW[n] - value.mw) > mw_tol)
            {
                Util::ParmParseException(INFO, pp.getPrefix(),
                    "thermo.type=rocfire requires all gas mw entries to be equal.");
            }
        }

        Util::Assert(INFO, TEST(value.mw > 0.0), "gas mw entries must be positive.");
        Util::Assert(INFO, TEST(value.cp_mass > 0.0), "gas.thermo.rocfire.cp_mass must be positive.");
        Util::Assert(INFO, TEST(value.Tref >= 0.0), "gas.thermo.rocfire.Tref must be greater than or equal to zero.");
    }

public:
    AMREX_FORCE_INLINE
    Set::Scalar cp_mol(Set::Scalar /*T*/, Set::Patch<const Set::Scalar>& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return cp_mol_species(0.0, 0);
    }

    AMREX_FORCE_INLINE
    Set::Scalar enthalpy_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return enthalpy_mol_species(T, 0);
    }

    AMREX_FORCE_INLINE
    Set::Scalar entropy_mol(Set::Scalar /*T*/, Set::Patch<const Set::Scalar>& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return 0.0;
    }

    AMREX_FORCE_INLINE
    Set::Scalar cp_mol_species(Set::Scalar /*T*/, int /*n*/) const
    {
        return cp_mass * mw;
    }

    AMREX_FORCE_INLINE
    Set::Scalar enthalpy_mol_species(Set::Scalar T, int /*n*/) const
    {
        return cp_mass * mw * (T - Tref);
    }

    AMREX_FORCE_INLINE
    Set::Scalar entropy_mol_species(Set::Scalar /*T*/, int /*n*/) const
    {
        return 0.0;
    }

}; // class Rocfire

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif
