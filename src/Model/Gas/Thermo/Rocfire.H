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
    Set::Scalar cp_mass;
    Set::Scalar Tref;

public:
    Rocfire() {} ;

    static void Parse(Rocfire & value, IO::ParmParse & pp)
    {
        // molecular weight
        pp.query_default("mw", value.mw, "26.0_g/mol", Unit::Mass()/Unit::Amount());
        // gravimetric specific heat
        pp.query_default("cp_mass", value.cp_mass, "0.3_cal/g/K", Unit::SpecificHeatCapacity());
        // reference temperature
        pp.query_default("Tref", value.Tref, "300.0_K", Unit::Temperature());

        Util::Assert(INFO, TEST(value.mw > 0.0), "gas.thermo.rocfire.mw must be positive.");
        Util::Assert(INFO, TEST(value.cp_mass > 0.0), "gas.thermo.rocfire.cp_mass must be positive.");
        Util::Assert(INFO, TEST(value.Tref >= 0.0), "gas.thermo.rocfire.Tref must be greater than or equal to zero.");
    }

public:
    Set::Scalar cp_mol(Set::Scalar /*T*/, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        Set::Scalar cp = 0.0;
        for (int n = 0; n < NSPECIES; ++n)
        {
            cp += X(i,j,k,n) * cp_mol_species(0.0, n);
        }
        return cp;
    }
    Set::Scalar enthalpy_mol(Set::Scalar T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        Set::Scalar h = 0.0;
        for (int n = 0; n < NSPECIES; ++n)
        {
            h += X(i,j,k,n) * enthalpy_mol_species(T, n);
        }
        return h;
    }
    Set::Scalar entropy_mol(Set::Scalar /*T*/, Set::Patch<const Set::Scalar>& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return 0.0;
    }
    Set::Scalar cp_mol_species(Set::Scalar /*T*/, int /*n*/) const
    {
        return cp_mass * mw;
    }
    Set::Scalar enthalpy_mol_species(Set::Scalar T, int /*n*/) const
    {
        return cp_mass * mw * (T - Tref);
    }
    Set::Scalar entropy_mol_species(Set::Scalar /*T*/, int /*n*/) const
    {
        return 0.0;
    }

}; // class Rocfire

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif
