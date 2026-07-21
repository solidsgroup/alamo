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
    int nspecies = 1;
    Set::Scalar mw;
    Set::Scalar cp_mass;
    Set::Scalar Tref;

public:
    Rocfire() {} ;

    static void Parse(Rocfire & value, IO::ParmParse & pp, int a_nspecies)
    {
        value.nspecies = a_nspecies;
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
    template<class Composition>
    Set::Scalar cp_mol(Set::Scalar /*T*/, const Composition& X, int i, int j, int k) const
    {
        Set::Scalar cp = 0.0;
        for (int n = 0; n < nspecies; ++n)
        {
            cp += X(i,j,k,n) * cp_mol_species(0.0, n);
        }
        return cp;
    }
    template<class Composition>
    Set::Scalar enthalpy_mol(Set::Scalar T, const Composition& X, int i, int j, int k) const
    {
        Set::Scalar h = 0.0;
        for (int n = 0; n < nspecies; ++n)
        {
            h += X(i,j,k,n) * enthalpy_mol_species(T, n);
        }
        return h;
    }
    template<class Composition>
    Set::Scalar entropy_mol(Set::Scalar /*T*/, const Composition& /*X*/, int /*i*/, int /*j*/, int /*k*/) const
    {
        return 0.0;
    }
    Set::Scalar cp_mol_species(Set::Scalar /*T*/, int /*n*/) const
    {
        return cp_mass * mw;
    }
    Set::Scalar dcp_mol_species_dT(Set::Scalar /*T*/, int /*n*/) const
    {
        return 0.0;
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
