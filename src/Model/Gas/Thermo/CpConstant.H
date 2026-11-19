#ifndef MODEL_GAS_THERMO_CPCONSTANT_H_
#define MODEL_GAS_THERMO_CPCONSTANT_H_

#include <vector>
#include <memory>
#include <cmath>
#include "IO/ParmParse.H"
#include "Model/Gas/Thermo/Thermo.H"

namespace Model {
namespace Gas {
namespace Thermo {

class CpConstant 
{
public:
    static constexpr const char* name = "cpconstant";

private:
    int nspecies = 1;
    std::vector<double> cp_moles;
    std::vector<double> h0;
    std::vector<double> s0;
    std::vector<double> Tref;
    
public:
    CpConstant() {} ;


    static void Parse(CpConstant & value, IO::ParmParse & pp, int a_nspecies)
    {
        value.nspecies = a_nspecies;

        // heat capacity in moles (one per species)
        pp.queryarr_required("cp_moles",value.cp_moles,Unit::Energy()/Unit::Amount()/Unit::Temperature());
        // molar enthalpy (one per species)
        pp.queryarr_required("h0",value.h0,Unit::Energy()/Unit::Amount());
        // molar entropy (one per species)
        pp.queryarr_required("s0",value.s0, Unit::Energy()/Unit::Amount()/Unit::Temperature());
        // reference temperature (one per species)
        pp.queryarr_required("Tref",value.Tref, Unit::Temperature());

        Util::Assert(INFO, TEST(value.cp_moles.size() == (size_t)value.nspecies));
        Util::Assert(INFO, TEST(value.h0.size() == (size_t)value.nspecies));
        Util::Assert(INFO, TEST(value.s0.size() == (size_t)value.nspecies));
        Util::Assert(INFO, TEST(value.Tref.size() == (size_t)value.nspecies));
    }

public:
    double cp_mol(double /*T*/, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        // Specific heat (constant pressure), J/(kmol-K)
        double cp = 0.0;
        for (int n=0; n<nspecies; ++n) cp += X(i,j,k,n) * cp_moles[n];
        return cp;
    }        
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        // Specific enthalpy, J/kmol)
        double h = 0.0;
        for (int n=0; n<nspecies; ++n) h += X(i,j,k,n) * (cp_moles[n] * (T - Tref[n]) + h0[n]);
        return h;
    }
    double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        // specific entropy, J/(kmol-K)
        double s = 0.0;
        for (int n=0; n<nspecies; ++n) s += X(i,j,k,n) * (cp_moles[n] * log(T / Tref[n]) + s0[n]);
        return s;
    }

}; // class CpConstant

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif
