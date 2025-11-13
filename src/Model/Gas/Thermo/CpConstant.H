#ifndef MODEL_GAS_THERMO_CPCONSTANT_H_
#define MODEL_GAS_THERMO_CPCONSTANT_H_

#include <vector>
#include "Model/Gas/Gas.H"

namespace Model {
namespace Gas {
namespace Thermo {

class CpConstant : public virtual Gas {
public:
    std::vector<double> cp_moles;
    std::vector<double> h0;
    std::vector<double> s0;
    std::vector<double> Tref;
    
    CpConstant(
        std::vector<double>& cp_in,
        std::vector<double>& h_in,
        std::vector<double>& s_in,
        std::vector<double>& Tref_in)
        : Gas(0, cp_in) /* Dummy constructor*/, cp_moles(cp_in), h0(h_in), s0(s_in), Tref(Tref_in) {};
    virtual ~CpConstant() = default;

public:
    double cp_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Specific heat (constant pressure), J/(kmol-K)
        double cp = 0.0;
        for (int n=0; n<nspecies; ++n) cp += X(i,j,k,n) * cp_moles[n];
        return cp;
    }        
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Specific enthalpy, J/kmol)
        double h = 0.0;
        for (int n=0; n<nspecies; ++n) h += X(i,j,k,n) * (cp_moles[n] * (T - Tref[n]) + h0[n]);
        return h;
    }
    double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // specific entropy, J/(kmol-K)
        double s = 0.0;
        for (int n=0; n<nspecies; ++n) s += X(i,j,k,n) * (cp_moles[n] * log(T / Tref[n]) + s0[n]);
        return s;
    }

} // class CpConstant

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif
