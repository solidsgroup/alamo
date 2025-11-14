#ifndef MODEL_GAS_THERMO_CPCONSTANT_H_
#define MODEL_GAS_THERMO_CPCONSTANT_H_

#include <vector>
#include <memory>
#include <cmath>
#include "Model/Gas/Thermo/Thermo.H"

namespace Model {
namespace Gas {
namespace Thermo {

class CpConstant : public Thermo {
private:
    std::vector<double> cp_moles;
    std::vector<double> h0;
    std::vector<double> s0;
    std::vector<double> Tref;
    
public:
    CpConstant(
        std::shared_ptr<GasData> data_in,
        const std::vector<double>& cp_in,
        const std::vector<double>& h_in,
        const std::vector<double>& s_in,
        const std::vector<double>& Tref_in)
        : Thermo(data_in), cp_moles(cp_in), h0(h_in), s0(s_in), Tref(Tref_in) {}
    ~CpConstant() override = default;

public:
    double cp_mol(double /*T*/, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Specific heat (constant pressure), J/(kmol-K)
        double cp = 0.0;
        for (int n=0; n<data->nspecies; ++n) cp += X(i,j,k,n) * cp_moles[n];
        return cp;
    }        
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Specific enthalpy, J/kmol)
        double h = 0.0;
        for (int n=0; n<data->nspecies; ++n) h += X(i,j,k,n) * (cp_moles[n] * (T - Tref[n]) + h0[n]);
        return h;
    }
    double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // specific entropy, J/(kmol-K)
        double s = 0.0;
        for (int n=0; n<data->nspecies; ++n) s += X(i,j,k,n) * (cp_moles[n] * log(T / Tref[n]) + s0[n]);
        return s;
    }

}; // class CpConstant

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif
