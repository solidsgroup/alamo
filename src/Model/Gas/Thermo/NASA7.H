#ifndef MODEL_GAS_THERMO_NASA7_H_
#define MODEL_GAS_THERMO_NASA7_H_

#include <vector>
#include <cmath>
#include <memory>
#include "Model/Gas/Thermo/Thermo.H"

namespace Model {
namespace Gas {
namespace Thermo {

class NASA7 : public Thermo {
private:
    std::vector<std::vector<double>> Tranges;
    std::vector<std::vector<std::vector<double>>> coeffs;
    
public:
    NASA7(
        std::shared_ptr<GasData> data_in,
        const std::vector<std::vector<double>>& Tranges_in,
        const std::vector<std::vector<std::vector<double>>>& coeffs_in)
        :   Thermo(data_in), Tranges(Tranges_in), coeffs(coeffs_in) {}
    ~NASA7() override = default;

private:
    // Compute nondimensional quantities (cp/R, h/(RT), s/R) from NASA7 coefficients
    double cp_R(double T, const std::vector<double>& Trange, const std::vector<std::vector<double>>& coeff) const {
        int hilo = (T <= Trange[1]) ? 0 : 1;
        const auto& a = coeff[hilo];
        return a[0] + a[1]*T + a[2]*T*T + a[3]*T*T*T + a[4]*T*T*T*T;
    }   
    
    double h_RT(double T, const std::vector<double>& Trange, const std::vector<std::vector<double>>& coeff) const {
        int hilo = (T <= Trange[1]) ? 0 : 1;
        const auto& a = coeff[hilo];
        return a[0] + 0.5*a[1]*T + (1.0/3.0)*a[2]*T*T + 0.25*a[3]*T*T*T + 0.2*a[4]*T*T*T*T + a[5]/T;
    }

    double s_R(double T, const std::vector<double>& Trange, const std::vector<std::vector<double>>& coeff) const {
        int hilo = (T <= Trange[1]) ? 0 : 1;
        const auto& a = coeff[hilo];
        return a[0]*log(T) + a[1]*T + 0.5*a[2]*T*T + (1.0/3.0)*a[3]*T*T*T + 0.25*a[4]*T*T*T*T + a[6];
    }

public:
    double cp_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Specific heat (constant pressure), J/(kmol-K)
        double cp = 0.0;
        for (int n=0; n<data->nspecies; ++n) cp += X(i,j,k,n) * cp_R(T, Tranges[n], coeffs[n]) * GasData::Ru;
        return cp;
    }        
    double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Specific enthalpy, J/kmol)
        double h = 0.0;
        for (int n=0; n<data->nspecies; ++n) h += X(i,j,k,n) * h_RT(T, Tranges[n], coeffs[n]) * GasData::Ru * T;
        return h;
    }
    double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // specific entropy, J/(kmol-K)
        double s = 0.0;
        for (int n=0; n<data->nspecies; ++n) s += X(i,j,k,n) * s_R(T, Tranges[n], coeffs[n]) * GasData::Ru;
        return s;
    }

}; // class NASA7

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif
