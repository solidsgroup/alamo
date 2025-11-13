#ifndef MODEL_GAS_EOS_CPG_H_
#define MODEL_GAS_EOS_CPG_H_

namespace Model {
namespace Gas {
namespace EOS {

class CPG : public Gas {
public:
    CPG() {};
    ~CPG() = default;

public:
    double ComputeT(double density, double momentumx, double momentumy, double E, double Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, double /*rtol=1e-12*/) const override {
        // Temperature, K
        // Since gamma is not a function of temperature, but is a function of composition, we can
        // compute gamma at any dummy temperature
        double P = (E - 0.5*(momentumx*momentumx + momentumy*momentumy)/density) * (gamma(Tguess, X, i, j, k) - 1.0);
        double T = P / density / R(X, i, j, k);
        return T;
    }
    double ComputeP(double density, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Pressure, Pa
        double P = density * R(X, i, j, k) * T;
        return P;
    }
    double ComputeE(double density, double momentumx, double momentumy, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Energy, J/m^3
        double P = density * R(X, i, j, k) * T;
        double rhoE = P / (gamma(T, X, i, j, k) - 1.0);
        double E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy)/density;
        return E;
    }

} // class CPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif
