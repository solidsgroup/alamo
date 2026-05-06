// Calorically Perfect Gas

#include "Model/Gas/Gas.H"
#include "Model/Gas/EOS/CPG.H"

namespace Model {
namespace Gas {
namespace EOS {

double CPG::ComputeT(double density, double momentumx, double momentumy, double E, double Tguess,
                Set::Patch<const Set::Scalar>& X, int i, int j, int k, double /*rtol=1e-12*/) const {
    // Temperature, K
    // Since gamma is not a function of temperature, but is a function of composition, we can
    // compute gamma at any dummy temperature
    double P = (E - 0.5*(momentumx*momentumx + momentumy*momentumy)/density) * (gas->gamma(Tguess, X, i, j, k) - 1.0);
    double T = P / density / gas->R(X, i, j, k);
    return T;
}
double CPG::ComputeT_from_primitives(double pressure, double density,
                Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Temperature, K
    double T = pressure / density / gas->R(X, i, j, k);
    return T;
}
double CPG::ComputeP(double density, double T,
                Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Pressure, Pa
    double P = density * gas->R(X, i, j, k) * T;
    return P;
}
double CPG::ComputeE(double density, double momentumx, double momentumy, double T,
                Set::Patch<const Set::Scalar>& X, int i, int j, int k) const {
    // Energy, J/m^3
    double P = density * gas->R(X, i, j, k) * T;
    double rhoE = P / (gas->gamma(T, X, i, j, k) - 1.0);
    double E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy)/density;
    return E;
}

} // namespace EOS
} // namespace Gas
} // namespace Model
