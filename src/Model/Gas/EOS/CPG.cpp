#include "Model/Gas/Gas.H"
#include "Model/Gas/EOS/CPG.H"

namespace Model {
namespace Gas {
namespace EOS {

double CPG::ComputeT(double density, double momentumx, double momentumy, double E, 
                double R, double gamma, double /*rtol=1e-12*/) const {
    // Temperature, K
    // Since gamma is not a function of temperature, but is a function of composition, we can
    // compute gamma at any dummy temperature
    double P = (E - 0.5*(momentumx*momentumx + momentumy*momentumy)/density) * (gamma - 1.0);
    double T = P / density / R;
    return T;
}
double CPG::ComputeT(double pressure, double density,
                double R) const {
    // Temperature, K
    double T = pressure / density / R;
    return T;
}
double CPG::ComputeP(double density, double T,
                double R) const {
    // Pressure, Pa
    double P = density * R * T;
    return P;
}
double CPG::ComputeE(double density, double momentumx, double momentumy, double T,
                double R, double gamma) const {
    // Energy, J/m^3
    double P = density * R * T;
    double rhoE = P / (gamma - 1.0);
    double E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy)/density;
    return E;
}

} // namespace EOS
} // namespace Gas
} // namespace Model
