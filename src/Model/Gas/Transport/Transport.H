#ifndef MODEL_GAS_TRANSPORT_H_
#define MODEL_GAS_TRANSPORT_H_

#include <vector>
#include <memory>
#include <cmath>

namespace Model {
namespace Gas {
namespace Transport {

class Transport {
protected:
    std::shared_ptr<GasData> data;
    std::shared_ptr<Thermo::Thermo> thermo;

public:
    Transport(
            std::shared_ptr<GasData> data_in,
            std::shared_ptr<Thermo::Thermo> thermo_in) 
        : data(data_in), thermo(thermo_in) {}
    virtual ~Transport() = default;

    // Mixture dynamic viscosity, Pa-s
    virtual double dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;
    // Mixture thermal conductivity coefficient, W/(m-K)
    virtual double thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;
    // Mixture species diffusion coefficients, m^2/s
    virtual void diffusion_coeffs(
        Set::Patch<Set::Scalar>& DKM, double T, double P,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) = 0;

}; // class Transport

} // namespace Transport
} // namespace Gas
} // namespace Model

#endif
