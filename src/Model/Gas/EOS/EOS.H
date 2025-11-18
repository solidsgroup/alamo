#ifndef MODEL_GAS_EOS_H_
#define MODEL_GAS_EOS_H_

#include <vector>
#include <memory>
#include <cmath>

namespace Model {
namespace Gas {
namespace EOS {

class EOS {
protected:
    std::shared_ptr<GasData> data;
    std::shared_ptr<Gas> gas;

public:
    EOS(std::shared_ptr<GasData> data_in, std::shared_ptr<Gas> gas_in) : data(data_in), gas(gas_in) {}
    virtual ~EOS() = default;

    // Temperature, K
    virtual double ComputeT(
        double density, double momentumx, double momentumy, double E, double Tguess,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const = 0;
    virtual double ComputeT(
        double pressure, double density,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;
    // Pressure, Pa
    virtual double ComputeP(
        double density, double T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;
    // Energy, J/m^3
    virtual double ComputeE(
        double density, double momentumx, double momentumy, double T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;

}; // class EOS

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif
