#ifndef MODEL_FLUID_H_
#define MODEL_FLUID_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Set/Set.H"
#include <cmath>
#include <limits>  // For NaN

namespace Model {
namespace Fluid {

class Fluid {
public:
    Fluid() {}

    virtual ~Fluid() {}

    template <typename RealType>
    inline RealType ComputePressureFromDensityAndInternalEnergy(
        RealType rho, RealType e, RealType gamma) const {
        if (rho <= 0 || e <= 0 || gamma <= 1.0) {
            Util::Message(INFO, "Invalid values in ComputePressureFromDensityAndInternalEnergy: " + std::to_string(rho));
            return 0.0;
        }
        return (gamma - 1.0) * e;
    }

    template <typename RealType>
    inline RealType ComputeWaveSpeed(
        RealType rho, RealType p, RealType gamma) const {
        if (rho <= 0 || p < 0 || gamma <= 1.0) {
            Util::Message(INFO, "Invalid values in ComputeWaveSpeed: " + std::to_string(p));
            return 0.0;
        }
        return std::sqrt(gamma * p / rho);
    }

    template <typename RealType>
    inline RealType ComputeInternalEnergyFromDensityAndPressure(
        RealType rho, RealType p, RealType gamma) const {
        if (rho <= 0 || gamma <= 1.0) {
            Util::Message(INFO, "Invalid values in ComputeInternalEnergyFromDensityAndPressure: " + std::to_string(rho));
            return 0.0;
        }
        return p / (gamma - 1.0);
    }

    virtual bool ContainsNan() {
        Util::Message(INFO, "Checking for NaN values in Fluid: " + std::to_string(std::numeric_limits<double>::quiet_NaN()));
        return false;  // Default behavior assumes no NaNs unless overridden.
    }
};

}  // namespace Fluid
}  // namespace Model

#endif  // MODEL_FLUID_H_

