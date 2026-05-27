#ifndef MODEL_GAS_EOS_H_
#define MODEL_GAS_EOS_H_

#include <vector>
#include <memory>
#include <cmath>

#include "Set/Set.H"

namespace Model {
namespace Gas {
namespace EOS {

class EOS {
protected:

public:
    EOS() {} ;
    virtual ~EOS() = default;
    virtual const char* model_name() const = 0;

    // Temperature, K
    #if AMREX_SPACEDIM == 2
    virtual double ComputeT(
        double density, double momentumx, double momentumy, double E, double Tguess,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const = 0;
    // Energy, J/m^3
    virtual double ComputeE(
        double density, double momentumx, double momentumy, double T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;
    #elif AMREX_SPACEDIM == 3
    virtual double ComputeT(
        double density, double momentumx, double momentumy, double momentumz, double E, double Tguess,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const = 0;
    // Energy, J/m^3
    virtual double ComputeE(
        double density, double momentumx, double momentumy, double momentumz, double T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;
    #endif
    virtual double ComputeT_from_primitives(
        double pressure, double density,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;
    // Pressure, Pa
    virtual double ComputeP(
        double density, double T,
        Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;

}; // class EOS

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif
