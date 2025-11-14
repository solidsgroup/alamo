#ifndef MODEL_GAS_THERMO_H_
#define MODEL_GAS_THERMO_H_

#include <vector>
#include <memory>
#include <cmath>
#include "Model/Gas/Gas.H"

namespace Model {
namespace Gas {
namespace Thermo {

class Thermo {
protected:
    std::shared_ptr<GasData> data;

public:
    Thermo(std::shared_ptr<GasData> data_in) : data(data_in) {}
    virtual ~Thermo() = default;

    // Specific heat (constant pressure), J/(kmol-K)
    virtual double cp_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;
    // Specific enthalpy, J/kmol)
    virtual double enthalpy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;
    // specific entropy, J/(kmol-K)
    virtual double entropy_mol(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const = 0;

}; // class Thermo

} // namespace Thermo
} // namespace Gas
} // namespace Model

#endif
