#ifndef MODEL_GAS_EOS_TPG_H_
#define MODEL_GAS_EOS_TPG_H_

#include <cmath>
#include <memory>
#include "IO/ParmParse.H"
#include "Model/Gas/Gas.H"
#include "Model/Gas/EOS/EOS.H"
#include "Util/Util.H"

namespace Model {
namespace Gas {
class Gas; // Forward declare Gas

namespace EOS {

class TPG : public EOS {
public:
    static constexpr const char* name = "tpg";

private:
    Gas const * const gas;

public:
    TPG() = delete;
    TPG(Gas *a_gas, IO::ParmParse& pp, std::string name) : gas(a_gas)
    {
        pp.queryclass(name, *this);
    }

    ~TPG() override = default;

    static void Parse(TPG & /*value*/, IO::ParmParse & /*pp*/) {}

public:
    double ComputeT(double density, double momentumx, double momentumy, double E, double Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const override;
    double ComputeT(double pressure, double density,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override;
    double ComputeP(double density, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override;
    double ComputeE(double density, double momentumx, double momentumy, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override;

}; // class TPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif
