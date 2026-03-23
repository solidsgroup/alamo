#ifndef MODEL_GAS_EOS_CPG_H_
#define MODEL_GAS_EOS_CPG_H_

#include "IO/ParmParse.H"
#include "Model/Gas/EOS/EOS.H"
#include <memory>

namespace Model {
namespace Gas {
class Gas; // Forward declare Gas

namespace EOS {

class CPG : public EOS {
public:
    static constexpr const char* name = "cpg";
    virtual const char* model_name() const override { return name; }

private:
    Gas const * const gas;

public:
    CPG() = delete;
    CPG(Gas *a_gas, IO::ParmParse& pp, std::string name) : gas(a_gas)
    {
        pp.queryclass(name, *this);
    }

    ~CPG() override = default;

    static void Parse(CPG & /*value*/, IO::ParmParse & /*pp*/) {}

public:
    double ComputeT(double density, double momentumx, double momentumy, double E, double Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, double /*rtol=1e-12*/) const override;
    double ComputeT(double pressure, double density,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override;
    double ComputeP(double density, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override;
    double ComputeE(double density, double momentumx, double momentumy, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override;

}; // class CPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif
