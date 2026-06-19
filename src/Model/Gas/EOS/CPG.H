#ifndef MODEL_GAS_EOS_CPG_H_
#define MODEL_GAS_EOS_CPG_H_

#include "IO/ParmParse.H"
#include "Model/Gas/EOS/EOS.H"
#include <memory>

namespace Model {
namespace Gas {
class Gas; // Forward declare Gas

namespace EOS {

class CPG  // : public EOS 
{
public:
    static constexpr const char* name = "cpg";
    //virtual const char* model_name() const override { return name; }

private:
    //Gas const * const gas;

public:
    // CPG() = default;
    // CPG(IO::ParmParse& pp, std::string name)
    // {
    //     pp.queryclass(name, *this);
    // }

    //~CPG() override = default;

    static void Parse(CPG & /*value*/, IO::ParmParse & /*pp*/) {}

public:
    AMREX_GPU_HOST_DEVICE
    double ComputeT(double density, double momentumx, double momentumy, double E, double Tguess,
                    double R, double gamma, double /*rtol=1e-12*/) const;

    AMREX_GPU_HOST_DEVICE
    double ComputeT(double pressure, double density,
                    double R) const;

    AMREX_GPU_HOST_DEVICE
    double ComputeP(double density, double T,
                    double R) const;

    AMREX_GPU_HOST_DEVICE
    double ComputeE(double density, double momentumx, double momentumy, double T,
                    double R, double gamma) const;

}; // class CPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif
