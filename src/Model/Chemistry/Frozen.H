#ifndef MODEL_CHEMISTRY_FROZEN_H
#define MODEL_CHEMISTRY_FROZEN_H

#include "AMReX_Extension.H"
#include "Model/Chemistry/Chemistry.H"

namespace Model
{
namespace Chemistry
{

struct Frozen
{
    static constexpr const char* name = "frozen";

    Frozen() = default;
    Frozen(IO::ParmParse &pp, std::string name)
    {
        pp_queryclass(name,*this);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    void compute(
        const Set::Scalar /*P*/,
        const Set::Scalar /*T*/,
        const std::array<double, NSPECIES>& /*rhoY*/,
        std::array<double, NSPECIES>& wdot,
        Set::Scalar& qdot)
    {
        // No chemistry: species are frozen

        for (int n = 0; n < NSPECIES; ++n)
        {
            wdot[n] = 0.0;
        }

        qdot = 0.0;
    }

    static void Parse(Frozen& /*value*/, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());
    }

};

} // namespace Chemistry
} // namespace Model

#endif
