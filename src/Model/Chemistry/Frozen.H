//
// :code:Model::Chemistry::Frozen:code: fixes chemical production rates and heat
// release due to chemical reactions to 0.0.  Valid for flows where momentum time
// scales are sufficiently fast that chemical reactions are negligle (Da << 1).
//

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
    std::pair<std::array<Set::Scalar,NSPECIES>,Set::Scalar>
    ComputeChemistrySources(
        const Set::Scalar /*P*/,
        const Set::Scalar /*T*/,
        const std::array<Set::Scalar, NSPECIES>& /*rhoY*/,
        const Set::Scalar /*dt*/,
        const Model::Gas::Gas* /*gas*/) const
    {
        // No chemistry: species are frozen
        std::array<Set::Scalar, NSPECIES> wdot;
        for (int n = 0; n < NSPECIES; ++n)
        {
            wdot[n] = 0.0;
        }

        return {wdot, 0.0};
    }

    static void Parse(Frozen& /*value*/, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());
    }

};

} // namespace Chemistry
} // namespace Model

#endif
