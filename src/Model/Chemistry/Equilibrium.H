#ifndef MODEL_CHEMISTRY_EQUILIBRIUM_H
#define MODEL_CHEMISTRY_EQUILIBRIUM_H

#include "AMReX_Extension.H"
#include "Model/Chemistry/Chemistry.H"
#include "IO/YamlParse.H"

using VecD = std::vector<double>; 
using MatD = std::vector<VecD>;

namespace Model
{
namespace Chemistry
{

struct Equilibrium
{
    static constexpr const char* name = "equilibrium";

    Equilibrium() = default; 
    Equilibrium(IO::ParmParse &pp, std::string name) 
    { 
        pp_queryclass(name,*this); 
    } 

    std::string yaml;
    CanteraYAML::CanteraKinetics kinetics;

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    void compute(
        const Set::Scalar /*P*/,
        const Set::Scalar /*T*/,
        const std::array<double, NSPECIES>& /*rhoY*/,
        std::array<double, NSPECIES>& wdot,
        Set::Scalar& qdot)
    {
        // Equilibrium Chemistry

        for (int n=0; n<NSPECIES; ++n) {
            wdot[n] = 0.0;
        }

        qdot = 0.0;
    }

    static void Parse(Equilibrium& value, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());
        pp_query_required("yaml", value.yaml);
        value.kinetics = CanteraYAML::ParseYaml(value.yaml);
        Util::Message(INFO, value.yaml);
        for (int n=0; n<NSPECIES; ++n) Util::Message(INFO, value.kinetics.species[n].name);
        
    }

};

} // namespace Chemistry
} // namespace Model

#endif
