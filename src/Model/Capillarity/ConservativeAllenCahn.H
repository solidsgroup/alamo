#ifndef MODEL_CAPILLARITY_CONSERVATIVEALLENCAHN_H
#define MODEL_CAPILLARITY_CONSERVATIVEALLENCAHN_H

#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Set/Set.H"
#include "Unit/Unit.H"
#include "Util/Util.H"

namespace Model
{
namespace Capillarity
{
// Locally conservative second-order profile restoration
//
//   d eta/dt = div[D (grad eta - 4 eta(1-eta)n/ell)],
//   D = ell^2/tau.
//
// The diffusion is treated implicitly by LowMach and the bounded compression
// flux is lagged.  This restores the declared tanh profile without changing
// the integral of a liquid phase.
class ConservativeAllenCahn
{
public:
    static constexpr const char* name = "conservative_allen_cahn";

    static void Parse(ConservativeAllenCahn& value, IO::ParmParse& pp)
    {
        pp.query_required("relaxation_time", value.relaxation_time,
                          Unit::Time());
        if (!(value.relaxation_time > 0.0))
            Util::Exception(INFO, "relaxation_time must be positive");
    }

    Set::Scalar Diffusivity(const Set::Scalar interface_thickness) const
    {
        return interface_thickness * interface_thickness / relaxation_time;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar CompressionMagnitude(
        const Set::Scalar eta, const Set::Scalar interface_thickness)
    {
        return 4.0 * eta * (1.0 - eta) / interface_thickness;
    }

private:
    Set::Scalar relaxation_time = 0.0;
};
}
}

#endif
