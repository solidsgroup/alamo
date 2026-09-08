#ifndef MODEL_CAPILLARITY_SINGLYDEGENERATECAHNHILLIARD_H
#define MODEL_CAPILLARITY_SINGLYDEGENERATECAHNHILLIARD_H

#include <string>

#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Set/Set.H"
#include "Unit/Unit.H"
#include "Util/Util.H"

namespace Model
{
namespace Capillarity
{
// Variational, locally conservative Cahn--Hilliard kinetics
//
//   d eta/dt = div[M0 g(eta) grad(mu)],
//   g(eta) = 4 eta(1-eta),
//   M0 = ell^3/(sigma_ref tau).
//
// The single zeros of g suppress diffusion in pure phases while retaining
// finite interfacial mobility.  LowMach treats the fourth-order chemical-
// potential stiffness with a linearly implicit two-factor Helmholtz solve.
// By default, mechanical capillarity remains in the momentum equation and
// retains its resolved capillary-wave timestep restriction.  The optional
// overdamped coupling uses this gradient flow alone to relax morphology.  It
// conserves liquid volume but deliberately omits capillary acceleration and
// the associated mechanical timestep restriction.
class SinglyDegenerateCahnHilliard
{
public:
    static constexpr const char* name =
        "singly_degenerate_cahn_hilliard";

    static void Parse(SinglyDegenerateCahnHilliard& value,
                      IO::ParmParse& pp)
    {
        pp.query_required("relaxation_time", value.relaxation_time,
                          Unit::Time());
        if (!(value.relaxation_time > 0.0))
            Util::Exception(INFO, "relaxation_time must be positive");

        std::string momentum_coupling;
        pp.query_default("momentum_coupling", momentum_coupling,
                         "korteweg");
        if (momentum_coupling == "korteweg")
            value.couple_capillary_momentum = true;
        else if (momentum_coupling == "overdamped")
            value.couple_capillary_momentum = false;
        else
            Util::Exception(INFO, "momentum_coupling must be korteweg or "
                "overdamped");
    }

    bool CouplesCapillaryMomentum() const
    {
        return couple_capillary_momentum;
    }

    Set::Scalar ReferenceMobility(
        const Set::Scalar interface_thickness,
        const Set::Scalar reference_surface_energy) const
    {
        return interface_thickness * interface_thickness *
            interface_thickness /
            (reference_surface_energy * relaxation_time);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar PairMobilityWeight(
        const Set::Scalar eta, const Set::Scalar other_eta)
    {
        return 4.0 * eta * other_eta;
    }

private:
    Set::Scalar relaxation_time = 0.0;
    bool couple_capillary_momentum = true;
};
}
}

#endif
