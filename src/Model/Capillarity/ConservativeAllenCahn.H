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
// Volume-constrained Allen--Cahn gradient flow of the configured multiphase
// free energy.  For every mobile pair (i,j), LowMach advances
//
//   d eta_i/dt = -L w_ij [mu_i-mu_j-lambda_ij],
//   d eta_j/dt = -d eta_i/dt,
//
// where w_ij=4 eta_i eta_j and lambda_ij is its weighted composite-domain
// average.  The pair exchange conserves both phase volumes, preserves the
// Gibbs sum locally, and dissipates free energy:
//
//   dF/dt = -L integral w_ij
//                    (mu_i-mu_j-lambda_ij)^2 dV <= 0.
//
// LowMach treats the stiff local and gradient curvatures implicitly and uses
// a free-energy line search for the nonlinear remainder.  Unlike the former
// profile-compression equation, every force driving this model is a variation
// of the same free energy used by the capillary stress.
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

    Set::Scalar Mobility(
        const Set::Scalar interface_thickness,
        const Set::Scalar reference_surface_energy) const
    {
        Util::Assert(INFO, TEST(interface_thickness > 0.0));
        Util::Assert(INFO, TEST(reference_surface_energy > 0.0));
        return interface_thickness /
            (reference_surface_energy * relaxation_time);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar PairMobilityWeight(
        const Set::Scalar eta, const Set::Scalar other)
    {
        return 4.0 * eta * other;
    }

private:
    Set::Scalar relaxation_time = 0.0;
};
}
}

#endif
