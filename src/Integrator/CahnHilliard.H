//
// This implements a basic two-phase field model with Cahn-Hilliard kinetics.
//
// The free energy is
//
// .. math::
//
//    F[\eta] = \int_\Omega \Big[\frac{1}{4}(\eta^2 - 1)^2 +
//                               \frac{1}{2}\gamma |\nabla\eta|^2\Big] d\mathbf{x}
//
// The corresponding governing equation under conservative kinetics is
//
// .. math::
//
//    \frac{\partial\eta}{\partial t} = L\nabla^2\Big(\eta^3 - \eta - \gamma\nabla^2\eta\Big)
//
// Constant mobility is integrated using forward Euler in realspace or a
// semi-implicit spectral scheme. With :code:`mobility = singly_degenerate`,
// the equation is instead
//
// .. math::
//
//    \partial_t\eta = \nabla\cdot(M(\eta)\nabla\mu),\qquad
//    M(\eta) = L\alpha^2(1-\alpha)^2 + M_0,\qquad
//    \alpha = \mathrm{clip}((\eta+1)/2,0,1),
//
// where :math:`\mu=\eta^3-\eta-\gamma\nabla^2\eta` and
// :math:`M_0` is :code:`mobility_floor` (default zero).
// The realspace method uses conservative face fluxes; the spectral method
// stabilizes the flux divergence with :code:`spectral_stabilization`, default
// :math:`\gamma(L/16+M_0)`. SDCH requires periodic boundaries. Its spectral
// method assembles a uniform finest-level grid and restricts the update to
// the AMR hierarchy; realspace AMR has no flux refluxing at coarse/fine interfaces.
// The solution is not clipped, so timesteps must resolve the interface dynamics.
//
// This is tested in :ref:`CahnHilliard` and :ref:`SinglyDegenerateCahnHilliard`

#ifndef INTEGRATOR_CAHNHILLIARD_H
#define INTEGRATOR_CAHNHILLIARD_H

#include <AMReX.H>
#include <AMReX_MLMG.H>

#include "IC/IC.H"
#include "BC/BC.H"
#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"

namespace Integrator
{

class CahnHilliard : public Integrator
{
public:
    static constexpr const char* name = "cahnhilliard";

    /// Basic constructor (don't use)
    CahnHilliard();
    
    /// Destroy pointers defined in Parse
    ~CahnHilliard();

    /// Use this constructor
    CahnHilliard(IO::ParmParse& pp) : CahnHilliard()
    { Parse(*this, pp); }

    /// Scan input values and initialize fields
    static void Parse(CahnHilliard& value, IO::ParmParse& pp);

protected:

    /// Set values in fields
    void Initialize (int lev) override;
    /// Integrate eta over one timestep on lev
    void Advance (int lev, Set::Scalar time, Set::Scalar dt) override;
    /// spectral version of advance
    void AdvanceSpectral (int lev, Set::Scalar time, Set::Scalar dt);
    /// realspace version of advance
    void AdvanceReal (int lev, Set::Scalar time, Set::Scalar dt);
    /// Conservative realspace update with variable mobility
    void AdvanceDegenerateReal (int lev, Set::Scalar time, Set::Scalar dt);
    /// Stabilized spectral update with variable mobility
    void AdvanceDegenerateSpectral (int lev, Set::Scalar time, Set::Scalar dt);
    /// Mark any cells that need to be refined
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

private:

    Set::Field<Set::Scalar> etanew_mf;     /// The new value for eta this timestep
    Set::Field<Set::Scalar> etaold_mf;     /// Last timestep's value for eta
    Set::Field<Set::Scalar> intermediate;  /// Intermediate field used for CH kinetics

    BC::BC<Set::Scalar> *bc; /// eta's bc object
    IC::IC<Set::Scalar> *ic; /// eta's ic object
    
    Set::Scalar gamma = NAN; 
    Set::Scalar L     = NAN; 
    Set::Scalar refinement_threshold = NAN; 

    std::string method;
    std::string mobility = "constant";
    Set::Scalar mobility_floor = 0.0;
    Set::Scalar spectral_stabilization = NAN;
};
}
#endif
