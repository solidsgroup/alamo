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
// Both SDCH methods use conservative face fluxes and the same discrete Laplacian.
// The spectral method applies the Fourier multiplier :math:`1/(1+\Delta t S D^2)`
// to the increment, where :math:`D` is the negative discrete Laplacian symbol
// and :math:`S` is :code:`spectral_stabilization`, default
// :math:`\gamma(L/16+M_0)`. SDCH requires periodic boundaries. Its spectral
// method assembles a uniform finest-level grid and restricts the update to
// the AMR hierarchy; realspace AMR has no flux refluxing at coarse/fine interfaces.
// By default the solution is not clipped, so timesteps must resolve the
// interface dynamics. Clipping can change the conserved mean.
// Optional input and field names are :code:`input_name`, :code:`field_name`,
// :code:`old_field_name`, and :code:`intermediate_name`. SDCH can optionally
// clip its updated field with :code:`clipping` (default false). Derived integrators can override
// :code:`FillMobilityMask` to apply a spatial multiplier to the face mobility.
// The hook receives each level's native box layout and valid cells; the caller
// fills periodic ghosts. Set :code:`diagnostics = true` for SDCH mobility, mask, flux, RHS, raw update,
// and update delta fields alongside the solution.
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

class CahnHilliard : virtual public Integrator
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

    /// Fill a cell-centered multiplier for the singly-degenerate mobility.
    /// The default is one everywhere; derived integrators may override it.
    virtual void FillMobilityMask(int lev, amrex::MultiFab& mask);
    void InitializeDiagnostics(int lev);

protected:

    Set::Field<Set::Scalar> etanew_mf;     /// The new value for eta this timestep
    Set::Field<Set::Scalar> etaold_mf;     /// Last timestep's value for eta
    Set::Field<Set::Scalar> intermediate;  /// Intermediate field used for CH kinetics
    Set::Field<Set::Scalar> diagnostic_mobility;
    Set::Field<Set::Scalar> diagnostic_mobility_mask;
    Set::Field<Set::Scalar> diagnostic_flux;
    Set::Field<Set::Scalar> diagnostic_rhs;
    Set::Field<Set::Scalar> diagnostic_raw;
    Set::Field<Set::Scalar> diagnostic_delta;

    BC::BC<Set::Scalar> *bc = nullptr; /// eta's bc object
    IC::IC<Set::Scalar> *ic = nullptr; /// eta's ic object
    
    Set::Scalar gamma = NAN; 
    Set::Scalar L     = NAN; 
    Set::Scalar refinement_threshold = NAN; 

    std::string method;
    std::string mobility = "constant";
    Set::Scalar mobility_floor = 0.0;
    Set::Scalar spectral_stabilization = NAN;
    Set::Scalar tstart = 1E100;
    Set::Scalar tstart_timestep = -1.0;
    bool clipping = false;
    bool output_diagnostics = false;
    std::string input_name = "eta";
    std::string field_name = "eta";
    std::string old_field_name = "eta_old";
    std::string intermediate_name = "int";
};
}
#endif
