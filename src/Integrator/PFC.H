//
// Basic implementation of the phase field crystal model of Elder et al, 2002.
//
// Free energy functional is:
//
// .. math::
//
//    \mathcal{F} = \int \Big[ \frac{1}{2}\eta\Big((q^2+\nabla^2)^2 - \epsilon\Big)\eta + \frac{1}{4}\eta^4 \Big]d\mathbf{x}
//
// Order parameter evolves with Cahn-Hilliard kinetics:
//
// .. math::
//
//    \frac{\partial \eta}{\partial t} = \nabla^2 \frac{\delta\mathcal{F}}{\delta\eta}
//
// The variational derivative is
//
// .. math::
//
//      \frac{\delta\mathcal{F}}{\delta\eta} = \eta^3 + (q^4-\epsilon)\eta + 2q^2\nabla^2\eta + \nabla^4\eta
//
// The semi-implicit spectral update is
//
// .. math::
//
//    \eta_{n+1} = \frac{\hat{\eta}_n - dt\,\mathbf{\omega}^2 \mathcal{F}[\eta^3_n]}
//                      {1 + dt\,[(q^4-\epsilon)\mathbf{\omega}^2 - 2q^2\mathbf{\omega}^4 + \mathbf{\omega}^6]}
//

#ifndef INTEGRATOR_PFC_H
#define INTEGRATOR_PFC_H

#include <AMReX.H>
#include <AMReX_MLMG.H>

#include "IC/IC.H"
#include "BC/BC.H"
#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"

namespace Integrator
{

class PFC : public Integrator
{
public:

    /// Basic constructor (don't use)
    PFC();
    
    /// Destroy pointers defined in Parse
    ~PFC();

    /// Use this constructor
    PFC(IO::ParmParse& pp) : PFC()
    { Parse(*this, pp); }

    /// Scan input values and initialize fields
    static void Parse(PFC& value, IO::ParmParse& pp);

protected:

    /// Set values in fields
    void Initialize (int lev) override;
    /// Integrate eta over one timestep on lev
    void Advance (int lev, Set::Scalar time, Set::Scalar dt) override;
    /// Mark any cells that need to be refined
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

private:

    Set::Field<Set::Scalar> eta_mf;           /// Order parameter field
    Set::Field<Set::Scalar> grad_chempot_mf;  /// Field to calculate FFT of nonlinar part 

    BC::BC<Set::Scalar> *bc; /// eta's bc object
    IC::IC<Set::Scalar> *ic; /// eta's ic object
    
    Set::Scalar q0 = NAN; 
    Set::Scalar eps  = NAN; 

};
}
#endif
