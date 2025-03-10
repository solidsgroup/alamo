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
// which is integrated using a forward Euler scheme.
//
// This is tested in :ref:`CahnHilliard`

#ifndef INTEGRATOR_CAHNHILLIARD_H
#define INTEGRATOR_CAHNHILLIARD_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include <AMReX_MLMG.H>

#include "IC/Random.H"
#include "Integrator/Integrator.H"
#include "BC/Nothing.H"

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

};
}
#endif
