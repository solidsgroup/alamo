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

    Set::Field<Set::Scalar> eta_mf;     /// The new value for eta this timestep
    Set::Field<Set::Scalar> grad_chempot_mf;  /// Intermediate field used for CH kinetics

    BC::BC<Set::Scalar> *bc; /// eta's bc object
    IC::IC<Set::Scalar> *ic; /// eta's ic object
    
    Set::Scalar q0 = NAN; 
    Set::Scalar r  = NAN; 
    Set::Scalar refinement_threshold = NAN; 

};
}
#endif
