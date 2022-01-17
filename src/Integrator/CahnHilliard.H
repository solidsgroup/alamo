///
/// \file CahnHilliard.H
///
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
#include "Operator/Implicit/Implicit.H"

namespace Integrator
{

class CahnHilliard : public Integrator
{
public:
    CahnHilliard();

protected:

    void Initialize (int lev) override;
    void TimeStepBegin(amrex::Real /*time*/, int /*iter*/) override;
    void Advance (int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

private:

    Set::Field<Set::Scalar> etanewmf; 
    Set::Field<Set::Scalar> etaoldmf; 
    Set::Field<Set::Scalar> intermediate; 

    const int nghost = 1;
    const int ncomp = 1;
    BC::BC<Set::Scalar> *bc;
    IC::IC *ic;
    
    const Set::Scalar gamma = 0.0005;

    Operator::Implicit::Implicit op;
};
}
#endif
