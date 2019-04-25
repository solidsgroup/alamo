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

namespace Integrator
{

class CahnHilliard : public Integrator
{
public:
	CahnHilliard();

protected:

	void Advance (int lev, Set::Scalar time, Set::Scalar dt) override;
	void Initialize (int lev) override;
	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

private:

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > etanewmf; 
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > etaoldmf; 
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > intermediate; 

	const int nghost = 1;
	const int ncomp = 1;
	BC::BC *bc;
	IC::IC *ic;
	
	const Set::Scalar gamma = 0.0005;
};
}
#endif
