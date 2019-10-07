
///
/// \file StressRelaxation.H
///
#ifndef STRESSRELAXATION_H
#define STRESSRELAXATION_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include <AMReX_MLMG.H>

#include "Integrator/Integrator.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "IC/IC.H"
#include "IC/Constant.H"
#include "Numeric/Interpolator/Linear.H"
#include "Util/Util.H"
#include "IO/FileNameParse.H"


#include "Operator/Elastic.H"
#include "Model/Solid/Viscoelastic/Isotropic.H"
#include "Model/Solid/Viscoelastic/Viscoelastic.H"

namespace Integrator
{

///
/// \class PolymerDegradation
/// \brief Phase field degradation model for hydro-thermal induced degradation of polymers.
///
/// Solve damage evolution laws for damage variable \f$ \eta \f$ or tensor \f$ \mathbf{F}_d \f$,
/// degrade material modulus tensor based on damage variable, and perform elasticity tests.
///
class StressRelaxation : public Integrator::Integrator
{
public:
	StressRelaxation();
	~StressRelaxation() {};

protected:

	/// \fn    Advance
	/// \brief Evolve phase field in time
	void Advance (int lev, Real time, Real dt);

	void Initialize (int lev);

	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

	void TimeStepBegin(amrex::Real time, int iter);

	void TimeStepComplete(amrex::Real time, int iter);

	//void DegradeMaterial(int lev);

private:

	// For mechanical degradation to be studied later
	amrex::Vector<amrex::MultiFab> displacement; ///< Multicomponent field variable storing pointwise displacements
	amrex::Vector<amrex::MultiFab> strain; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<amrex::MultiFab> stress; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<amrex::MultiFab> stress_vm; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<amrex::MultiFab> rhs; ///< Multicomponent field variable storing pointwise displacements

	amrex::Vector<amrex::MultiFab> energy; ///< Single component variable storing strain energy
	amrex::Vector<amrex::MultiFab> energies; ///< Multicomponent field variable storing individual strain energies

	int nlevels = 1;

	std::string ic_type, gb_type, filename;


	// Elasticity parameters
	bool        elastic_on 				= false;
	int         elastic_int 			= 1;
	std::string elastic_type;
	int         elastic_max_iter 		= 200;
	int         elastic_max_fmg_iter 	= 0;
	int 		elastic_bottom_max_iter	= 200;
	int         elastic_verbose 		= 3;
	int         elastic_cgverbose 		= 3;
	Set::Scalar elastic_tol_rel 		= 1.0E-5;
	Set::Scalar elastic_tol_abs 		= 1.0E-5;
	Set::Scalar elastic_tstart 			= 0.0;
	Set::Scalar elastic_tend 			= 0.0;
	std::string bottom_solver			= "bicgstab";
	int 		linop_maxorder 	     	= 2;
	bool	    elastic_use_fsmooth 	= false;
	int 		max_coarsening_level	= 0;
	bool agglomeration 	  = true;
	bool consolidation 	  = false;

	// Elastic BC
	std::map<std::string,Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC > 			bc_map;
	std::array<Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC,AMREX_SPACEDIM> 	bc_x_lo, bc_x_hi;
#if AMREX_SPACEDIM > 1
	std::array<Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC,AMREX_SPACEDIM> 	bc_y_lo, bc_y_hi;
#if AMREX_SPACEDIM > 2
	std::array<Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC,AMREX_SPACEDIM> 	bc_z_lo, bc_z_hi;
#endif
#endif

	amrex::Vector<amrex::BoxArray> ngrids;

	amrex::Vector<Set::Vector> elastic_bc_left;
	amrex::Vector<Set::Vector> elastic_bc_right;
	amrex::Vector<Set::Vector> elastic_bc_bottom;
	amrex::Vector<Set::Vector> elastic_bc_top;
	amrex::Vector<Set::Vector> elastic_bc_back;
	amrex::Vector<Set::Vector> elastic_bc_front;
	amrex::Vector<Set::Scalar> elastic_bc_left_t;
	amrex::Vector<Set::Scalar> elastic_bc_right_t;
	amrex::Vector<Set::Scalar> elastic_bc_bottom_t;
	amrex::Vector<Set::Scalar> elastic_bc_top_t;
	amrex::Vector<Set::Scalar> elastic_bc_back_t;
	amrex::Vector<Set::Scalar> elastic_bc_front_t;


	amrex::Vector<Set::Scalar> body_force = {AMREX_D_DECL(0.0,0.0,0.0)};
	Operator::Elastic<Model::Solid::Viscoelastic::Isotropic> elastic_operator;
	LPInfo info;

	// Material model
	amrex::Vector<amrex::FabArray<amrex::BaseFab<Model::Solid::Viscoelastic::Isotropic> > > model;
	Model::Solid::Viscoelastic::Isotropic modeltype;
	std::string 	input_material = "isotropic";

	std::vector<std::string> PlotFileNameNode (std::string plot_file_name, int lev) const;

	amrex::Vector<std::unique_ptr<amrex::MultiFab> >water_conc;
	IC::IC			*water_ic;
	BC::BC			*water_bc;

};
}
#endif
