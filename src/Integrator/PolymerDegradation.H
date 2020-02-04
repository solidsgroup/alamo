///
/// \file PolymerDegradation.H
///
#ifndef POLYMERDEGRADATION_H
#define POLYMERDEGRADATION_H

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
#include "BC/Operator/Elastic.H"
#include "IC/IC.H"
#include "IC/Constant.H"
#include "Numeric/Interpolator/Linear.H"
#include "Util/Util.H"
#include "IO/FileNameParse.H"


#include "Operator/Elastic.H"
#include "Model/Solid/LinearElastic/Degradable/Isotropic.H"
#include "Model/Solid/LinearElastic/LinearElastic.H"
//#include "Model/Solid/LinearElastic/Laplacian.H"


namespace Integrator
{
using pd_model_type = Model::Solid::LinearElastic::Degradable::Isotropic;
//using pd_model_type = Model::Solid::LinearElastic::Isotropic;
//using pd_model_type = Model::Solid::LinearElastic::Laplacian;

///
/// \class PolymerDegradation
/// \brief Phase field degradation model for hydro-thermal induced degradation of polymers.
///
/// Solve damage evolution laws for damage variable \f$ \eta \f$ or tensor \f$ \mathbf{F}_d \f$,
/// degrade material modulus tensor based on damage variable, and perform elasticity tests.
///
class PolymerDegradation : public Integrator::Integrator
{
public:
	PolymerDegradation();
	~PolymerDegradation();

protected:

	/// \fn    Advance
	/// \brief Evolve phase field in time
	void Advance (int lev, Real time, Real dt);

	void Initialize (int lev);

	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

	void TimeStepBegin(amrex::Real time, int iter);

	void TimeStepComplete(amrex::Real time, int iter);

	void DegradeMaterial(int lev,amrex::FabArray<amrex::BaseFab<pd_model_type> > &model);

private:

	int number_of_ghost_cells = 2;

	// Degradation variable
	amrex::Vector<std::unique_ptr<amrex::MultiFab> >eta_new; 		///< Degradation variable for the __current__timestep
	amrex::Vector<std::unique_ptr<amrex::MultiFab> >eta_old;		///< Degradation variable for the __previous__timestep

	amrex::Vector<std::unique_ptr<amrex::MultiFab> >eta_w_new; 		///< Degradation variable for the __current__timestep
	amrex::Vector<std::unique_ptr<amrex::MultiFab> >eta_w_old;		///< Degradation variable for the __previous__timestep

	amrex::Vector<std::unique_ptr<amrex::MultiFab> >eta_T_new; 		///< Degradation variable for the __current__timestep
	amrex::Vector<std::unique_ptr<amrex::MultiFab> >eta_T_old;		///< Degradation variable for the __previous__timestep

	// For water induced degradation
	amrex::Vector<std::unique_ptr<amrex::MultiFab> >water_conc;		///< Water concentration for the __current__timestep
	amrex::Vector<std::unique_ptr<amrex::MultiFab> >water_conc_old;	///< Water concentration for the __previous__timestep

	// For temperature induced degradation
	amrex::Vector<std::unique_ptr<amrex::MultiFab> >Temp;			///< Temperature for the __current__timestep
	amrex::Vector<std::unique_ptr<amrex::MultiFab> >Temp_old;		///< Temperature for the __previous__timestep

	// For mechanical degradation to be studied later
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > displacement; ///< Multicomponent field variable storing pointwise displacements
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > strain; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress_vm; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > rhs; ///< Multicomponent field variable storing pointwise displacements
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > residual;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > energy; ///< Single component variable storing strain energy
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > energies; ///< Multicomponent field variable storing individual strain energies

	int nlevels;

	std::string ic_type, gb_type, filename;


	amrex::Vector<amrex::BoxArray> ngrids,cgrids;
	amrex::Vector<amrex::DistributionMapping> ndmap;

	// Water diffusion parameters
	struct{
		bool 			on 						=	false;
		Set::Scalar 	diffusivity				=	1.0;
		Set::Scalar 	refinement_threshold 	=	0.01;
		std::string 	ic_type;
		IC::IC			*ic;
		BC::BC			*bc;
	} water;

	// Thermal diffusion parameters
	struct{
		bool			on 						=	false;
		Set::Scalar 	diffusivity 			=	1.0;
		Set::Scalar 	refinement_threshold 	=	0.01;
		std::string		ic_type;
		IC::IC			*ic;
		BC::BC			*bc;
	} thermal;

	// Damage parameters
	struct{
		std::string		type;
		bool			anisotropy 					= false;
		int				number_of_eta 				= 1;
		Set::Scalar		refinement_threshold		= 0.01;
		std::string 	ic_type;
		IC::IC			*ic;
		BC::BC			*bc;
	} damage;

	// Damage model: water relaxation
	struct{
		int 						number_of_terms	= 4;
		amrex::Vector<Set::Scalar> 	d_i;
		Set::Scalar 				d_final 		= 1.0;
		amrex::Vector<Set::Scalar> 	tau_i;
		amrex::Vector<Set::Scalar> 	t_start_i;
	} damage_w;

	// Damage model: temperature
	struct{
		Set::Scalar c0;
		Set::Scalar c1;
		Set::Scalar c2;
		Set::Scalar c3;
		Set::Scalar tau_T;
	} damage_T;

	// Elasticity parameters
	struct{
		Set::Scalar bottom_tol 				= 1.0E-6;
		bool        on 						= false;
		int         interval	 			= 1;
		std::string type;
		int         max_iter 				= 200;
		int         max_fmg_iter 			= 0;
		int 		bottom_max_iter			= 200;
		int 		max_fixed_iter 			= 5;
		int         verbose 				= 3;
		int         cgverbose 				= 3;
		Set::Scalar tol_rel 				= 1.0E-5;
		Set::Scalar tol_abs 				= 1.0E-5;
		Set::Scalar cg_tol_rel 				= 1.0E-5;
		Set::Scalar cg_tol_abs 				= 1.0E-5;
		Set::Scalar tstart 					= 0.0;
		Set::Scalar tend 					= 0.0;
		std::string bottom_solver			= "bicgstab";
		int 		linop_maxorder 	     	= 2;
		bool	    use_fsmooth 			= false;
		int 		max_coarsening_level	= 0;
		bool 		agglomeration 	  		= true;
		bool 		consolidation 	  		= false;

		// Elastic BC
		std::array<BC::Operator::Elastic<pd_model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xlo, bc_ylo, bc_zlo);
		std::array<BC::Operator::Elastic<pd_model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xhi, bc_yhi, bc_zhi);

		Set::Vector AMREX_D_DECL(bc_left,bc_bottom,bc_back);
		Set::Vector AMREX_D_DECL(bc_right,bc_top,bc_front);
		//amrex::Vector<Set::Scalar> AMREX_D_DECL(bc_left_t,bc_bottom_t,bc_back_t);
		//amrex::Vector<Set::Scalar> AMREX_D_DECL(bc_right_t,bc_top_t,bc_front_t);
		
		amrex::Vector<Set::Scalar> body_force = {AMREX_D_DECL(0.0,0.0,0.0)};
	} elastic;

	pd_model_type 		*modeltype;
	std::string 	input_material = "isotropic";

	std::vector<std::string> PlotFileNameNode (std::string plot_file_name, int lev) const;

	Numeric::Interpolator::Linear<Set::Vector> AMREX_D_DECL(interpolate_left,interpolate_bottom,interpolate_back);
	Numeric::Interpolator::Linear<Set::Vector> AMREX_D_DECL(interpolate_right,interpolate_top,interpolate_front);

};
}
#endif
