#ifndef DUCTILEFRACTURE_INTEGRATOR_H
#define DUCTILEFRACTURE_INTEGRATOR_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/BC.H"

#include "IC/Ellipsoid.H"
#include "IC/Notch.H"

#include "BC/Operator/Elastic.H"
#include "BC/Constant.H"

#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "Model/Solid/LinearElastic/Degradable/Isotropic.H"
#include "Model/Solid/LinearElastic/LinearElastic.H"
#include "Model/Solid/LinearElastic/Laplacian.H"

#include "Model/Interface/Crack/Crack.H"
#include "Model/Interface/Crack/Constant.H"

#include "Numeric/Stencil.H"

/// \class Fracture::Integrator
/// \brief Solve the heat diffusion equation in 2D or 3D.
///
/// This is a nontrivial example class that demonstrates how to write an explicit AMR
/// code using the #Integrator virtual class that abstracts the AmrBase class
/// from Amrex.
///
/// For more details:
///    - See documentation on #Initialize for input parameters
///    - See documentation on #Advance for equations and discretization
///    - See documentation on #TagCellsForRefinement for cell refinement criteria
/// For boundary conditions:
///    - See #BC
/// For initial conditions:
///    - See #Cylinder or #Constant
namespace Integrator
{
using ductile_fracture_model_type = Model::Solid::LinearElastic::Degradable::Isotropic;
class DuctileFracture : public Integrator
{
public:

	/// \brief Read in parameters and register field variables
	DuctileFracture();
	~DuctileFracture();

protected:

	/// \brief Use the #ic object to initialize #Temp
	void Initialize (int lev) override;

	void TimeStepBegin(amrex::Real /*time*/, int /*iter*/) override;

	void TimeStepComplete(amrex::Real time,int iter) override;

	/// \brief Integrate the heat equation
	void Advance (int lev, amrex::Real /*time*/, amrex::Real dt) override;

	/// \brief Tag cells for mesh refinement based on temperature gradient
	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;

	void ScaledModulus(int lev,amrex::FabArray<amrex::BaseFab<ductile_fracture_model_type> > &model);

	void Integrate(int amrlev, Set::Scalar time, int step,const amrex::MFIter &mfi, const amrex::Box &box) override;

private:
	int number_of_ghost_cells = 2;				///< Number of ghost cells

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_c;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_c_old;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_driving_force;
	
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_disp;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_strain;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_strain_p;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_p;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_lambda;

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_strain_pold;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_pold;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_lambdaold;

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_stress;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_stressdev;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_rhs;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_residual;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_energy;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_energy_pristine;

	struct{
		Set::Scalar bottom_tol 				= 1.0E-6;
		int         interval	 			= 1;
		std::string type					= "single";
		amrex::Vector<Set::Scalar> test_time;
		Set::Scalar test_rate				= 1.; 
		//Set::Scalar test_duration			= 2.;
		Set::Scalar test_max				= 1.0;
		Set::Scalar test_init 				= 0.01;
		int 		test_step				= 0;
		int			current_test			= 0;
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
		std::array<BC::Operator::Elastic<ductile_fracture_model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xlo, bc_ylo, bc_zlo);
		std::array<BC::Operator::Elastic<ductile_fracture_model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xhi, bc_yhi, bc_zhi);

		Set::Vector AMREX_D_DECL(bc_left,bc_bottom,bc_back);
		Set::Vector AMREX_D_DECL(bc_right,bc_top,bc_front);
		
		amrex::Vector<Set::Scalar> body_force = {AMREX_D_DECL(0.0,0.0,0.0)};
	} elastic;

	ductile_fracture_model_type 		*modeltype;
	amrex::Vector<amrex::FabArray<amrex::BaseFab<ductile_fracture_model_type> > > model;

	std::string 	input_material = "isotropic";
	Set::Scalar yield_strength = 1.e3, hardening_modulus = 1., epscrit = 1.0;

	Model::Interface::Crack::Crack *boundary;

	Set::Scalar crack_err_norm, crack_err_norm_init, crack_err_temp_norm, crack_err_temp_norm_init, c_new_norm;
	Set::Scalar disp_err_norm, disp_err_norm_init;

	bool err_crack_init = true, err_crack_temp_init = true, err_disp_init = true;

	Set::Scalar eta_epsilon, mobility, scaleModulusMax;
	Set::Scalar max_error;

	Set::Scalar refinement_threshold = 0.001;			///< Criterion for cell refinement
	Set::Scalar tol_crack = 1e-7, tol_step = 1e-3;
	int nlevels;

	std::string ic_type;						///< String to determine what kind of initial condition to use
	IC::IC *ic;            					///< Pointer to abstract IC object
	BC::BC *mybc;                   				///< Stock generic boundary condition object

	bool newCrackProblem=true, solveCrack=false, solveElasticity=false;

	bool crackStressTest = false;
	Set::Matrix identity = Set::Matrix::Zero();
};
}
#endif
