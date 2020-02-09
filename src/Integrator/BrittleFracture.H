#ifndef BRITTLEFRACTURE_INTEGRATOR_H
#define BRITTLEFRACTURE_INTEGRATOR_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/BC.H"
#include "BC/Constant.H"

#include "IC/Ellipsoid.H"
#include "IC/Notch.H"

#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

#include "Model/Solid/Solid.H"
#include "Model/Solid/Linear/IsotropicDegradable.H"
//#include "Model/Solid/LinearElastic/LinearElastic.H"
//#include "Model/Solid/LinearElastic/Laplacian.H"

#include "Model/Interface/Crack/Crack.H"
#include "Model/Interface/Crack/Constant.H"
#include "Model/Interface/Crack/Sin.H"

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
using fracture_model_type = Model::Solid::Linear::IsotropicDegradable;
class BrittleFracture : public Integrator
{
public:

	/// \brief Read in parameters and register field variables
	BrittleFracture();
	~BrittleFracture();

protected:

	/// \brief Use the #ic object to initialize #Temp
	void Initialize (int lev) override;

	void TimeStepBegin(amrex::Real /*time*/, int /*iter*/) override;

	void TimeStepComplete(amrex::Real time,int iter) override;

	/// \brief Integrate the heat equation
	void Advance (int lev, amrex::Real /*time*/, amrex::Real dt) override;

	/// \brief Tag cells for mesh refinement based on temperature gradient
	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;

	void ScaledModulus(int lev,amrex::FabArray<amrex::BaseFab<fracture_model_type> > &model);

	void Integrate(int amrlev, Set::Scalar time, int step,const amrex::MFIter &mfi, const amrex::Box &box) override;

	void CrackProblem(int lev, amrex::Real /*time*/, amrex::Real dt);
	void ElasticityProblem(amrex::Real /*time*/);

private:
	int number_of_ghost_cells = 2;				///< Number of ghost cells

	Set::Field<Set::Scalar> m_c;
	Set::Field<Set::Scalar> m_c_old;
	Set::Field<Set::Scalar> m_driving_force;
	
	Set::Field<Set::Scalar> m_disp;
	Set::Field<Set::Scalar> m_strain;
	Set::Field<Set::Scalar> m_stress;
	Set::Field<Set::Scalar> m_stressvm;
	Set::Field<Set::Scalar> m_rhs;
	Set::Field<Set::Scalar> m_residual;
	Set::Field<Set::Scalar> m_energy;
	Set::Field<Set::Scalar> m_energy_pristine;
	Set::Field<Set::Scalar> m_energy_pristine_old;

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
		//std::array<BC::Operator::Elastic<fracture_model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xlo, bc_ylo, bc_zlo);
		//std::array<BC::Operator::Elastic<fracture_model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xhi, bc_yhi, bc_zhi);

		//Set::Vector AMREX_D_DECL(bc_left,bc_bottom,bc_back);
		//Set::Vector AMREX_D_DECL(bc_right,bc_top,bc_front);
		BC::Operator::Elastic<fracture_model_type> 	bc;
		Operator::Elastic<fracture_model_type> op;
		Set::Vector body_force 						= Set::Vector::Zero();
		Set::Scalar bc_top 							= 0.;

	} elastic;

	struct{
		fracture_model_type modeltype;
		Set::Field<fracture_model_type> model;
		std::string input_material = "isotropic";
	} material;

	struct
	{
		Model::Interface::Crack::Crack *boundary;
		std::string ic_type;
		IC::IC *ic;
		BC::BC *mybc;
		BC::BC *mybcdf;

		bool newCrackProblem	=true;
		bool solveCrack			=false; 
		bool solveElasticity	=false;

		bool crackStressTest 	=false;
	} crack;

	

	Set::Scalar crack_err_norm, crack_err_norm_init, crack_err_temp_norm, crack_err_temp_norm_init, c_new_norm;
	Set::Scalar disp_err_norm, disp_err_norm_init;

	bool err_crack_init = true, err_crack_temp_init = true, err_disp_init = true;

	Set::Scalar eta_epsilon, mobility, scaleModulusMax;
	Set::Scalar max_error;

	Set::Scalar refinement_threshold = 0.001;			///< Criterion for cell refinement
	Set::Scalar tol_crack = 1e-7, tol_step = 1e-3;
	int nlevels;

	
};
}
#endif
