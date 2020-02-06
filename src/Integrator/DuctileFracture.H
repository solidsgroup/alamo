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
#include "Solver/Nonlocal/Newton.H"
#include "Model/Solid/Affine/IsotropicDegradablePlastic.H"
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
using ductile_fracture_model_type = Model::Solid::Affine::IsotropicDegradablePlastic;
class DuctileFracture : public Integrator
{
public:

	/// \brief Read in parameters and register field variables
	DuctileFracture();
	~DuctileFracture();

protected:

	/// \brief Use the #ic object to initialize #Temp
	void Initialize (int lev) override;

	/// \brief This primarily solves the inhomogeneous elasticity problem.
	void TimeStepBegin(amrex::Real /*time*/, int /*iter*/) override;

	/// \brief This function checks for whether the solution has converged or not.
	void TimeStepComplete(amrex::Real time,int iter) override;

	/// \brief This function advances the crack field, plastic variable fields.
	void Advance (int lev, amrex::Real /*time*/, amrex::Real dt) override;

	/// \brief Tag cells for mesh refinement based on crack field gradient
	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;

	/// \brief Scales the elastic modulus in the domain based on crack and plasticity field variables.
	void ScaledModulus(int lev,amrex::FabArray<amrex::BaseFab<ductile_fracture_model_type> > &model);

	/// \brief Computes the error norm at every time step
	void Integrate(int amrlev, Set::Scalar time, int step,const amrex::MFIter &mfi, const amrex::Box &box) override;

private:
	int number_of_ghost_cells = 2;											///< Number of ghost cells

	Set::Field<Set::Scalar> m_c;					///< Crack field fab
	Set::Field<Set::Scalar> m_c_old;				///< Crack field fab at previous time step
	Set::Field<Set::Scalar> m_driving_force;		///< Crack driving force
	
	Set::Field<Set::Scalar> m_disp;				///< Displacement variable
	Set::Field<Set::Scalar> m_strain;				///< Total strain variable
	Set::Field<Set::Scalar> m_strain_p;			///< Plastic strain variable
	Set::Field<Set::Scalar> m_p;					///< Plastic exponent for crack function g(c)
	Set::Field<Set::Scalar> m_lambda;				///< Evolving hardening variable

	Set::Field<Set::Scalar> m_strain_pold;			///< Plastic strain at previous time step
	Set::Field<Set::Scalar> m_pold;				///< Plastic exponent at previous time step
	Set::Field<Set::Scalar> m_lambdaold;			///< Hardening variable at previous time step

	Set::Field<Set::Scalar> m_stress;				///< Stress field
	Set::Field<Set::Scalar> m_stressdev;			///< Deviatoric Stress field
	Set::Field<Set::Scalar> m_rhs;					///< RHS field for inhomogeneous problem
	Set::Field<Set::Scalar> m_residual;			///< Implicit problem residual field
	Set::Field<Set::Scalar> m_energy;				///< Elastic energy field
	Set::Field<Set::Scalar> m_energy_pristine;		///< Elastic energy field as if crack doesn't exist.

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
		//std::array<BC::Operator::Elastic<ductile_fracture_model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xlo, bc_ylo, bc_zlo);
		//std::array<BC::Operator::Elastic<ductile_fracture_model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xhi, bc_yhi, bc_zhi);

		//Set::Vector AMREX_D_DECL(bc_left,bc_bottom,bc_back);
		//Set::Vector AMREX_D_DECL(bc_right,bc_top,bc_front);
		
		Set::Vector body_force = Set::Vector::Zero();
		Operator::Elastic<ductile_fracture_model_type> op;
		BC::Operator::Elastic<ductile_fracture_model_type> bc;
		Set::Scalar bc_top 							= 0.;
	} elastic;

	struct{
		Model::Interface::Crack::Crack *boundary;
		std::string ic_type;
		IC::IC *ic;            						///< Pointer to abstract IC object
		BC::BC *mybc;                   			///< Stock generic boundary condition object

		bool newCrackProblem=true;
		bool solveCrack=false;
		bool solveElasticity=false;
		bool crackStressTest = false;
	}crack;

	struct{
		ductile_fracture_model_type 		*modeltype;
		Set::Field<ductile_fracture_model_type> model;
		std::string input_material = "isotropic";
		Set::Scalar yield_strength = 1.e3;
		Set::Scalar hardening_modulus = 1.;
		Set::Scalar epscrit = 1.0;
	}material;

	Set::Scalar crack_err_norm, c_new_norm;
	bool err_crack_init = true;
	bool err_crack_temp_init = true;
	bool err_disp_init = true;

	Set::Scalar eta_epsilon, mobility, scaleModulusMax;
	Set::Scalar max_error;

	Set::Scalar refinement_threshold = 0.001;			///< Criterion for cell refinement
	Set::Scalar tol_crack = 1e-7, tol_step = 1e-3;
	int nlevels;
};
}
#endif
