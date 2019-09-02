#ifndef FRACTURE_INTEGRATOR_H
#define FRACTURE_INTEGRATOR_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/BC.H"

#include "IC/HeatConduction/Cylinder.H"
#include "IC/HeatConduction/Constant.H"
#include "IC/Cuboid.H"
#include "IC/Sphere.H"
#include "IC/Ellipsoid.H"

#include "BC/Operator/Elastic.H"
#include "BC/Constant.H"

#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "Model/Solid/LinearElastic/Degradable/Isotropic.H"
#include "Model/Solid/LinearElastic/Degradable/Isotropic2.H"
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
using model_type = Model::Solid::LinearElastic::Degradable::Isotropic;
class Fracture : public Integrator
{
public:

	/// \brief Read in parameters and register field variables
	Fracture();
	~Fracture();

protected:

	/// \brief Use the #ic object to initialize #Temp
	void Initialize (int lev) override;

	void TimeStepBegin(amrex::Real /*time*/, int /*iter*/) override;

	void TimeStepComplete(amrex::Real /*time*/,int /*iter*/) override;

	/// \brief Integrate the heat equation
	void Advance (int lev, amrex::Real /*time*/, amrex::Real dt) override;

	/// \brief Tag cells for mesh refinement based on temperature gradient
	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;

	void ScaledModulus(int lev,amrex::FabArray<amrex::BaseFab<model_type> > &model);

	void Integrate(int amrlev, Set::Scalar time, int step,const amrex::MFIter &mfi, const amrex::Box &box) override;

	void CrackProblem(int lev, amrex::Real /*time*/, amrex::Real dt);
	void ElasticityProblem(amrex::Real /*time*/);

private:
	int number_of_ghost_cells = 2;				///< Number of ghost cells

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_c;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_c_old;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_c_conv;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_c_temp;

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_disp;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_disp_old;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_disp_conv;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_disp_temp;

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_strain;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_stress;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_stressvm;
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
		Set::Scalar test_dt					= 0.01;
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
		std::array<BC::Operator::Elastic<model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xlo, bc_ylo, bc_zlo);
		std::array<BC::Operator::Elastic<model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bc_xhi, bc_yhi, bc_zhi);

		Set::Vector AMREX_D_DECL(bc_left,bc_bottom,bc_back);
		Set::Vector AMREX_D_DECL(bc_right,bc_top,bc_front);
		
		amrex::Vector<Set::Scalar> body_force = {AMREX_D_DECL(0.0,0.0,0.0)};
	} elastic;

	model_type 		*modeltype;
	std::string 	input_material = "isotropic";
	Model::Interface::Crack::Crack *boundary;

	Set::Scalar crack_norm, crack_norm_old, crack_norm_conv, crack_norm_temp;
	Set::Scalar disp_norm, disp_norm_old, disp_norm_conv; 
	Set::Scalar eta_epsilon, mobility;
	Set::Scalar max_error;

	Set::Scalar refinement_threshold = 0.01;			///< Criterion for cell refinement
	int nlevels;

	std::string ic_type;						///< String to determine what kind of initial condition to use
	IC::IC *ic;            					///< Pointer to abstract IC object
	BC::BC *mybc;                   				///< Stock generic boundary condition object

	bool newCrackProblem=true, solveCrack=false, solveElasticity=false;
};
}
#endif
