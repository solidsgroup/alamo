#ifndef FRACTURE_INTEGRATOR_H
#define FRACTURE_INTEGRATOR_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Operator/Elastic/Constant.H"

#include "IC/IC.H"
#include "IC/Ellipsoid.H"
#include "IC/Notch.H"

#include "Operator/Operator.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

#include "Model/Solid/Solid.H"
#include "Model/Solid/Linear/IsotropicDegradable.H"
#include "Model/Solid/Affine/J2PlasticDegradable.H"
#include "Model/Solid/Affine/CrystalPlasticDegradable.H"

#include "Model/Interface/Crack/Crack.H"
#include "Model/Interface/Crack/Constant.H"
#include "Model/Interface/Crack/Sin.H"

#include "Numeric/Stencil.H"

/// \class Fracture::Integrator
/// \brief Solves ductile and brittle fracture in 2D. 3D has not been implemented yet
///
/// For more details:
///     - See documentation on #Initialize for input parameters
///     - See documentation on #Advance for equations and discretization
///     - See documentation on #TagCellsForRefinement for cell refinement criteria
/// For boundary conditions:
///     - See #BC
/// For initial conditions:
///    - See #Notch or #Ellipsoidal

namespace Integrator
{
using brittle_fracture_model_type = Model::Solid::Linear::IsotropicDegradable;
using ductile_fracture_model_type = Model::Solid::Affine::J2PlasticDegradable;
//using ductile_fracture_model_type = Model::Solid::Affine::CrystalPlasticDegradable;

class Fracture : public Integrator
{
public:

	/// \brief Read in parameters and register field variables
	Fracture();
	~Fracture();
	enum FractureType {Brittle, Ductile};
	enum ModeType {ModeI, ModeII, ModeIII};
	enum LoadType {Force, Displacement};

protected:

	/// \brief Use the #ic object to initialize field variables
	void Initialize (int lev) override;

    /// \brief Contains methods to solve for elastic fields
	void TimeStepBegin(amrex::Real /*time*/, int /*iter*/) override;

	void TimeStepComplete(amrex::Real time,int iter) override;

	/// \brief Time marching to advance crack related fields
	void Advance (int lev, amrex::Real /*time*/, amrex::Real dt) override;

	/// \brief Tag cells for mesh refinement based on crack field gradient
	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;

    /// \brief Perform integration of field variables for error norm calculations
	void Integrate(int amrlev, Set::Scalar time, int step,const amrex::MFIter &mfi, const amrex::Box &box) override;

private:
    int number_of_ghost_cells = 3;				///< Number of ghost cells
	int number_of_ghost_nodes = 2;				///< Number of ghost nodes
    int nlevels;

    FractureType fracture_type = FractureType::Brittle;   ///< Specifies the type of problem we are solving
	
	struct{
        Set::Field<Set::Scalar> field;                ///< crack field at current time step
        Set::Field<Set::Scalar> field_old;            ///< crack field at previous time step
        Set::Field<Set::Scalar> driving_force;         ///< crack driving forces.

        Model::Interface::Crack::Crack *cracktype;       ///< type of crack. See Crack/Constant or Crack/Sin
		std::string ic_type;                            ///< crack IC type. See IC/Notch and IC/Ellipsoid
		IC::IC *ic;                                     ///< crack IC. See IC/Notch and IC/Ellipsoid
		BC::BC<Set::Scalar> *bc;                      ///< BC objects for crack fields
		BC::BC<Set::Scalar> *bcdf;                     ///< BC object for crack driving force

        Set::Scalar mobility;                       ///< crack mobility
        Set::Scalar scaleModulusMax = 0.02;         ///< material modulus ratio inside crack (default = 0.02).
        Set::Scalar refinement_threshold = 0.001;

		Set::Scalar norm;							///< L2 norm of current crack field
		Set::Scalar error_norm;						///< L2 norm of different between current and old field
		Set::Scalar tol_abs = 1.e-6;
		Set::Scalar tol_rel = 1.e-6;

		Set::Scalar mult_df_Gc = 1.0;				///< Multiplier for Gc/zeta term
		Set::Scalar mult_df_lap = 1.0;				///< Multiplier for the laplacian term
    } crack;

    struct
	{
		bool on = false;
		Set::Scalar tstart = 0.0;
		Set::Scalar timestep = 0.0;
		Set::Scalar beta = 1.0e-2;
		int plot_int = 1;
		int elastic_int = 1;
		Set::Scalar plot_dt = 0.0;
	} anisotropy;
	
    struct{
        Set::Field<Set::Scalar> disp;             ///< displacement field
        Set::Field<Set::Scalar> strain;           ///< total strain field (gradient of displacement)
        Set::Field<Set::Scalar> stress;           ///< stress field
        Set::Field<Set::Scalar> stressvm;         ///< von-Mises stress
        Set::Field<Set::Scalar> rhs;              ///< rhs fab for elastic solution
        Set::Field<Set::Scalar> residual;         ///< residual field for solver
        Set::Field<Set::Scalar> energy;           ///< total elastic energy
        Set::Field<Set::Scalar> energy_pristine;      ///< energy of the prisitne material as if no crack is present
        Set::Field<Set::Scalar> energy_pristine_old;  ///< energy of the pristine material for previous time step.

        BC::Operator::Elastic::Constant 	brittlebc;  ///< elastic BC if using brittle fracture
        BC::Operator::Elastic::Constant 	ductilebc;  ///< elastic BC if using ductile fracture
		Set::Scalar df_mult = 1.0;              ///< mulitplier for elastic driving force.
    } elastic;

    struct{
        Set::Vector body_force 				= Set::Vector::Zero();
        ModeType 	mode                    = ModeType::ModeI;
        LoadType 	load                    = LoadType::Displacement;
        Set::Scalar val 					= 0.;
        Set::Scalar rate				= 1.; 
		Set::Scalar max				    = 1.0;
		Set::Scalar init 				= 0.01;
		int 		step				= 0;
		int			current_test			= 0;
		
    } loading;

    struct{
        Set::Field<Set::Scalar> strain;           ///< plastic strains
		Set::Field<Set::Scalar> strain_old;			///< plastic strain at previous time step.
		Set::Scalar norm;
		Set::Scalar error_norm;
		Set::Scalar tol_rel = 1.e-6;
		Set::Scalar tol_abs = 1.e-6;
    } plastic;

    struct{
        brittle_fracture_model_type brittlemodeltype;
        ductile_fracture_model_type ductilemodeltype;
		Set::Field<brittle_fracture_model_type> brittlemodel;
        Set::Field<ductile_fracture_model_type> ductilemodel;
		std::string input_material = "isotropic";
    } material;
    
    struct{
		Set::Scalar bottom_tol 				= 1.0E-6;
		int         interval	 			= 1;
		std::string type					= "single";
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
	} sol;

};
}

#endif