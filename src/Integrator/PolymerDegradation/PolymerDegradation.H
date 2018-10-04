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
#include "IC/IC.H"
#include "IC/Constant/Constant.H"
#include "Numeric/Interpolator/Linear/Linear.H"
#include "Util/Util.H"
#include "IO/FileNameParse.H"


#include "Operator/Elastic/Elastic.H"
#include "Model/Solid/Elastic/Isotropic/Degradable/Isotropic.H"


#define TEMP_OLD(i,j,k) Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))
#define TEMP(i,j,k) Temp_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))
#define WATER_OLD(i,j,k) water_conc_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))
#define WATER(i,j,k) water_conc_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))
#define ETA_OLD(i,j,k,m) eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m)
#define ETA_NEW(i,j,k,m) eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m)

namespace Integrator
{

///
/// \class PhaseFieldMicrostructure
/// \brief Microstructure evolution with grain boundary anisotropy
///
/// Solve the Allen-Cahn evolution equation for microstructure with parameters \f$\eta_1\ldots\eta_n\f$,
/// where n corresponds to the number of grains.
///
class PolymerDegradation : public Integrator::Integrator
{
public:
	PolymerDegradation();

protected:

	/// \fn    Advance
	/// \brief Evolve phase field in time
	void Advance (int lev, Real time, Real dt);

	void Initialize (int lev);

	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

	void TimeStepBegin(amrex::Real time, int iter);

private:

	int number_of_ghost_cells = 2;

	// Degradation variable
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > eta_new; 		///< Degradation variable for the __current__timestep
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > eta_old;		///< Degradation variable for the __previous__timestep

	// For water induced degradation
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > water_conc;		///< Water concentration for the __current__timestep
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > water_conc_old;	///< Water concentration for the __previous__timestep

	// For temperature induced degradation
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > Temp;			///< Temperature for the __current__timestep
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > Temp_old;		///< Temperature for the __previous__timestep		

	// For mechanical degradation to be studied later
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > displacement; ///< Multicomponent field variable storing pointwise displacements
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > strain; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress_vm; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > rhs; ///< Multicomponent field variable storing pointwise displacements
	
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > energy; ///< Single component variable storing strain energy
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > energies; ///< Multicomponent field variable storing individual strain energies

	Set::Scalar M, mu, gamma, sigma0, l_gb, beta;
	int anisotropy = 0;

	amrex::Real anisotropy_tstart;

	std::string ic_type, gb_type, filename;

	amrex::Real damp;

	// Water diffusion parameters
	bool 		water_diffusion_on =		false;
	Set::Scalar 	water_diffusivity =		1.0;
	Set::Scalar 	water_refinement_threshold =	0.01;
	std::string 	water_ic_type;
	IC::IC		*water_ic;
	BC::BC		*water_bc;

	// Thermal diffusion parameters
	bool		heat_diffusion_on =		false;
	Set::Scalar 	thermal_diffusivity =		1.0;
	Set::Scalar 	thermal_refinement_threshold =	0.01;
	std::string	thermal_ic_type;
	IC::IC		*thermal_ic;
	BC::BC		*thermal_bc;

	// Damage parameters
	std::string	damage_type;
	bool		damage_anisotropy = 0;
	int		number_of_eta = 1;
	Set::Scalar	damage_refinement_threshold = 0.01;
	std::string 	eta_ic_type;
	IC::IC		*eta_ic;
	BC::BC		*eta_bc;

	// Damage model: water relaxation
	int number_of_terms = 4;
	amrex::Vector<Set::Scalar> d_i;
	Set::Scalar d_final = 1.0;
	amrex::Vector<Set::Scalar> tau_i;
	amrex::Vector<Set::Scalar> t_start_i;

	// Elasticity parameters
	bool        elastic_on = false;
	int         elastic_int = 1;
	std::string elastic_type;
	int         elastic_max_iter = 200;
	int         elastic_max_fmg_iter = 0;
	int         elastic_verbose = 3;
	int         elastic_cgverbose = 3;
	Set::Scalar elastic_tol_rel = 1.0E-5;
	Set::Scalar elastic_tol_abs = 1.0E-5;
	Set::Scalar elastic_tstart = 0.0;
	Set::Scalar elastic_tend = 0.0;
	std::string bottom_solver = "bicgstab";
	int linop_maxorder 	     = 2;
	bool	    elastic_use_fsmooth = false;

	// Elastic BC
	BC::BC *elastic_bc;
	std::map<std::string,Operator::Elastic::BC> bc_map;
	std::array<Operator::Elastic::BC,AMREX_SPACEDIM> bc_x_lo, bc_x_hi;
#if AMREX_SPACEDIM > 1
	std::array<Operator::Elastic::BC,AMREX_SPACEDIM> bc_y_lo, bc_y_hi;
#if AMREX_SPACEDIM > 2
	std::array<Operator::Elastic::BC,AMREX_SPACEDIM> bc_z_lo, bc_z_hi;
#endif
#endif
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
	Operator::Elastic::Elastic<Model::Solid::Elastic::Degradable::Isotropic::Isotropic> *elastic_operator;

	// Material model
	amrex::Vector<amrex::FabArray<amrex::BaseFab<Model::Solid::Elastic::Degradable::Isotropic::Isotropic> > > model;
	Model::Solid::Elastic::Degradable::Isotropic::Isotropic modeltype;
	std::string 	input_material = "isotropic";

};
}
#endif
