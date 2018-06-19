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
#include "IC/IC.H"
#include "IC/Constant/Constant.H"
#include "IC/PerturbedInterface/PerturbedInterface.H"

#include "Operator/Elastic/Degradation/Degradation.H"
#include "Operator/Elastic/Degradation/Isotropic/Isotropic.H"
#include "Operator/Elastic/Degradation/Cubic/Cubic.H"
#include "Operator/Elastic/PolyCrystal/PolyCrystal.H"
#include "Operator/Elastic/Isotropic/Isotropic.H"
#include "Operator/Elastic/PolyCrystal/Isotropic/Isotropic.H"
#include "Operator/Elastic/PolyCrystal/Cubic/Cubic.H"
#include "Operator/Elastic/Cubic/Cubic.H"


#define TEMP_OLD(i,j,k) Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))
#define TEMP(i,j,k) Temp_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))
#define WATER_OLD(i,j,k) water_conc_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))
#define WATER(i,j,k) water_conc_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))

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
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > body_force; ///< Multicomponent field variable storing pointwise displacements
	
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > energy; ///< Single component variable storing strain energy
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > energies; ///< Multicomponent field variable storing individual strain energies

	BC::BC *mybc;

	amrex::Real M, mu, gamma, sigma0, l_gb, beta;
	int anisotropy = 0;

	amrex::Real anisotropy_tstart;

	std::string ic_type, gb_type, filename;

	amrex::Real damp;

	// Water diffusion parameters
	bool 		water_diffusion_on =		false;
	amrex::Real 	water_diffusivity =		1.0;
	amrex::Real 	water_refinement_threshold =	0.01;
	std::string 	water_ic_type;
	IC::IC		*water_ic;
	BC::BC		*water_bc;

	// Thermal diffusion parameters
	bool		heat_diffusion_on =		false;
	amrex::Real 	thermal_diffusivity =		1.0;
	amrex::Real 	thermal_refinement_threshold =	0.01;
	std::string	thermal_ic_type;
	IC::IC		*thermal_ic;
	BC::BC		*thermal_bc;

	// Damage parameters
	std::string	damage_type;
	int		damage_anisotropy = 0;
	int		number_of_eta;
	std::string 	eta_ic_type;
	IC::IC		*eta_ic;
	BC::BC		*eta_bc;

	// Damage model: relaxation
	int number_of_terms = 4;
	amrex::Vector<amrex::Real> E_i;
	amrex::Real E0;
	amrex::Vector<amrex::Real> tau_i;
	amrex::Vector<amrex::Real> t_start_i;

	// Elasticity parameters
	bool        elastic_on = false;
	int         elastic_int = 1;
	std::string elastic_type;
	int         elastic_max_iter = 200;
	int         elastic_max_fmg_iter = 0;
	int         elastic_verbose = 3;
	int         elastic_cgverbose = 3;
	amrex::Real elastic_tol_rel = 0.0;
	amrex::Real elastic_tol_abs = 1.0E-10;
	amrex::Real elastic_tstart = 0.0;
	amrex::Vector<amrex::Real> elastic_load_t;
	amrex::Vector<amrex::Real> elastic_load_disp;
	Operator::Elastic::Degradation::Degradation *elastic_operator;

};
}
#endif
