///
/// \file PhaseFieldMicrostructure.H
///
#ifndef INTEGRATOR_PHASEFIELDMICROSTRUCTURE_H
#define INTEGRATOR_PHASEFIELDMICROSTRUCTURE_H

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
#include "IC/TabulatedInterface.H"
#include "IC/PerturbedInterface.H"
#include "IC/Voronoi.H"
#include "IC/Sphere.H"

#include "Model/Interface/GB/GB.H"
#include "Model/Interface/GB/Sin.H"
#include "Model/Interface/GB/AbsSin.H"
#include "Model/Interface/GB/Read.H"

#include "Model/Solid/Linear/Cubic.H"
#include "Model/Solid/Affine/Cubic.H"
//#include "Model/Solid/LinearElastic/Cubic.H"
//#include "Model/Solid/LinearElastic/Isotropic.H"
//#include "Model/Solid/Linear/Laplacian.H"
#include "Operator/Elastic.H"

namespace Integrator
{
//using model_type = Model::Solid::LinearElastic::Isotropic;
//using model_type = Model::Solid::LinearElastic::Cubic;
//using model_type = Model::Solid::Linear::Cubic;
using model_type = Model::Solid::Affine::Cubic;

enum RegularizationType{
	Wilmore = 1,
	K12 = 2};

///
/// \class PhaseFieldMicrostructure
/// \brief Microstructure evolution with grain boundary anisotropy
///
/// Solve the Allen-Cahn evolution equation for microstructure with parameters \f$\eta_1\ldots\eta_n\f$,
/// where n corresponds to the number of grains.
///
class PhaseFieldMicrostructure : public Integrator
{
public:
	PhaseFieldMicrostructure();

protected:

	/// \fn    Advance
	/// \brief Evolve phase field in time
	void Advance (int lev, Real time, Real dt) override;

	void Initialize (int lev) override;

	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

	void TimeStepBegin(amrex::Real time, int iter) override;
	void TimeStepComplete(amrex::Real time, int iter) override;
	void Integrate(int amrlev, Set::Scalar time, int step,
		       const amrex::MFIter &mfi, const amrex::Box &box) override;

private:

	int number_of_grains = 2;
	int number_of_ghost_cells = 3;
	Set::Scalar ref_threshold = 0.1;

	// Cell fab
	Set::Field<Set::Scalar> eta_new_mf; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __current__ timestep
	Set::Field<Set::Scalar> eta_old_mf; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __previous__ timestep
	// Node fab
	Set::Field<Set::Scalar> disp_mf; 
	Set::Field<Set::Scalar> rhs_mf; 
	Set::Field<Set::Scalar> stress_mf;
	Set::Field<Set::Scalar> energy_mf;

	BC::BC *mybc;

	//amrex::Real M, mu, gamma, sigma0, l_gb, beta;
	RegularizationType regularization = RegularizationType::K12;
	struct {
		Set::Scalar M;
		Set::Scalar mu;
		Set::Scalar gamma;
		Set::Scalar sigma0;
		Set::Scalar l_gb;
		Set::Scalar elastic_mult = 1.0;
		Set::Scalar elastic_threshold = 0.0;
	} pf;

	struct {
		int on = 0;
		Set::Scalar beta;
		Set::Scalar timestep;
		Set::Scalar tstart;
		int plot_int = -1;
		Set::Scalar plot_dt = -1.0;
		int thermo_int = -1, thermo_plot_int = -1;
		Set::Scalar theta0,sigma0,sigma1;
		Set::Scalar phi0 = 0.0;
		int elastic_int = -1;
	} anisotropy;
	
	struct { 
		int on = 0;
		Set::Scalar tstart = 0.0;
		Set::Scalar vol0 = 0.0;
		Set::Scalar lambda = 0.0;
	} lagrange;

	std::string ic_type, gb_type, filename;

	Model::Interface::GB::GB *boundary;

	IC::IC *ic;

	Set::Scalar volume = 5;
	Set::Scalar area = 0.0;
	Set::Scalar gbenergy = 0.0;
	Set::Scalar realgbenergy = 0.0;
	Set::Scalar regenergy = 0.0;

	struct
	{
		bool        on       = false;
		int         interval = -1;
		int         max_coarsening_level = -1;
		amrex::Real tol_rel = 0.0;
		amrex::Real tol_abs = 1.0E-10;
		amrex::Real tstart = 0.0;
		amrex::Vector<amrex::Real> load_t;
		amrex::Vector<amrex::Real> load_disp;
		std::vector<model_type> model;

		amrex::Vector<amrex::Real> AMREX_D_DECL(bc_xlo,bc_ylo,bc_zlo);
		amrex::Vector<amrex::Real> AMREX_D_DECL(bc_xhi,bc_yhi,bc_zhi);
		std::array<BC::Operator::Elastic<model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bctype_xlo, bctype_ylo, bctype_zlo);
		std::array<BC::Operator::Elastic<model_type>::Type,AMREX_SPACEDIM> AMREX_D_DECL(bctype_xhi, bctype_yhi, bctype_zhi);
		Operator::Elastic<model_type> op;

		BC::Operator::Elastic<model_type> bc;

		Set::Scalar strainenergy = 0.0;
		Set::Scalar force = 0.0;
		Set::Scalar disp = 0.0;
	} elastic;

};
}
#endif
