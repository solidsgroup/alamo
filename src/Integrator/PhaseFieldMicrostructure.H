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

#include "Operator/CellElastic/PolyCrystal/PolyCrystal.H"
#include "Operator/CellElastic/Isotropic/Isotropic.H"
#include "Operator/CellElastic/PolyCrystal/Isotropic/Isotropic.H"
#include "Operator/CellElastic/PolyCrystal/Cubic/Cubic.H"
#include "Operator/CellElastic/Cubic/Cubic.H"
#include "Model/Solid/LinearElastic/Cubic.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Model/Solid/LinearElastic/Laplacian.H"
#include "Operator/Elastic.H"

namespace Integrator
{

using model_type = Model::Solid::LinearElastic::Isotropic;


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
	int number_of_ghost_cells = 2;

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > eta_new_mf; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __current__ timestep
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > eta_old_mf; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __previous__ timestep
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > etas_mf; ///< Color code order params
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > displacement; ///< Multicomponent field variable storing pointwise displacements
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > residual; ///< Multicomponent field variable storing pointwise displacements
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > strain; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress_vm; ///< Multicomponent field variable storing pointwise strains (6 components)
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > body_force; ///< Multicomponent field variable storing pointwise displacements
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > energy; ///< Single component variable storing strain energy
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > energies; ///< Multicomponent field variable storing individual strain energies

	BC::BC *mybc;

	amrex::Real M, mu, gamma, sigma0, l_gb, beta;
	int anisotropy = 0;
	amrex::Real anisotropy_timestep;
	amrex::Real anisotropy_tstart;

	std::string ic_type, gb_type, filename;

	Model::Interface::GB::GB *boundary;

	IC::IC *ic;

	bool sdf_on = false;
	amrex::Vector<amrex::Real> sdf;

	Set::Scalar volume = 5;

	struct
	{
		bool        on       = false;
		Grid        grid = Grid::Cell;
		int         interval = 1;
		//std::string type;
		int         max_iter = 20;
		int         fixed_iter = 0;
		int         max_fmg_iter = 20;
		int         bottom_max_iter = 20;
		int         verbose = 3;
		int         cgverbose = 3;
		amrex::Real tol_rel = 0.0;
		amrex::Real tol_abs = 1.0E-10;
		amrex::Real tstart = 0.0;
		amrex::Vector<amrex::Real> load_t;
		amrex::Vector<amrex::Real> load_disp;

		// Node only
		amrex::Vector<amrex::Real> AMREX_D_DECL(bc_xlo,bc_ylo,bc_zlo);
		amrex::Vector<amrex::Real> AMREX_D_DECL(bc_xhi,bc_yhi,bc_zhi);
		std::array<Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> AMREX_D_DECL(bctype_xlo, bctype_ylo, bctype_zlo);
		std::array<Operator::Elastic<model_type>::BC,AMREX_SPACEDIM> AMREX_D_DECL(bctype_xhi, bctype_yhi, bctype_zhi);


		Operator::Elastic<model_type> *op;
		//BC::BC *bc;

	} elastic;
	//#endif

	Operator::CellElastic::PolyCrystal::PolyCrystal *cellelastic_operator;
	Operator::Elastic<Model::Solid::LinearElastic::Cubic> *elastic_operator;
	BC::BC *elastic_bc;

	std::vector<Operator::CellElastic::PolyCrystal::PolyCrystalModel *> models;
	
	const std::array<amrex::IntVect,AMREX_SPACEDIM> dx = {{AMREX_D_DECL(amrex::IntVect(AMREX_D_DECL(1,0,0)),
									    amrex::IntVect(AMREX_D_DECL(0,1,0)),
									    amrex::IntVect(AMREX_D_DECL(0,0,1)))}};

};
}
#endif
