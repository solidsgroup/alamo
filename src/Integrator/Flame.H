#ifndef INTEGRATOR_FLAME_H
#define INTEGRATOR_FLAME_H

#include <string>
#include <limits>
#include <memory>


#include "Util/Util.H"
#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "IC/Voronoi.H"
#include "IC/Wedge.H"
#include "IC/Constant.H"

namespace Integrator
{
class Flame : public Integrator::Integrator
{
public:
	Flame ();

protected:
	void Initialize (int lev) ;
	void Advance (int lev, amrex::Real time, amrex::Real dt);
	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/);
	void Regrid(int lev, Set::Scalar time) override;
private:

	amrex::Vector<std::unique_ptr<amrex::MultiFab> > Temp;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > Temp_old;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > Eta;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > Eta_old;
	amrex::Vector<std::unique_ptr<amrex::MultiFab> > FlameSpeedFab;
	BC::BC *TempBC;
	BC::BC *EtaBC;
	IC::Voronoi *VoronoiIC;

	// int number_of_fabs = 1;
	// int number_of_components = 2;

	amrex::Real M;                  // Burn speed parameter
	amrex::Real kappa;              // Burn width parameter
	amrex::Real w1, w12, w0;        // Energy (before, during, after)
	amrex::Real rho1, rho0;         // Density (before and after)
	amrex::Real k1, k0;             // Thermal conductivity (before and after)
	amrex::Real cp1, cp0;           // Specific heat (before and after)
	amrex::Real qdotburn;           // Heat lost in flame

	int fs_number = 200;
	amrex::Real fs_min = -1.0;
	amrex::Real fs_max = 1.0;
	
	IC::IC *EtaIC;

};
}

#endif
