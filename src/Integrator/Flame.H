#ifndef INTEGRATOR_FLAME_H
#define INTEGRATOR_FLAME_H

#include <string>
#include <limits>
#include <memory>


#include "Util/Util.H"
#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "IC/Voronoi.H"
#include "IC/PackedSpheres.H"
#include "IC/IC.H"
#include "IC/Wedge.H"
#include "IC/Constant.H"
#include "Set/Set.H"

 //extern std::vector<Set::Scalar> alpha;
 //  extern int pack_size;


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

	Set::Field<Set::Scalar> Temp_mf;
	Set::Field<Set::Scalar> Temp_old_mf;
	Set::Field<Set::Scalar> Eta_mf;
	Set::Field<Set::Scalar> Eta_old_mf;
	Set::Field<Set::Scalar> FlameSpeed_mf;
	Set::Field<Set::Scalar> field;
	BC::BC<Set::Scalar> *TempBC;
	BC::BC<Set::Scalar> *EtaBC;
	IC::Voronoi *VoronoiIC;
	IC::PackedSpheres *PackedSpheresIC;

	// int number_of_fabs = 1;
	// int number_of_components = 2;

	Set::Scalar M;                  // Burn speed parameter
	Set::Scalar kappa;              // Burn width parameter
	Set::Scalar w1, w12, w0;        // Energy (before, during, after)
	Set::Scalar rho1, rho0;         // Density (before and after)
	Set::Scalar ka, kh, k0;             // Thermal conductivity (before and after)
	Set::Scalar cp1, cp0;           // Specific heat (before and after)
	Set::Scalar qdotburn;           // Heat lost in flame
	Set::Scalar delA, delH;       // Thermal flux of each material


    Set::Scalar fs_ap;               
	Set::Scalar fs_htpb; 
	Set::Scalar P;
	Set::Scalar r;
	Set::Scalar n;
	
	int fs_number = 1;
	Set::Scalar fs_min = -1.0;
	Set::Scalar fs_max = 1.0;
	//Set::Scalar mean = 0.1;
	Set::Scalar volume_fraction = 0.1;
	//Set::Scalar std_deviation = 0.05;
	Set::Scalar R_min = 0.1;
	Set::Scalar R_max = 0.2;
 
   
	Set::Scalar m_refinement_criterion = 0.001;
	
	IC::IC *EtaIC;

};
}

#endif
