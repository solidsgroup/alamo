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
#include "Model/Solid/Affine/Isotropic.H"
#include "IO/ParmParse.H"
#include "BC/Operator/Elastic/Constant.H"
#include "Solver/Nonlocal/Newton.H"
//extern std::vector<Set::Scalar> alpha;
//  extern int pack_size;


namespace Integrator
{
class Flame : public Integrator::Integrator
{

	using model_type = Model::Solid::Affine::Isotropic;
	
public:
    Flame ();
 

protected:
    void Initialize (int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void Advance (int lev, amrex::Real time, amrex::Real dt) override;
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
private:

    Set::Field<Set::Scalar> Temp_mf;
    Set::Field<Set::Scalar> Temp_old_mf;
    Set::Field<Set::Scalar> Eta_mf;
    Set::Field<Set::Scalar> Eta_old_mf;
    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> field;
    BC::BC<Set::Scalar> *TempBC;
    BC::BC<Set::Scalar> *EtaBC;
    //IC::Voronoi *VoronoiIC;
    IC::IC *PhiIC;

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
    Set::Scalar fmod_ap;               
    Set::Scalar fmod_htpb; 
    Set::Scalar fs_comb;               
    Set::Scalar fmod_comb; 
    Set::Scalar P;
    Set::Scalar r_ap;
    Set::Scalar n_ap;
    Set::Scalar r_htpb;
    Set::Scalar n_htpb;
    Set::Scalar r_comb;
    Set::Scalar n_comb;

    int fs_number = 1;
    //Set::Scalar fs_min = -1.0;
    //Set::Scalar fs_max = 1.0;
    //Set::Scalar mean = 0.1;
    Set::Scalar volume_fraction = 0.1;
    //Set::Scalar std_deviation = 0.05;
    Set::Scalar R_min = 0.1;
    Set::Scalar R_max = 0.2;
    Set::Scalar R_mean = 0.075;
    Set::Scalar verbosity = 0.1;
 
   
    Set::Scalar m_refinement_criterion = 0.001;
    Set::Scalar t_refinement_criterion = 0.001;

    IC::IC *EtaIC;

    struct {
        Set::Field<Set::Scalar> disp_mf;
        Set::Field<Set::Scalar> rhs_mf;
        Set::Field<Set::Scalar> res_mf;
        Set::Field<Set::Scalar> stress_mf;
        Set::Field<model_type>  model_mf;
        model_type model_ap, model_htpb, model_void;
        int interval = 0;
        BC::Operator::Elastic::Constant bc;
        Solver::Nonlocal::Newton<model_type> *solver;//(elastic.op);
    } elastic;
};
}

#endif
