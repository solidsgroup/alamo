#ifndef INTEGRATOR_FLAME_H
#define INTEGRATOR_FLAME_H

#include <string>
#include <limits>
#include <memory>


#include "Util/Util.H"
#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "IC/Voronoi.H"
#include "IC/IC.H"
#include "IC/Wedge.H"
#include "IC/Constant.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/Isotropic.H"
#include "IO/ParmParse.H"
#include "BC/Operator/Elastic/Constant.H"
#include "Solver/Nonlocal/Newton.H"



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
    Set::Field<Set::Scalar> Mob_mf; 
    Set::Field<Set::Scalar> Mob_old_mf;
    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> field;
    BC::BC<Set::Scalar> *TempBC;
    BC::BC<Set::Scalar> *EtaBC;
    IC::IC *PhiIC;


    int fs_number = 1;
 
   
    Set::Scalar m_refinement_criterion = 0.001;
    Set::Scalar t_refinement_criterion = 0.001;
    
    IC::IC *EtaIC;

    struct {
       // Set::Scalar M = 0.0;
        Set::Scalar gamma = 1.0;
        Set::Scalar eps = 0.0;
        Set::Scalar lambda = 0.0;
        Set::Scalar kappa = 0.0;
        Set::Scalar w1 = 0.0, w12 = 0.0, w0 = 0.0;
        Set::Scalar P = 0.0;
        Set::Scalar r_ap = 0.0;
        Set::Scalar n_ap = 0.0;
        Set::Scalar r_htpb = 0.0;
        Set::Scalar n_htpb = 0.0;
        Set::Scalar r_comb = 0.0;
        Set::Scalar n_comb = 0.0;
    } pf;

    struct {
        int on = 0;
        Set::Scalar rho_ap, rho_htpb;         
        Set::Scalar k_ap, k_htpb, k0, k_comb;             
        Set::Scalar cp_ap, cp_htpb;           
        Set::Scalar q_ap, q_htpb, q_comb;
        Set::Scalar temperature_delay = 0.01;
        Set::Scalar ae_ap, ae_htpb, ae_comb;
        Set::Scalar bound, addtemp, k_set; 
        Set::Scalar reaction_temperature;      
    } thermal;

    struct {
        //Set::Scalar tol_rel = 1E-8;
        //Set::Scalar tol_abs = 1E-8;
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
