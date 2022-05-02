#ifndef INTEGRATOR_FLAME_H
#define INTEGRATOR_FLAME_H

#include <string>
#include <limits>
#include <memory>


#include "Util/Util.H"
#include "Integrator/Integrator.H"
#include "Integrator/Mechanics.H"
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

class Flame : virtual public Integrator::Integrator, 
              virtual public MechanicsBase<Model::Solid::Affine::Isotropic>
{

    using model_type = Model::Solid::Affine::Isotropic;

public:
    Flame ();
    Flame (IO::ParmParse &pp);

    static void Parse(Flame &value, IO::ParmParse &pp);


protected:
    void Initialize (int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void Advance (int lev, amrex::Real time, amrex::Real dt) override;
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
    virtual void UpdateModel(int a_step) override;
protected:

    Set::Field<Set::Scalar> temp_mf;
    Set::Field<Set::Scalar> temp_old_mf;
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;
    Set::Field<Set::Scalar> Mob_mf; 
    Set::Field<Set::Scalar> Mob_old_mf;
    Set::Field<Set::Scalar> mdot_mf;
    Set::Field<Set::Scalar> heat_mf;
    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> field;
    Set::Field<Set::Scalar> alpha_mf;
    
    BC::BC<Set::Scalar> *bc_temp;
    BC::BC<Set::Scalar> *bc_eta;
    IC::IC *ic_phi;


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

        Set::Scalar m0 = 1.0;
        Set::Scalar Ea = 0.1;
        Set::Scalar T0 = 0.1;
    } pf;

    struct {
        int on = 0;
        Set::Scalar rho_ap, rho_htpb;         
        Set::Scalar k_ap, k_htpb, k0, k_comb;             
        Set::Scalar cp_ap, cp_htpb;           
        Set::Scalar q_ap, q_htpb, q_comb;
        Set::Scalar q0=1.0;
        Set::Scalar temperature_delay = 0.01;
        Set::Scalar ae_ap, ae_htpb, ae_comb;
        Set::Scalar bound=0.0, addtemp;
        Set::Scalar beta1, beta2, MinTemp; 
        Set::Scalar A_ap, A_htpb, A_comb;    
    } thermal;

    struct {
        Model::Solid::Affine::Isotropic model_ap, model_htpb, model_void;
    } elastic;

    using MechanicsBase<Model::Solid::Affine::Isotropic>::m_type;
};
}

#endif
