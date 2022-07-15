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
    void Integrate(int amrlev, Set::Scalar time, int step,
                    const amrex::MFIter &mfi, const amrex::Box &box) override;
private:

    Set::Field<Set::Scalar> Temp_mf;
    Set::Field<Set::Scalar> Temp_old_mf;
    Set::Field<Set::Scalar> Eta_mf;
    Set::Field<Set::Scalar> Eta_old_mf;
    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> field;
    BC::BC<Set::Scalar> *TempBC;
    BC::BC<Set::Scalar> *EtaBC;
    IC::IC *EtaIC;
    IC::IC *PhiIC;


    int fs_number = 1;
 
   
    Set::Scalar m_refinement_criterion = 0.001;
    Set::Scalar t_refinement_criterion = 0.001;
    
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
        Set::Scalar rho1, rho0;         
        Set::Scalar ka, kh, k0;             
        Set::Scalar cp1, cp0;           
        Set::Scalar delA, delH;       
    } thermal;

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


    // Integrated variables
    Set::Scalar volume = 0.0;
    Set::Scalar area = 0.0;
};
}

#endif
