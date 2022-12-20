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
#include "IC/Expression.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/Isotropic.H"
#include "IO/ParmParse.H"
#include "BC/Operator/Elastic/Constant.H"
#include "Solver/Nonlocal/Newton.H"



namespace Integrator
{

class Flame : virtual public Integrator::Integrator, 
              virtual public Base::Mechanics<Model::Solid::Affine::Isotropic>
{

    using model_type = Model::Solid::Affine::Isotropic;

public:
    Flame ();
    Flame (IO::ParmParse &pp);

    static void Parse(Flame &value, IO::ParmParse &pp);


protected:
    void Initialize (int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void Advance (int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
                    const amrex::MFIter &mfi, const amrex::Box &box) override;
    void UpdateModel(int a_step) override;
private:

    Set::Field<Set::Scalar> temp_mf;
    Set::Field<Set::Scalar> temp_old_mf;
    Set::Field<Set::Scalar> temps_mf;
    Set::Field<Set::Scalar> temps_old_mf;
  
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;
    Set::Field<Set::Scalar> mob_mf; 
    Set::Field<Set::Scalar> mdot_mf;
    Set::Field<Set::Scalar> deta_mf;
    Set::Field<Set::Scalar> deta2_mf; 

    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> field;
    Set::Field<Set::Scalar> alpha_mf;
    Set::Field<Set::Scalar> heatflux_mf;  

    Set::Field<Set::Scalar> laser_mf;
  
    BC::BC<Set::Scalar> *bc_temp;
    BC::BC<Set::Scalar> *bc_eta;
    IC::IC *ic_phi;
    IC::IC *ic_laser;
 
    Set::Scalar m_refinement_criterion = 0.001;
    Set::Scalar t_refinement_criterion = 0.001;
    Set::Scalar t_refinement_restriction = 0.1;
    Set::Scalar zeta = 0.05;
    Set::Scalar zeta_0 = 0.05;
    Set::Scalar small = 1.0e-8;
    Set::Scalar base_time =1.0e-4;
    IC::IC *ic_eta;
   
    struct {
      Set::Scalar evolve = 1;
      Set::Scalar round = 1.0;
    } conditional; 

  
   struct {
      Set::Scalar gamma = 1.0;
      Set::Scalar eps = 0.0;
      Set::Scalar lambda = 0.0;
      Set::Scalar kappa = 0.0;
      Set::Scalar w1 = 0.0, w12 = 0.0, w0 = 0.0;
      Set::Scalar min_eta = 0.001; 
      int evolve = 1;
    } pf;

    struct{

      Set::Scalar P=1.0;
      Set::Scalar a1 = 0.0, a2 = 0.0, a3 = 0.0;
      Set::Scalar b1 = 0.0, b2 = 0.0, b3 = 0.0;
      Set::Scalar c1 = 0.0;
      Set::Scalar E1, E2;
    } pressure;


    struct {
      int on = 0;
      Set::Scalar rho_ap, rho_htpb;         
      Set::Scalar k_ap , k_htpb;
      Set::Scalar cp_ap, cp_htpb;           
      Set::Scalar q0=0.0;
      Set::Scalar bound = 0.0;
      Set::Scalar m_ap, m_htpb, m_comb;  
      Set::Scalar E_ap, E_htpb;
      Set::Scalar hc = 1.0, w1 = 1.0;
      Set::Scalar mlocal_ap = 0.0, mlocal_htpb = 0.0, mlocal_comb = 0.0;
      Set::Scalar T_fluid = 300.0;
      int controlj1 = 0, controlj2 = 0, laseron = 0;
      int lasershut = 1;
      Set::Scalar laseroff = 1.0;
    } thermal;

    struct {
      Model::Solid::Affine::Isotropic model_ap, model_htpb, model_void;
    } elastic;

    // Integrated variables
    Set::Scalar volume = 0.0;
    Set::Scalar area = 0.0;
};
}

#endif
