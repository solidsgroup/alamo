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
#include "IC/PSRead.H" 
#include "Set/Set.H"
#include "Model/Solid/Affine/Isotropic.H"
#include "Model/Solid/Elastic/NeoHookeanPredeformed.H"
#include "IO/ParmParse.H"
#include "BC/Operator/Elastic/Constant.H"
#include "Solver/Nonlocal/Newton.H"



namespace Integrator
{

class Flame : virtual public Integrator::Integrator,
    virtual public Base::Mechanics<Model::Solid::Elastic::NeoHookeanPredeformed>
{

    using model_type = Model::Solid::Elastic::NeoHookeanPredeformed;

public:
    Flame();
    Flame(IO::ParmParse& pp);

    static void Parse(Flame& value, IO::ParmParse& pp);


protected:
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
        const amrex::MFIter& mfi, const amrex::Box& box) override;
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

    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> phicell_mf;
    Set::Field<Set::Scalar> field;
    Set::Field<Set::Scalar> alpha_mf;
    Set::Field<Set::Scalar> heatflux_mf;

    Set::Field<Set::Scalar> laser_mf;
    //Set::Field<Set::Scalar> rhs_mf;
   
  Set::Field<Set::Scalar> pressure_mf;
  
    BC::BC<Set::Scalar>* bc_temp;
    BC::BC<Set::Scalar>* bc_eta;
    IC::IC* ic_phi;
    IC::IC* ic_phicell;
    IC::IC* ic_laser;

    Set::Scalar phi_refinement_criterion = std::numeric_limits<Set::Scalar>::infinity();
    Set::Scalar m_refinement_criterion = 0.001;
    Set::Scalar t_refinement_criterion = 0.001;
    Set::Scalar t_refinement_restriction = 0.1;
    Set::Scalar zeta = 1.0e-5;
    Set::Scalar zeta_0 = 1.0e-5;
    Set::Scalar small = 1.0e-8;
    Set::Scalar base_time = 1.0e-4;
    IC::IC* ic_eta;
    int ghost_count = 2;
    int homogeneousSystem = 0;
    bool plot_field = false;
    int variable_pressure = 0;
  
    struct {
        Set::Scalar gamma = 1.0;
        Set::Scalar eps = 0.0;
        Set::Scalar lambda = 0.0;
        Set::Scalar kappa = 0.0;
        Set::Scalar w1 = 0.0, w12 = 0.0, w0 = 0.0;
        Set::Scalar min_eta = 0.001;
        int evolve = 1;
    } pf;

    struct {
        Set::Scalar P = 1.0;
        Set::Scalar gamma = 1.4;
        Set::Scalar Energy = 1;

        struct {
            Set::Scalar a1 = 0.0, a2 = 0.0, a3 = 0.0;
            Set::Scalar b1 = 0.0, b2 = 0.0, b3 = 0.0;
            Set::Scalar c1 = 0.0;
            int mob_ap = 0;
            int dependency = 1;
            Set::Scalar h1 = 1.81, h2 = 1.34;

        } arrhenius;
        struct {
            Set::Scalar r_ap = 0.0, r_htpb = 0.0, r_comb = 0.0;
            Set::Scalar n_ap = 0.0, n_htpb = 0.0, n_comb = 0.0;
        } power;
    } pressure;

    struct {
        bool on = 0;
        Set::Scalar rho_ap, rho_htpb;
        Set::Scalar k_ap, k_htpb;
        Set::Scalar cp_ap, cp_htpb;
        Set::Scalar q0 = 0.0;
        Set::Scalar bound = 0.0;
        Set::Scalar m_ap, m_htpb, m_comb;
        Set::Scalar E_ap, E_htpb;
        Set::Scalar hc = 1.0, w1 = 1.0;
        Set::Scalar modeling_ap = 1.0, modeling_htpb = 1.0;
        Set::Scalar mlocal_ap = 0.0, mlocal_htpb = 0.0, mlocal_comb = 0.0;
        Set::Scalar T_fluid = 300.0;
        Set::Scalar massfraction = 0.8;
        Set::Scalar disperssion1 = 0.93;
        Set::Scalar disperssion2 = 920.;
        Set::Scalar disperssion3 = 2418.29;
    } thermal;

    struct {
        int on = 0;
        Set::Scalar Tref = 300.0;
        model_type model_ap, model_htpb;
        Set::Scalar traction = 0.0;
        int phirefinement = 1;
    } elastic;

    // Integrated variables
    Set::Scalar volume = 0.0;
    Set::Scalar area = 0.0;

    BC::BC<Set::Scalar>* bc_psi = nullptr;



};
}

#endif
