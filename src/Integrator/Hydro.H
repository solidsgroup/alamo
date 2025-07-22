#ifndef INTEGRATOR_HYDRO_H
#define INTEGRATOR_HYDRO_H

#include "Integrator/Integrator.H"
#include "IO/ParmParse.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Nothing.H"
#include "IC/IC.H"

#include "Set/Base.H"
#include "Solver/Local/Riemann/Roe.H"

#if AMREX_SPACEDIM==2

namespace Integrator
{
class Hydro : virtual public Integrator
{
public:
    static constexpr const char* name = "hydro";

    Hydro() {};

    Hydro(  Set::Field<Set::Scalar> &a_eta,
            Set::Field<Set::Scalar> &a_etaold,
            bool a_invert) :
        eta_mf(&a_eta), eta_old_mf(&a_etaold), invert(a_invert)
    {
        this->managed = true;
    }

    Hydro(IO::ParmParse& pp); 
    static void Parse(Hydro& value, IO::ParmParse& pp);
    virtual ~Hydro()
    {
        if (!managed)
        {
            delete eta_mf;
            delete eta_old_mf;
        }

        delete density_bc;
        delete momentum_bc;
        delete energy_bc;
        delete eta_bc;
        delete tracer_bc;

        delete density_ic;
        delete velocity_ic;
        delete pressure_ic;
        delete tracer_ic;

        delete ic_m0;
        delete ic_u0;
        delete ic_q;

        delete eta_ic;
        delete etadot_ic;

        delete solid.momentum_ic;
        delete solid.density_ic;
        delete solid.energy_ic;

        delete riemannsolver;
    }

protected:
    void Initialize(int lev) override;
    void Mix(int lev);
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;

    virtual void UpdateEta(int lev, Set::Scalar time);
    virtual void UpdateFluxes(int lev, Set::Scalar time);
protected:
    
    void RHS(   int lev, Set::Scalar time, 
                amrex::MultiFab &rho_rhs_mf, 
                amrex::MultiFab &M_rhs_mf, 
                amrex::MultiFab &E_rhs_mf,
                const amrex::MultiFab &rho_mf,
                const amrex::MultiFab &M_mf,
                const amrex::MultiFab &E_mf);

    Set::Field<Set::Scalar> density_mf;
    Set::Field<Set::Scalar> density_old_mf;

    Set::Field<Set::Scalar> energy_mf;
    Set::Field<Set::Scalar> energy_old_mf;

    Set::Field<Set::Scalar> momentum_mf;
    Set::Field<Set::Scalar> momentum_old_mf;

    Set::Field<Set::Scalar> velocity_mf;
    Set::Field<Set::Scalar> pressure_mf;


    Set::Field<Set::Scalar> tracer_mf;
    Set::Field<Set::Scalar> tracer_old_mf;

    struct {
        Set::Field<Set::Scalar> density_mf;
        Set::Field<Set::Scalar> momentum_mf;
        Set::Field<Set::Scalar> energy_mf;
        IC::IC<Set::Scalar>* density_ic   = nullptr;
        IC::IC<Set::Scalar>* momentum_ic  = nullptr;
        IC::IC<Set::Scalar>* energy_ic    = nullptr;
    } solid;

    Set::Field<Set::Scalar>* eta_mf;
    Set::Field<Set::Scalar>* eta_old_mf;
    Set::Field<Set::Scalar> etadot_mf;

    Set::Field<Set::Scalar> vorticity_mf;

    Set::Field<Set::Scalar> m0_mf;
    Set::Field<Set::Scalar> u0_mf;
    Set::Field<Set::Scalar> q_mf;
    Set::Field<Set::Scalar> Source_mf;
    Set::Field<Set::Scalar> flux_mf;

    BC::BC<Set::Scalar>* density_bc = nullptr;
    BC::BC<Set::Scalar>* momentum_bc = nullptr;
    BC::BC<Set::Scalar>* energy_bc = nullptr;
    BC::BC<Set::Scalar>* eta_bc = nullptr;
    BC::BC<Set::Scalar>* tracer_bc = nullptr;
 
    IC::IC<Set::Scalar>* velocity_ic = nullptr;
    IC::IC<Set::Scalar>* pressure_ic = nullptr;
    IC::IC<Set::Scalar>* density_ic  = nullptr;
    IC::IC<Set::Scalar>* tracer_ic  = nullptr;

    IC::IC<Set::Scalar>* ic_m0 = nullptr;
    IC::IC<Set::Scalar>* ic_u0 = nullptr;
    IC::IC<Set::Scalar>* ic_q = nullptr;

    IC::IC<Set::Scalar>* eta_ic = nullptr;
    IC::IC<Set::Scalar>* etadot_ic = nullptr;

    BC::Nothing bc_nothing;

    BC::Constant neumann_bc_1 = BC::Constant::ZeroNeumann(1);
    BC::Constant neumann_bc_D = BC::Constant::ZeroNeumann(AMREX_SPACEDIM);

    Set::Scalar c_max = 0.0;
    Set::Scalar vx_max = 0.0;
    Set::Scalar vy_max = 0.0;

    Set::Scalar eta_refinement_criterion=NAN, omega_refinement_criterion=NAN, gradu_refinement_criterion=NAN, p_refinement_criterion=NAN, rho_refinement_criterion;
    Set::Scalar gamma=NAN, cfl=NAN, cfl_v=NAN, mu=NAN, k=NAN, pref=NAN;

    Set::Scalar small=NAN;
    Set::Scalar cutoff=NAN;
    Set::Scalar lagrange=NAN;
    Set::Scalar lagrange_m0=NAN;

    Set::Scalar nspecies=0;
    Set::Scalar temperature=0.0;
    std::vector<std::string> species;
    std::vector<double> species_mw;
    std::vector<double> species_massf;
    std::vector<double> species_molef;
    std::vector<double> species_mu;
    std::vector<double> species_k;
    std::vector<double> species_cp;

    Solver::Local::Riemann::Riemann *riemannsolver = nullptr;

    enum IntegrationScheme {
        ForwardEuler, SSPRK3, RK4
    };
    IntegrationScheme scheme;

    enum PrescribedFlowMode {
        Absolute, Relative
    };
    PrescribedFlowMode prescribedflowmode;


    Set::Vector g;

    bool managed = false;
    bool invert = false;

    std::vector<bool> mixed;

};
}

#endif
#endif
