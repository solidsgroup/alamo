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
    Hydro(IO::ParmParse& pp);
    static void Parse(Hydro& value, IO::ParmParse& pp);
    virtual ~Hydro()
    {
        delete density_bc;
        delete momentum_bc;
        delete energy_bc;
        delete eta_bc;

        delete density_ic;
        delete velocity_ic;
        delete pressure_ic;

        delete ic_m0;
        delete ic_u0;
        delete ic_q;

        delete eta_ic;
        delete etadot_ic;

        delete solid.momentum_ic;
        delete solid.density_ic;
        delete solid.energy_ic;

        delete roesolver;
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
private:

    Set::Field<Set::Scalar> density_mf;
    Set::Field<Set::Scalar> density_old_mf;

    Set::Field<Set::Scalar> energy_mf;
    Set::Field<Set::Scalar> energy_old_mf;

    Set::Field<Set::Scalar> momentum_mf;
    Set::Field<Set::Scalar> momentum_old_mf;

    Set::Field<Set::Scalar> velocity_mf;
    Set::Field<Set::Scalar> pressure_mf;

    struct {
        Set::Field<Set::Scalar> density_mf;
        Set::Field<Set::Scalar> momentum_mf;
        Set::Field<Set::Scalar> energy_mf;
        IC::IC<Set::Scalar>* density_ic   = nullptr;
        IC::IC<Set::Scalar>* momentum_ic  = nullptr;
        IC::IC<Set::Scalar>* energy_ic    = nullptr;
    } solid;

protected:
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;
private:
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

    IC::IC<Set::Scalar>* velocity_ic = nullptr;
    IC::IC<Set::Scalar>* pressure_ic = nullptr;
    IC::IC<Set::Scalar>* density_ic  = nullptr;

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
    Set::Scalar gamma=NAN, cfl=NAN, cfl_v=NAN, mu=NAN,  pref=NAN;

    Set::Scalar small=NAN;
    Set::Scalar cutoff=NAN;
    Set::Scalar lagrange=NAN;

    Solver::Local::Riemann::Roe *roesolver = nullptr;

};
}

#endif
#endif
