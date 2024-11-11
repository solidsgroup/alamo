#ifndef INTEGRATOR_HYDRO_H
#define INTEGRATOR_HYDRO_H

#include <string>
#include <limits>
#include <memory>

#include "Integrator/Integrator.H"
#include "Integrator/AllenCahn.H"
#include "IO/ParmParse.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Nothing.H"
#include "IC/IC.H"

#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/Laminate.H"
#include "IC/Expression.H"

#if AMREX_SPACEDIM==2

namespace Integrator
{
class Hydro : virtual public AllenCahn, virtual public Integrator
{
public:

    Hydro() {};
    Hydro(IO::ParmParse& pp);
    static void Parse(Hydro& value, IO::ParmParse& pp);
    virtual ~Hydro()
    {
        delete density_bc;
        delete momentum_bc;
        delete energy_bc;
        delete eta_bc;

        delete velocity_ic;
        delete pressure_ic;

        delete ic_rho_injected;
        delete ic_mdot;
        delete ic_q;

        delete eta_ic;
        delete etadot_ic;
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
        IC::IC* density_ic   = nullptr;
        IC::IC* momentum_ic  = nullptr;
        IC::IC* energy_ic    = nullptr;
    } solid;

protected:
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;
private:
    Set::Field<Set::Scalar> etadot_mf;

    Set::Field<Set::Scalar> vorticity_mf;

    Set::Field<Set::Scalar> rho_injected_mf;
    Set::Field<Set::Scalar> mdot_mf;
    Set::Field<Set::Scalar> q_mf;
    Set::Field<Set::Scalar> Source_mf;
    Set::Field<Set::Scalar> flux_mf;

    BC::BC<Set::Scalar>* density_bc = nullptr;
    BC::BC<Set::Scalar>* momentum_bc = nullptr;
    BC::BC<Set::Scalar>* energy_bc = nullptr;
    BC::BC<Set::Scalar>* eta_bc = nullptr;

    IC::IC* velocity_ic = nullptr;
    IC::IC* pressure_ic = nullptr;
    IC::IC* density_ic  = nullptr;

    IC::IC* ic_rho_injected = nullptr;
    IC::IC* ic_mdot = nullptr;
    IC::IC* ic_q = nullptr;

    IC::IC* eta_ic = nullptr;
    IC::IC* etadot_ic = nullptr;

    BC::Nothing bc_nothing;

    BC::Constant neumann_bc_1 = BC::Constant::ZeroNeumann(1);
    BC::Constant neumann_bc_D = BC::Constant::ZeroNeumann(AMREX_SPACEDIM);

    Set::Scalar c_max = 0.0;
    Set::Scalar vx_max = 0.0;
    Set::Scalar vy_max = 0.0;

    Set::Scalar r_refinement_criterion=NAN, e_refinement_criterion=NAN, m_refinement_criterion=NAN, eta_refinement_criterion=NAN, omega_refinement_criterion=NAN, gradu_refinement_criterion=NAN, p_refinement_criterion=NAN, rho_refinement_criterion;
    Set::Scalar gamma=NAN, cfl=NAN, cfl_v=NAN, mu=NAN, Lfactor=NAN, Pfactor=NAN, pref=NAN;

    Set::Scalar small=NAN;
};
}

#endif
#endif
