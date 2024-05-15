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
        delete bc_rho;
        delete bc_E;
        delete bc_M;
        delete bc_v;
        delete bc_p;
        delete bc_eta;

        delete ic_Velocity;
        delete ic_Pressure;

        delete ic_rhoInterface;
        delete ic_vInjected;
        delete ic_q;

        delete ic_eta;
        delete ic_etadot;
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

    Set::Field<Set::Scalar> Density_mf;
    Set::Field<Set::Scalar> Density_old_mf;

    Set::Field<Set::Scalar> Energy_mf;
    Set::Field<Set::Scalar> Energy_old_mf;

    Set::Field<Set::Scalar> Momentum_mf;
    Set::Field<Set::Scalar> Momentum_old_mf;

    Set::Field<Set::Scalar> Velocity_mf;
    Set::Field<Set::Scalar> Pressure_mf;

protected:
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;
private:
    Set::Field<Set::Scalar> etadot_mf;

    Set::Field<Set::Scalar> Vorticity_mf;

    Set::Field<Set::Scalar> rhoInterface_mf;
    Set::Field<Set::Scalar> vInjected_mf;
    Set::Field<Set::Scalar> q_mf;
    Set::Field<Set::Scalar> Source_mf;

    BC::BC<Set::Scalar>* bc_rho = nullptr;
    BC::BC<Set::Scalar>* bc_E = nullptr;
    BC::BC<Set::Scalar>* bc_M = nullptr;
    BC::BC<Set::Scalar>* bc_v = nullptr;
    BC::BC<Set::Scalar>* bc_p = nullptr;
    BC::BC<Set::Scalar>* bc_eta = nullptr;

    IC::IC* ic_Velocity = nullptr;
    IC::IC* ic_Pressure = nullptr;

    IC::IC* ic_rhoInterface = nullptr;
    IC::IC* ic_vInjected = nullptr;
    IC::IC* ic_q = nullptr;

    IC::IC* ic_eta = nullptr;
    IC::IC* ic_etadot = nullptr;

    BC::Nothing bc_nothing;

    Set::Scalar c_max = 0.0;
    Set::Scalar vx_max = 0.0;
    Set::Scalar vy_max = 0.0;

    Set::Scalar r_refinement_criterion, e_refinement_criterion, m_refinement_criterion, eta_refinement_criterion, omega_refinement_criterion;
    Set::Scalar gamma, cfl, mu;
    Set::Scalar rho_fluid, rho_solid, v_solid;
    Set::Scalar epsilon, Ldot_active, peak_grad_eta;
};
}

#endif
