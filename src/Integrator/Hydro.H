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

    Set::Field<Set::Scalar> etaDensity_mf;
    Set::Field<Set::Scalar> etaDensity_old_mf;

    Set::Field<Set::Scalar> etaEnergy_mf;
    Set::Field<Set::Scalar> etaEnergy_old_mf;

    Set::Field<Set::Scalar> etaMomentum_mf;
    Set::Field<Set::Scalar> etaMomentum_old_mf;

    Set::Field<Set::Scalar> Velocity_mf;
    Set::Field<Set::Scalar> etaPressure_mf;

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

    BC::BC<Set::Scalar>* bc_rho;
    BC::BC<Set::Scalar>* bc_E;
    BC::BC<Set::Scalar>* bc_M;
    BC::BC<Set::Scalar>* bc_v;
    BC::BC<Set::Scalar>* bc_p;
    BC::BC<Set::Scalar>* bc_omega;
    BC::BC<Set::Scalar>* bc_eta;

    IC::IC* ic_Velocity;
    IC::IC* ic_etaPressure;

    IC::IC* ic_rhoInterface;
    IC::IC* ic_vInjected;
    IC::IC* ic_q;

    IC::IC* ic_eta;
    IC::IC* ic_etadot;

BC::Nothing bc_nothing;

    Set::Scalar c_max = 0.0;
    Set::Scalar vx_max = 0.0;
    Set::Scalar vy_max = 0.0;

    Set::Scalar r_refinement_criterion, e_refinement_criterion, m_refinement_criterion, eta_refinement_criterion, omega_refinement_criterion;
    Set::Scalar gamma, cfl, mu;
    Set::Scalar rho_fluid;
    Set::Scalar epsilon, Ldot_active;
};
}

#endif
