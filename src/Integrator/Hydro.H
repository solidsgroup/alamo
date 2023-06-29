#ifndef INTEGRATOR_HYDRO_H
#define INTEGRATOR_HYDRO_H

#include <string>
#include <limits>
#include <memory>

#include "Integrator/Integrator.H"
#include "IO/ParmParse.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "IC/IC.H"

#include "IC/Cylinder.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/Laminate.H"
#include "IC/Expression.H"

namespace Integrator
{
class Hydro: virtual public Integrator::Integrator
{
public:

    //Hydro();
    //Hydro (IO::ParmParse &pp);  
    Hydro(IO::ParmParse& pp) { pp.queryclass(*this); }
    static void Parse(Hydro& value, IO::ParmParse& pp)
    {
        BL_PROFILE("Integrator::Hydro::Hydro()");
        //General Variables Input Read:
        {
            //pp.query("r_refinement_criterion", value.r_refinement_criterion);
            //pp.query("e_refinement_criterion", value.e_refinement_criterion);
            //pp.query("m_refinement_criterion", value.m_refinement_criterion);
            pp.query("gamma", value.gamma);
            pp.query("cfl", value.cfl);
            pp.query("rho_solid", value.rho_solid);
            pp.query("rho_fluid", value.rho_fluid);
            pp.query("E_solid", value.E_solid);
            pp.query("E_fluid", value.E_fluid);

            value.bc_eta = new BC::Constant(1, pp, "pf.eta.bc");
            value.bc_rho = new BC::Constant(1, pp, "rho.bc");
            value.bc_E = new BC::Constant(1, pp, "E.bc");
            value.bc_M = new BC::Constant(2, pp, "M.bc");
            value.bc_v = new BC::Constant(2, pp, "v.bc");
            value.bc_p = new BC::Constant(1, pp, "p.bc");

        }
        // Register FabFields:
        {
	  int nghost = 2;
            value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, nghost, "eta", true);
            value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, nghost, "eta_old", false);

            value.RegisterNewFab(value.Density_mf, value.bc_rho, 1, nghost, "Density", true);
            value.RegisterNewFab(value.Density_old_mf, value.bc_rho, 1, nghost, "rho_old", false);

            value.RegisterNewFab(value.Energy_mf, value.bc_E, 1, nghost, "Energy", true);
            value.RegisterNewFab(value.Energy_old_mf, value.bc_E, 1, nghost, "E_old", false);

            value.RegisterNewFab(value.Momentum_mf, value.bc_M, 2, nghost, "Momentum", true);
            value.RegisterNewFab(value.Momentum_old_mf, value.bc_M, 2, nghost, "M_old", false);

            value.RegisterNewFab(value.Velocity_mf, value.bc_v, 2, nghost, "Velocity", true);

            value.RegisterNewFab(value.Pressure_mf, value.bc_p, 1, nghost, "Pressure", true);
        }
        {
            BL_PROFILE("Integrator::Hydro::Hydro()");
            std::string type = "constant";
            pp.query("eta.ic.type", type);
            if (type == "constant") value.ic_eta = new IC::Constant(value.geom, pp, "eta.ic.constant");
            else if (type == "laminate") value.ic_eta = new IC::Laminate(value.geom, pp, "eta.ic.laminate");
            else Util::Abort(INFO, "Invalid eta.ic: ", type);
        }

    }

protected:
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
    //void Integrate(int amrlev, Set::Scalar time, int step, const amrex::MFIter &mfi, const amrex::Box &box) override;
    //virtual void UpdateModel(int a_step);
private:

    Set::Field<Set::Scalar> Density_mf;
    Set::Field<Set::Scalar> Density_old_mf;

    Set::Field<Set::Scalar> Energy_mf;
    Set::Field<Set::Scalar> Energy_old_mf;

    Set::Field<Set::Scalar> Momentum_mf;
    Set::Field<Set::Scalar> Momentum_old_mf;

    Set::Field<Set::Scalar> Velocity_mf;

    Set::Field<Set::Scalar> Pressure_mf;

    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;

    BC::BC<Set::Scalar>* bc_rho;
    BC::BC<Set::Scalar>* bc_E;
    BC::BC<Set::Scalar>* bc_M;
    BC::BC<Set::Scalar>* bc_v;
    BC::BC<Set::Scalar>* bc_p;
    BC::BC<Set::Scalar>* bc_eta;

    IC::IC* ic_rho;
    IC::IC* ic_E;

    IC::IC* ic_M;

    IC::IC* ic_eta;

    Set::Scalar rho_solid = 1.2;
    Set::Scalar rho_fluid = 0.1;
    Set::Scalar E_solid = 1.0;
    Set::Scalar E_fluid = 1.0;

    Set::Scalar Mx_init = 0.0;
    Set::Scalar My_init = 0.0;

    Set::Scalar V_x = 0.0;
    Set::Scalar V_y = 0.0;

    Set::Scalar gamma = 1.4;

    Set::Scalar cfl = 0.6;
    Set::Scalar c_max = 0.0;
    Set::Scalar vx_max = 0.0;
    Set::Scalar vy_max = 0.0;

    Set::Scalar epsilon = 0.1;

    //refinement criterion variables
    Set::Scalar r_refinement_criterion = 0.1;
    Set::Scalar e_refinement_criterion = 0.1;
    Set::Scalar m_refinement_criterion = 0.1;
    Set::Scalar eta_refinement_criterion = 0.11;
};
}

#endif
