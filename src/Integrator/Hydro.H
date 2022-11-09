#ifndef INTEGRATOR_HYDRO_H
#define INTEGRATOR_HYDRO_H

#include <string>
#include <limits>
#include <memory>

#include "Integrator/Integrator.H"
#include "IO/ParmParse.H"

#include "BC/BC.H"
#include "IC/IC.H"

#include "IC/Cylinder.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/Expression.H"

namespace Integrator
{
class Hydro : virtual public Integrator::Integrator
{
public:
  Hydro() {};
    Hydro(IO::ParmParse &pp) : Hydro() {pp.queryclass(*this);}

    static void Parse(Hydro &value, IO::ParmParse &pp)
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
	//pp.query("rho0", value.rho0);
	//pp.query("V", value.V);

      }
      // Register FabFields:
      {
	value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, 2, "eta", true);
	value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 2, "eta_old", false);

	value.RegisterNewFab(value.Density_mf, value.bc_rho, 1, 2, "Density", true);
	value.RegisterNewFab(value.Density_old_mf, value.bc_rho, 1, 2, "rho_old", false);

	value.RegisterNewFab(value.Energy_mf, value.bc_E, 1, 2, "Energy", true);
	value.RegisterNewFab(value.Energy_old_mf, value.bc_E, 1, 2, "E_old", false);

	value.RegisterNewFab(value.Momentum_mf, value.bc_M, 1, 2, "Momentum", true);
	value.RegisterNewFab(value.Momentum_old_mf, value.bc_M, 1, 2, "M_old", false);

	value.RegisterNewFab(value.Velocity_mf, value.bc_M, 1, 2, "Velocity", true);
	value.RegisterNewFab(value.Pressure_mf, value.bc_M, 1, 2, "Pressure", true);

      }
    }

protected:
    void Initialize (int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void TimeStepEnd(int lev);
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
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

    Set::Field<Set::Scalar> flux_x_mf;
    Set::Field<Set::Scalar> flux_y_mf;

    BC::BC<Set::Scalar> *bc_rho;
    BC::BC<Set::Scalar> *bc_E;
    BC::BC<Set::Scalar> *bc_M;
  
    BC::BC<Set::Scalar> *bc_eta;

    IC::IC *ic_rho;
    IC::IC *ic_E;
  
    IC::IC *ic_M;
  
    IC::IC *ic_eta;

    Set::Scalar rho_solid = 1.2;
    Set::Scalar rho_fluid = 1.0;
    Set::Scalar E_solid = 5.0;
    Set::Scalar E_fluid = 0.25;
 
    Set::Scalar gamma = 1.4;
    Set::Scalar rho0 = 0.5;
    Set::Scalar V = 1.0;
    Set::Scalar cfl = 0.8;
    Set::Scalar c_max = 0.0;
    Set::Scalar vx_max = 0.0;
    Set::Scalar vy_max = 0.0;
 
    //refinement criterion variables
    Set::Scalar r_refinement_criterion = 0.1;
    Set::Scalar e_refinement_criterion = 0.1;
    Set::Scalar m_refinement_criterion = 0.1;
};
}

#endif
