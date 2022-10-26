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
    Hydro();
    Hydro(IO::ParmParse &pp);

    static void Parse(Hydro &value, IO::ParmParse &pp);

protected:
    void Initialize (int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void TimeStepEnd(Set::Scalar a_time, int a_iter) override;
    void Advance (int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
                    const amrex::MFIter &mfi, const amrex::Box &box) override;
    void UpdateModel(int a_step) override;
private:

    Set::Field<Set::Scalar> Density_mf;
    Set::Field<Set::Scalar> Density_old_mf;
    Set::Field<Set::Scalar> Energy_mf;
    Set::Field<Set::Scalar> Energy_old_mf;
  
    Set::Field<Set::Vector> Momentum_mf;
    Set::Field<Set::Vector> Momentum_old_mf;

    Set::Field<Set::Vector> Velocity_mf;
    Set::Field<Set::Scalar> Pressure_mf;
  
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;

    Set::Field<Set::Scalar> flux_x_mf;
    Set::Field<Set::Scalar> flux_y_mf;

    BC::BC<Set::Scalar> *bc_rho;
    BC::BC<Set::Scalar> *bc_E;
  
    BC::BC<Set::Scalar> *bc_P;
  
    BC::BC<Set::Scalar> *bc_eta;

    IC::IC *ic_rho;
    IC::IC *ic_E;
  
    IC::IC *ic_M;
  
    IC::IC *ic_eta;
 
    Set::Scalar gamma = 1.4;
    Set::Scalar cfl = 0.8;
    Set::Scalar c_max = 0.0;
    Set::Scalar vx_max = 0.0;
    Set::Scalar vy_max = 0.0;
 
    //refinement criterion variables
    Set::Scalar r_refinement_criterion = 0.1;
    Set::Scalar e_refinement_criterion = 0.1;
    Set::Scalar m_refinemetn_criterion = 0.1;
};
}

#endif
