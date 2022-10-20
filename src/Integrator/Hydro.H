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
    void Advance (int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
                    const amrex::MFIter &mfi, const amrex::Box &box) override;
    void UpdateModel(int a_step) override;
private:

    Set::Field<Set::Scalar> rho_mf;
    Set::Field<Set::Scalar> rho_old_mf;
    Set::Field<Set::Scalar> E_mf;
    Set::Field<Set::Scalar> E_old_mf;
  
    Set::Field<Set::Scalar> P_mf;
    Set::Field<Set::Scalar> P_old_mf;
  
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;

    Set::Field<Set::Scalar> deta_mf;
  
    BC::BC<Set::Scalar> *bc_rho;
    BC::BC<Set::Scalar> *bc_E;
  
    BC::BC<Set::Scalar> *bc_P;
  
    BC::BC<Set::Scalar> *bc_eta;

    IC::IC *ic_rho;
    IC::IC *ic_E;
  
    IC::IC *ic_P;
  
    IC::IC *ic_eta;
    
 
    Set::Scalar gamma = 1.4;
    Set::Scalar cfl = 0.8;
};
}

#endif
