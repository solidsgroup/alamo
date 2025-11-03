#ifndef INTEGRATOR_CVODE_H
#define INTEGRATOR_CVODE_H

#include "AMReX_MultiFab.H"
#include <nvector/nvector_serial.h>
#include <cvode/cvode.h>
#include <vector>

namespace Integrator {
class Hydro; // forward declaration
}

namespace hydro_cvode {

struct CVODEUserData {
    Integrator::Hydro *hydro;
    amrex::MultiFab *rho_mf;
    amrex::MultiFab *M_mf;
    amrex::MultiFab *E_mf;
    void *cvode_mem;
};

// Packing MultiFabs into a contiguous array (N_Vector)
long PackStateToArray(const amrex::MultiFab &rho_mf,
                      const amrex::MultiFab &momentum_mf,
                      const amrex::MultiFab &energy_mf,
                      realtype *array_ptr);

// Unpacking from N_Vector back into MultiFabs
long UnpackArrayToState(const realtype *array_ptr,
                        amrex::MultiFab &rho_mf,
                        amrex::MultiFab &momentum_mf,
                        amrex::MultiFab &energy_mf);

// CVODE RHS function
int CVodeRHS(realtype t, N_Vector y, N_Vector ydot, void *user_data);

// Configuration structure (optional, can extend)
struct CVODEConfig {
    int max_steps;
    double rel_tol;
    double abs_tol;
    double init_step;
    double max_step;
    double max_num_steps;
};

// Advance function
int AdvanceImplicit(Integrator::Hydro *hydro,
                    amrex::MultiFab &rho_mf,
                    amrex::MultiFab &M_mf,
                    amrex::MultiFab &E_mf,
                    int lev,
                    double t0,
                    double dt,
                    const CVODEConfig &cfg);

} // namespace hydro_cvode

#endif
