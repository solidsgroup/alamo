#include "Integrator/Hydro_CVODE.H"
#include "Integrator/Hydro.H"  // full Hydro definition
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spgmr.h> 
#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_context.h>
#include <iostream>

namespace hydro_cvode {

// ------------------------ Pack / Unpack ------------------------

long PackStateToArray(const amrex::MultiFab &rho_mf,
                      const amrex::MultiFab &M_mf,
                      const amrex::MultiFab &E_mf,
                      realtype *array_ptr)
{
    long idx = 0;
    for (amrex::MFIter mfi(rho_mf); mfi.isValid(); ++mfi) {
        const auto &rho_arr = rho_mf[mfi].array();
        const auto &M_arr   = M_mf[mfi].array();
        const auto &E_arr   = E_mf[mfi].array();
        const auto &bx      = mfi.validbox();
        for (int k = 0; k < rho_mf.nComp(); ++k)
        for (int kx = bx.smallEnd(0); kx <= bx.bigEnd(0); ++kx)
        for (int ky = bx.smallEnd(1); ky <= bx.bigEnd(1); ++ky)
        for (int kz = bx.smallEnd(2); kz <= bx.bigEnd(2); ++kz)
            array_ptr[idx++] = rho_arr(kx,ky,kz,k);

        for (int d = 0; d < M_mf.nComp(); ++d)
        for (int kx = bx.smallEnd(0); kx <= bx.bigEnd(0); ++kx)
        for (int ky = bx.smallEnd(1); ky <= bx.bigEnd(1); ++ky)
        for (int kz = bx.smallEnd(2); kz <= bx.bigEnd(2); ++kz)
            array_ptr[idx++] = M_arr(kx,ky,kz,d);

        for (int kx = bx.smallEnd(0); kx <= bx.bigEnd(0); ++kx)
        for (int ky = bx.smallEnd(1); ky <= bx.bigEnd(1); ++ky)
        for (int kz = bx.smallEnd(2); kz <= bx.bigEnd(2); ++kz)
            array_ptr[idx++] = E_arr(kx,ky,kz,0);
    }
    return idx;
}

long UnpackArrayToState(const realtype *array_ptr,
                        amrex::MultiFab &rho_mf,
                        amrex::MultiFab &M_mf,
                        amrex::MultiFab &E_mf)
{
    long idx = 0;
    for (amrex::MFIter mfi(rho_mf); mfi.isValid(); ++mfi) {
        auto rho_arr = rho_mf[mfi].array();
        auto M_arr   = M_mf[mfi].array();
        auto E_arr   = E_mf[mfi].array();
        const auto &bx = mfi.validbox();
        for (int k = 0; k < rho_mf.nComp(); ++k)
        for (int kx = bx.smallEnd(0); kx <= bx.bigEnd(0); ++kx)
        for (int ky = bx.smallEnd(1); ky <= bx.bigEnd(1); ++ky)
        for (int kz = bx.smallEnd(2); kz <= bx.bigEnd(2); ++kz)
            rho_arr(kx,ky,kz,k) = array_ptr[idx++];

        for (int d = 0; d < M_mf.nComp(); ++d)
        for (int kx = bx.smallEnd(0); kx <= bx.bigEnd(0); ++kx)
        for (int ky = bx.smallEnd(1); ky <= bx.bigEnd(1); ++ky)
        for (int kz = bx.smallEnd(2); kz <= bx.bigEnd(2); ++kz)
            M_arr(kx,ky,kz,d) = array_ptr[idx++];

        for (int kx = bx.smallEnd(0); kx <= bx.bigEnd(0); ++kx)
        for (int ky = bx.smallEnd(1); ky <= bx.bigEnd(1); ++ky)
        for (int kz = bx.smallEnd(2); kz <= bx.bigEnd(2); ++kz)
            E_arr(kx,ky,kz,0) = array_ptr[idx++];
    }
    return idx;
}

// ------------------------ CVODE RHS ------------------------

int CVodeRHS(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    CVODEUserData *ud = static_cast<CVODEUserData*>(user_data);
    Integrator::Hydro *hydro = ud->hydro;

    // Temporary MultiFabs for state
    amrex::MultiFab rho_tmp(ud->rho_mf->boxArray(), ud->rho_mf->DistributionMap(),
                            ud->rho_mf->nComp(), ud->rho_mf->nGrow());
    amrex::MultiFab M_tmp(ud->M_mf->boxArray(), ud->M_mf->DistributionMap(),
                          ud->M_mf->nComp(), ud->M_mf->nGrow());
    amrex::MultiFab E_tmp(ud->E_mf->boxArray(), ud->E_mf->DistributionMap(),
                          ud->E_mf->nComp(), ud->E_mf->nGrow());

    // Unpack N_Vector into MultiFabs
    UnpackArrayToState(NV_DATA_S(y), rho_tmp, M_tmp, E_tmp);

    // Call Hydro RHS via public wrapper
    hydro->RHSPublic(0, static_cast<Set::Scalar>(t),
                     rho_tmp, M_tmp, E_tmp,
                     *ud->rho_mf, *ud->M_mf, *ud->E_mf);

    // Pack back to ydot
    PackStateToArray(rho_tmp, M_tmp, E_tmp, NV_DATA_S(ydot));

    return 0;
}

// Jacobian-times-vector function for CVODE SPGMR
int CVodeJtv(N_Vector v, N_Vector Jv,
             realtype t, N_Vector y, N_Vector fy,
             void *user_data, N_Vector tmp)
{
    CVODEUserData *ud = static_cast<CVODEUserData*>(user_data);
    auto *hydro = ud->hydro;
    const realtype sigma = 1e-8;

    N_Vector y_tmp = N_VClone(y);
    N_Vector f_yv  = N_VClone(y);

    // f(y + sigma*v)
    N_VLinearSum(1.0, y, sigma, v, y_tmp);
    CVodeRHS(t, y_tmp, f_yv, user_data);

    // Jv ≈ (f(y + σv) - f(y)) / σ
    N_VLinearSum(1.0/sigma, f_yv, -1.0/sigma, fy, Jv);

    N_VDestroy(y_tmp);
    N_VDestroy(f_yv);

    return 0;
}

// ------------------------ AdvanceImplicit ------------------------

int AdvanceImplicit(Integrator::Hydro *hydro,
                    amrex::MultiFab &rho_mf,
                    amrex::MultiFab &M_mf,
                    amrex::MultiFab &E_mf,
                    int lev,
                    double t0,
                    double dt,
                    const CVODEConfig &cfg)
{
    long N_local = rho_mf.nComp() * rho_mf.boxArray().numPts()
                 + M_mf.nComp()   * M_mf.boxArray().numPts()
                 + E_mf.nComp()   * E_mf.boxArray().numPts();

    SUNContext sunctx;
    SUNContext_Create(nullptr, &sunctx);   // create context
    
    N_Vector y    = N_VNew_Serial(N_local, sunctx);
    N_Vector ydot = N_VNew_Serial(N_local, sunctx);

    // Pack initial state into y
    PackStateToArray(rho_mf, M_mf, E_mf, NV_DATA_S(y));

    void *cvode_mem = CVodeCreate(CV_BDF, sunctx);
    int flag = CVodeInit(cvode_mem, CVodeRHS, t0, y);
    if (flag != 0) {
        std::cerr << "CVodeInit failed with flag=" << flag << "\n";
        return flag;
    }
    CVodeSetInitStep(cvode_mem, 1e-12);
    CVodeSetMaxStep(cvode_mem, 0.01*dt);
    flag = CVodeSStolerances(cvode_mem, cfg.rel_tol, cfg.abs_tol);
    if (flag != CV_SUCCESS) {
        amrex::Abort("CVodeSStolerances failed!");
    }

    SUNMatrix A = SUNDenseMatrix(N_local, N_local, sunctx);
    if (A == nullptr) amrex::Abort("SUNDenseMatrix creation failed");
    //SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    //if (LS == nullptr) amrex::Abort("SUNLinSol_Dense creation failed");
    int pretype = PREC_LEFT;
    int maxl = 50;
    SUNLinearSolver LS = SUNLinSol_SPGMR(y, PREC_LEFT, maxl, sunctx);
    if (LS == nullptr) amrex::Abort("SUNLinSol_SPGMR creation failed");
    
    // Attach the linear solver to CVODE
    //flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    //if (flag != CV_SUCCESS) {
    //    amrex::Abort("CVodeSetLinearSolver failed");
    //}
    flag = CVodeSetLinearSolver(cvode_mem, LS, nullptr);
    if (flag != CV_SUCCESS) {
        amrex::Abort("CVodeSetLinearSolver failed");
    }

    // Attach the Jacobian-times-vector function
    flag = CVodeSetJacTimes(cvode_mem, nullptr, CVodeJtv);
    if (flag != CV_SUCCESS) amrex::Abort("CVodeSetJacTimes failed");
    
    // (Optional) Preconditioner interface stubs
    // If you set PREC_LEFT or PREC_RIGHT above, define these:
    auto prec_setup = [](realtype t, N_Vector y, N_Vector fy,
                         booleantype jok, booleantype *jcurPtr,
                         realtype gamma, void *user_data) {
        *jcurPtr = SUNFALSE;
        return 0; // success
    };
    auto prec_solve = [](realtype t, N_Vector y, N_Vector fy,
                         N_Vector r, N_Vector z,
                         realtype gamma, realtype delta,
                         int lr, void *user_data) {
        N_VScale(1.0, r, z);
        return 0; // identity preconditioner
    };
    
    // If using PREC_LEFT or PREC_RIGHT, enable:
    if (pretype != PREC_NONE) {
        flag = CVodeSetPreconditioner(cvode_mem, prec_setup, prec_solve);
        if (flag != CV_SUCCESS) amrex::Abort("CVodeSetPreconditioner failed");
    }

    // Set user data
    CVODEUserData ud;
    ud.hydro = hydro;
    ud.rho_mf = &rho_mf;
    ud.M_mf = &M_mf;
    ud.E_mf = &E_mf;
    CVodeSetUserData(cvode_mem, &ud);
    
    // Integrate one timestep
    double t_out = t0 + dt;
    flag = CVode(cvode_mem, t_out, y, &t_out, CV_NORMAL);
    if (flag < 0) {
        std::cerr << "CVode integration failed with flag=" << flag << "\n";
    }
    long nsteps;
    flag = CVodeGetNumSteps(cvode_mem, &nsteps);
    if (flag != CV_SUCCESS) {
        amrex::Abort("CVodeGetNumSteps failed");
    }
    amrex::Print() << "CVODE internal steps taken: " << nsteps << "\n";
    long nfevals;
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfevals);
    if (flag != CV_SUCCESS) {
        amrex::Abort("CVodeGetNumRhsEvals failed");
    }
    amrex::Print() << "CVODE RHS evaluations: " << nfevals << "\n";

    // Unpack result back into MultiFabs
    UnpackArrayToState(NV_DATA_S(y), rho_mf, M_mf, E_mf);

    // Free memory
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(y);
    N_VDestroy(ydot);
    SUNContext_Free(&sunctx);
    CVodeFree(&cvode_mem);

    return flag;
}

} // namespace hydro_cvode

