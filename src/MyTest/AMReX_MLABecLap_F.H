#ifndef AMREX_MLABECLAP_F_H_
#define AMREX_MLABECLAP_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_mlabeclap_flux (const int* lo, const int* hi,
                               amrex_real* fx, const int* fxlo, const int* fxhi,
#if (AMREX_SPACEDIM >= 2)
                               amrex_real* fy, const int* fylo, const int* fyhi,
#if (AMREX_SPACEDIM == 3)
                               amrex_real* fz, const int* fzlo, const int* fzhi,
#endif
#endif
                               const amrex_real* sol, const int* slo, const int* shi,
                               const amrex_real* bx, const int* bxlo, const int* bxhi,
#if (AMREX_SPACEDIM >= 2)
                               const amrex_real* by, const int* bylo, const int* byhi,
#if (AMREX_SPACEDIM == 3)
                               const amrex_real* bz, const int* bzlo, const int* bzhi,
#endif
#endif
                               const amrex_real* dxinv, const amrex_real beta, const int face_only);

#ifdef __cplusplus
}
#endif

#endif
