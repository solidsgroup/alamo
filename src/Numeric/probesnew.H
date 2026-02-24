#pragma once
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

static inline amrex::IntVect PickDeterministicValidPoint(amrex::MultiFab const& mf)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        amrex::Box vb = mfi.validbox();
        amrex::IntVect iv = vb.smallEnd();
        iv[0] = std::min(iv[0] + 1, vb.bigEnd(0));
        iv[1] = std::min(iv[1] + 1, vb.bigEnd(1));
#if (AMREX_SPACEDIM==3)
        iv[2] = std::min(iv[2] + 0, vb.bigEnd(2));
#endif
        return iv;
    }
    return amrex::IntVect(AMREX_D_DECL(0,0,0));
}

static inline void ProbeFace_AtValidPoint(
    char const* tag,
    amrex::MultiFab const& mf_face,
    int lev, int comp)
{
    amrex::Gpu::streamSynchronize();
    amrex::IntVect iv = PickDeterministicValidPoint(mf_face);

    for (amrex::MFIter mfi(mf_face); mfi.isValid(); ++mfi) {
        const amrex::Box& vb = mfi.validbox();
        if (!vb.contains(iv)) continue;

        auto const& a = mf_face.const_array(mfi);
        amrex::Print() << tag
                       << " lev="<<lev
                       << " face=("<<iv[0]<<","<<iv[1]<<","<<iv[2]<<")"
                       << " vb="<<vb
                       << " comp="<<comp
                       << " val="<<a(iv[0],iv[1],iv[2],comp)
                       << "\n";
        return;
    }
}
