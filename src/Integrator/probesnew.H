#pragma once
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

static inline amrex::IntVect PickDeterministicValidPoint(amrex::MultiFab const& mf)
{
    // Deterministic choice: first MFIter, first validbox, offset inside box.
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        amrex::Box vb = mfi.validbox();
        amrex::IntVect iv = vb.smallEnd();

        // Nudge inward to avoid boundary faces/cells
        iv[0] = std::min(iv[0] + 1, vb.bigEnd(0));
        iv[1] = std::min(iv[1] + 1, vb.bigEnd(1));
#if (AMREX_SPACEDIM==3)
        iv[2] = std::min(iv[2] + 0, vb.bigEnd(2));
#endif
        return iv;
    }
    return amrex::IntVect(AMREX_D_DECL(0,0,0)); // should never happen if mf non-empty
}

static inline void ProbeCC_AtValidPoint(
    char const* tag,
    amrex::MultiFab const& mf,
    int lev, int comp)
{
    amrex::Gpu::streamSynchronize();

    amrex::IntVect iv = PickDeterministicValidPoint(mf);

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& vb = mfi.validbox();
        if (!vb.contains(iv)) continue;

        auto const& a = mf.const_array(mfi);
        amrex::Print() << tag
                       << " lev="<<lev
                       << " iv=("<<iv[0]<<","<<iv[1]<<","<<iv[2]<<")"
                       << " vb="<<vb
                       << " comp="<<comp
                       << " val="<<a(iv[0],iv[1],iv[2],comp)
                       << "\n";
        return;
    }

    amrex::Print() << tag << " lev="<<lev << " valid-point probe failed\n";
}

static inline void ProbeFixedCellValid(
    char const* tag,
    amrex::MultiFab const& mf,
    int lev,
    int i, int j, int k,
    int comp)
{
    amrex::Gpu::streamSynchronize();

    amrex::IntVect iv(AMREX_D_DECL(i,j,k));
    bool printed = false;

    for (amrex::MFIter mfi(mf); mfi.isValid() && !printed; ++mfi) {
        const amrex::Box& vb = mfi.validbox();   // VALID ONLY (no ghosts)
        if (!vb.contains(iv)) continue;

        auto const& a = mf.const_array(mfi);

        amrex::AllPrint() << tag
            << " lev="<<lev
            << " iv=("<<i<<","<<j<<","<<k<<")"
            << " vb="<<vb
            << " comp="<<comp
            << " val="<<a(i,j,k,comp)
            << "\n";

        printed = true;
    }

    // Optional: if you want to know when it's NOT valid anywhere on this rank:
    // if (!printed) amrex::AllPrint() << tag << " iv not in validbox on this rank\n";
}

