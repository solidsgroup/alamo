#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include "Set/Set.H"

#include <AMReX_ArrayLim.H>

#include "Set/Set.H"
#include "Diagonal.H"



namespace Operator
{
Diagonal::Diagonal (const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info)
{define(a_geom, a_grids, a_dmap, a_info);}
Diagonal::~Diagonal () {}


void
Diagonal::Fapply (int amrlev,
        int mglev,
        MultiFab& f,
        const MultiFab& u
        ) const
{
    BL_PROFILE("Operator::Diagonal::Diagonal::Fapply()");
    std::cout << "in fapply, amrlev = " << amrlev << std::endl; 

    amrex::Box domain(m_geom[amrlev][mglev].Domain());

    for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
    {
        const amrex::FArrayBox &ufab    = u[mfi];
        if(ufab.contains_inf()) Util::Abort(INFO, "Inf in ufab [before update]");
        if(ufab.contains_nan()) Util::Abort(INFO, "Nan in ufab [before update]");
    }

    for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const amrex::FArrayBox &ufab    = u[mfi];
        amrex::FArrayBox       &ffab    = f[mfi];

        AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
                for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
                for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
        {
            amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));
            for (int i=0; i<AMREX_SPACEDIM; i++)
            {
                if (AMREX_D_TERM(m1 == domain.loVect()[0],
                        || m2 == domain.loVect()[1],
                        || m3 == domain.loVect()[2]) ||
                    AMREX_D_TERM(m1 == domain.hiVect()[0]+1,
                        || m2 == domain.hiVect()[1]+1,
                        || m3 == domain.hiVect()[2]+1))
                {
                    ffab(m,i) = 2*ufab(m,i);
                }
                else if (AMREX_D_TERM(m1 == bx.loVect()[0], 
                            || m2 == bx.loVect()[1],
                            || m3 == bx.loVect()[2]) || 
                    AMREX_D_TERM(m1 == bx.hiVect()[0],
                            || m2 == bx.hiVect()[1],
                            || m3 == bx.hiVect()[2]))
                    continue;
                else 
                    ffab(m,i) = 2*ufab(m,i);
            }
        }
    }
    for (MFIter mfi(f, true); mfi.isValid(); ++mfi)
    {
        amrex::FArrayBox       &ffab    = f[mfi];
        if(ffab.contains_inf()) Util::Abort(INFO, "Inf in ffab [after update]");
        if(ffab.contains_nan()) Util::Abort(INFO, "Nan in ffab [after update]");
    }
}


void
Diagonal::Fsmooth (int amrlev,
        int mglev,  
        MultiFab& u, 
        const MultiFab& rhs
        ) const
{
    for (int redblack = 0; redblack < 2; redblack++)
    {
        BL_PROFILE("Operator::Diagonal::Diagonal::Fsmooth()");

        amrex::Box domain(m_geom[amrlev][mglev].Domain());

        for (MFIter mfi(u,MFItInfo().EnableTiling().SetDynamic(true));
            mfi.isValid(); ++mfi)
        {
            const Box&       bx     = mfi.tilebox();
            FArrayBox&       ufab    = u[mfi];
            const FArrayBox& rhsfab  = rhs[mfi];

            AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++),
                    for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++),
                    for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++))
            {
                if ((AMREX_D_TERM(m1, + m2, + m3))%2 == redblack) continue;
                amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));

                for (int i=0; i<AMREX_SPACEDIM; i++)
                {
                    for (int k=0; k<AMREX_SPACEDIM; k++)
                    {
                        if (AMREX_D_TERM(m1 == domain.loVect()[0],
                                || m2 == domain.loVect()[1],
                                || m3 == domain.loVect()[2]) ||
                            AMREX_D_TERM(m1 == domain.hiVect()[0]+1,
                                || m2 == domain.hiVect()[1]+1,
                                || m3 == domain.hiVect()[2]+1))
                        {
                            ufab(m,k) = rhsfab(m,k) / 2.0;
                        }
                        else if (AMREX_D_TERM(m1 == bx.loVect()[0],
                                    || m2 == bx.loVect()[1],
                                    || m3 == bx.loVect()[2]) || 
                            AMREX_D_TERM(m1 == bx.hiVect()[0], 
                                    || m2 == bx.hiVect()[1],
                                    || m3 == bx.hiVect()[2]))
                            continue;
                        else 
                            ufab(m,k) = rhsfab(m,k) / 2.0;
                    }
                }
            }
        }
    }
}
}
