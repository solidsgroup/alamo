#include "Expression.H"

namespace BC
{

void
Expression::FillBoundary (amrex::BaseFab<Set::Scalar> &a_in,
                        const amrex::Box &a_box,
                        int ngrow, int /*dcomp*/, int /*ncomp*/, Set::Scalar time,
                        Orientation face, const amrex::Mask * /*mask*/)
{
    const auto DX = m_geom.CellSizeArray();
    const auto prob_lo = m_geom.ProbLoArray();

    Util::Assert(INFO,TEST(a_in.nComp() == (int)m_ncomp));

    amrex::Box box = a_box;
    box.grow(ngrow);
    const amrex::Dim3 lo= amrex::lbound(m_geom.Domain()), hi = amrex::ubound(m_geom.Domain());

    amrex::Array4<amrex::Real> const& in = a_in.array();

    amrex::IndexType type = amrex::IndexType::TheCellType();

    for (int n = 0; n < a_in.nComp(); n++)
    {
        const auto bc_type_xlo = m_bc_type[Face::XLO][n];
        const auto bc_type_xhi = m_bc_type[Face::XHI][n];
        const auto bc_type_ylo = m_bc_type[Face::YLO][n];
        const auto bc_type_yhi = m_bc_type[Face::YHI][n];
        const auto bc_func_xlo = m_bc_func[Face::XLO][n];
        const auto bc_func_xhi = m_bc_func[Face::XHI][n];
        const auto bc_func_ylo = m_bc_func[Face::YLO][n];
        const auto bc_func_yhi = m_bc_func[Face::YHI][n];
#if AMREX_SPACEDIM > 2
        const auto bc_type_zlo = m_bc_type[Face::ZLO][n];
        const auto bc_type_zhi = m_bc_type[Face::ZHI][n];
        const auto bc_func_zlo = m_bc_func[Face::ZLO][n];
        const auto bc_func_zhi = m_bc_func[Face::ZHI][n];
#endif
        amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Vector pos = Set::Position(i, j, k, prob_lo, DX, type);
            Set::Scalar x = 0.0, y=0.0, z=0.0, t=time;
            x = pos(0);
            #if AMREX_SPACEDIM > 1
            y = pos(1);
            #if AMREX_SPACEDIM > 2
            z = pos(2);
            #endif
            #endif


            amrex::IntVect glevel;
            AMREX_D_TERM(glevel[0] = std::max(std::min(0,i-lo.x),i-hi.x); ,
                        glevel[1] = std::max(std::min(0,j-lo.y),j-hi.y); ,
                        glevel[2] = std::max(std::min(0,k-lo.z),k-hi.z); );

            if (glevel[0]<0 && (face == Orientation::xlo || face == Orientation::All)) // Left boundary
            {
                if (BCUtil::IsDirichlet(bc_type_xlo))
                    in(i,j,k,n) = bc_func_xlo(x,y,z,t);
                else if(BCUtil::IsNeumann(bc_type_xlo))
                    in(i,j,k,n) = in(i-glevel[0],j,k,n) - bc_func_xlo(x,y,z,t)*DX[0];
                else if(BCUtil::IsReflectEven(bc_type_xlo))
                    in(i,j,k,n) = in(1-glevel[0],j,k,n);
                else if(BCUtil::IsReflectOdd(bc_type_xlo))
                    in(i,j,k,n) = -in(1-glevel[0],j,k,n);
                else if(BCUtil::IsPeriodic(bc_type_xlo)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if (glevel[0]>0 && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
            {
                if (BCUtil::IsDirichlet(bc_type_xhi))
                    in(i,j,k,n) = bc_func_xhi(x,y,z,t);
                else if(BCUtil::IsNeumann(bc_type_xhi))
                    in(i,j,k,n) = in(i-glevel[0],j,k,n) - bc_func_xhi(x,y,z,t)*DX[0];
                else if(BCUtil::IsReflectEven(bc_type_xhi))
                    in(i,j,k,n) = in(hi.x-glevel[0],j,k,n);
                else if(BCUtil::IsReflectOdd(bc_type_xhi))
                    in(i,j,k,n) = -in(hi.x-glevel[0],j,k,n);
                else if(BCUtil::IsPeriodic(bc_type_xhi)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if (glevel[1]<0 && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
            {
                if (BCUtil::IsDirichlet(bc_type_ylo))
                    in(i,j,k,n) = bc_func_ylo(x,y,z,t);
                else if (BCUtil::IsNeumann(bc_type_ylo))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n) - bc_func_ylo(x,y,z,t)*DX[1];
                else if (BCUtil::IsReflectEven(bc_type_ylo))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n);
                else if (BCUtil::IsReflectOdd(bc_type_ylo))
                    in(i,j,k,n) = -in(i,j-glevel[1],k,n);
                else if(BCUtil::IsPeriodic(bc_type_ylo)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if (glevel[1]>0 && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
            {
                if (BCUtil::IsDirichlet(bc_type_yhi))
                    in(i,j,k,n) = bc_func_yhi(x,y,z,t);
                else if (BCUtil::IsNeumann(bc_type_yhi))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n) - bc_func_yhi(x,y,z,t)*DX[1];
                else if (BCUtil::IsReflectEven(bc_type_yhi))
                    in(i,j,k,n) = in(i,hi.y-glevel[1],k,n);
                else if (BCUtil::IsReflectOdd(bc_type_yhi))
                    in(i,j,k,n) = -in(i,hi.y-glevel[1],k,n);
                else if(BCUtil::IsPeriodic(bc_type_yhi)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
#if AMREX_SPACEDIM>2
            else if (glevel[2]<0 && (face == Orientation::zlo || face == Orientation::All))
            {
                if (BCUtil::IsDirichlet(bc_type_zlo))
                    in(i,j,k,n) = bc_func_zlo(x,y,z,t);
                else if (BCUtil::IsNeumann(bc_type_zlo))
                    in(i,j,k,n) = in(i,j,k-glevel[2],n) - bc_func_zlo(x,y,z,t)*DX[2];
                else if (BCUtil::IsReflectEven(bc_type_zlo))
                    in(i,j,k,n) = in(i,j,1-glevel[2],n);
                else if (BCUtil::IsReflectOdd(bc_type_zlo))
                    in(i,j,k,n) = -in(i,j,1-glevel[2],n);
                else if(BCUtil::IsPeriodic(bc_type_zlo)) {}
                else Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if (glevel[2]>0 && (face == Orientation::zhi || face == Orientation::All))
            {
                if (BCUtil::IsDirichlet(bc_type_zhi))
                    in(i,j,k,n) = bc_func_zhi(x,y,z,t);
                else if(BCUtil::IsNeumann(bc_type_zhi))
                    in(i,j,k,n) = in(i,j,k-glevel[2],n) - bc_func_zhi(x,y,z,t)*DX[2];
                else if(BCUtil::IsReflectEven(bc_type_zhi))
                    in(i,j,k,n) = in(i,j,hi.z-glevel[2],n);
                else if(BCUtil::IsReflectOdd(bc_type_zhi))
                    in(i,j,k,n) = -in(i,j,hi.z-glevel[2],n);
                else if(BCUtil::IsPeriodic(bc_type_zhi)) {}
                else Util::Abort(INFO, "Incorrect boundary conditions");
            }
#endif

        });

        // Average physical corner ghost cells from the face ghosts, matching
        // Constant BC behavior and avoiding undefined corner values in centered stencils.
        amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            if (i < lo.x && j < lo.y) in(i,j,k,n) = 0.5*( in(i+1,j,k,n) + in(i,j+1,k,n) );
            if (i < lo.x && j > hi.y) in(i,j,k,n) = 0.5*( in(i+1,j,k,n) + in(i,j-1,k,n) );
            if (i > hi.x && j < lo.y) in(i,j,k,n) = 0.5*( in(i-1,j,k,n) + in(i,j+1,k,n) );
            if (i > hi.x && j > hi.y) in(i,j,k,n) = 0.5*( in(i-1,j,k,n) + in(i,j-1,k,n) );
        });
    }
}

amrex::BCRec
Expression::GetBCRec() 
{
    int bc_lo[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XLO][0],m_bc_type[Face::YLO][0],m_bc_type[Face::XLO][0])};
    int bc_hi[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XHI][0],m_bc_type[Face::YHI][0],m_bc_type[Face::XHI][0])};

    return amrex::BCRec(bc_lo,bc_hi);
}

amrex::Array<int,AMREX_SPACEDIM>
Expression::IsPeriodic()
{
    return {AMREX_D_DECL(BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))};
}
amrex::Periodicity Expression::Periodicity () const
{
    return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(m_geom.Domain().length(0) * BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                                                            m_geom.Domain().length(1) * BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                                                            m_geom.Domain().length(2) * BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))));
}
amrex::Periodicity Expression::Periodicity (const amrex::Box& b) {
    return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(b.length(0) * BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                                                        b.length(1) * BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                                                        b.length(2) * BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))));

}


}
