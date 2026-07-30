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
    const auto prob_hi = m_geom.ProbHiArray();

    Util::Assert(INFO,TEST(a_in.nComp() == (int)m_ncomp));

    amrex::Box box = a_box;
    box.grow(ngrow);
    amrex::Box domain = m_geom.Domain();
    const amrex::IndexType type = a_in.box().ixType();
    domain.convert(type);
    const bool nodal = type == amrex::IndexType::TheNodeType();
    const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);

    amrex::Array4<amrex::Real> const& in = a_in.array();

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
            Set::Vector pos =
                Set::Position(i, j, k, prob_lo, DX, type);
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

            if ((glevel[0]<0 || (nodal && i == lo.x)) && (face == Orientation::xlo || face == Orientation::All)) // Left boundary
            {
                if (BCUtil::IsDirichlet(bc_type_xlo))
                {
                    const Set::Scalar val = bc_func_xlo(prob_lo[0],y,z,t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[0]<0) in(i,j,k,n) = 2.0 * val - in(2 * lo.x - i - 1,j,k,n);
                }
                else if(glevel[0]<0 && BCUtil::IsNeumann(bc_type_xlo))
                    in(i,j,k,n) = in(i-glevel[0],j,k,n) - bc_func_xlo(x,y,z,t)*DX[0];
                else if(glevel[0]<0 && BCUtil::IsReflectEven(bc_type_xlo))
                    in(i,j,k,n) = in(1-glevel[0],j,k,n);
                else if(glevel[0]<0 && BCUtil::IsReflectOdd(bc_type_xlo))
                    in(i,j,k,n) = -in(1-glevel[0],j,k,n);
                else if(BCUtil::IsPeriodic(bc_type_xlo) || (nodal && i == lo.x)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[0]>0 || (nodal && i == hi.x)) && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
            {
                if (BCUtil::IsDirichlet(bc_type_xhi))
                {
                    const Set::Scalar val = bc_func_xhi(prob_hi[0],y,z,t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[0]>0) in(i,j,k,n) = 2.0 * val - in(2 * hi.x - i + 1,j,k,n);
                }
                else if(glevel[0]>0 && BCUtil::IsNeumann(bc_type_xhi))
                    in(i,j,k,n) = in(i-glevel[0],j,k,n) - bc_func_xhi(x,y,z,t)*DX[0];
                else if(glevel[0]>0 && BCUtil::IsReflectEven(bc_type_xhi))
                    in(i,j,k,n) = in(hi.x-glevel[0],j,k,n);
                else if(glevel[0]>0 && BCUtil::IsReflectOdd(bc_type_xhi))
                    in(i,j,k,n) = -in(hi.x-glevel[0],j,k,n);
                else if(BCUtil::IsPeriodic(bc_type_xhi) || (nodal && i == hi.x)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[1]<0 || (nodal && j == lo.y)) && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
            {
                if (BCUtil::IsDirichlet(bc_type_ylo))
                {
                    const Set::Scalar val = bc_func_ylo(x,prob_lo[1],z,t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[1]<0) in(i,j,k,n) = 2.0 * val - in(i,2 * lo.y - j - 1,k,n);
                }
                else if (glevel[1]<0 && BCUtil::IsNeumann(bc_type_ylo))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n) - bc_func_ylo(x,y,z,t)*DX[1];
                else if (glevel[1]<0 && BCUtil::IsReflectEven(bc_type_ylo))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n);
                else if (glevel[1]<0 && BCUtil::IsReflectOdd(bc_type_ylo))
                    in(i,j,k,n) = -in(i,j-glevel[1],k,n);
                else if(BCUtil::IsPeriodic(bc_type_ylo) || (nodal && j == lo.y)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[1]>0 || (nodal && j == hi.y)) && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
            {
                if (BCUtil::IsDirichlet(bc_type_yhi))
                {
                    const Set::Scalar val = bc_func_yhi(x,prob_hi[1],z,t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[1]>0) in(i,j,k,n) = 2.0 * val - in(i,2 * hi.y - j + 1,k,n);
                }
                else if (glevel[1]>0 && BCUtil::IsNeumann(bc_type_yhi))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n) - bc_func_yhi(x,y,z,t)*DX[1];
                else if (glevel[1]>0 && BCUtil::IsReflectEven(bc_type_yhi))
                    in(i,j,k,n) = in(i,hi.y-glevel[1],k,n);
                else if (glevel[1]>0 && BCUtil::IsReflectOdd(bc_type_yhi))
                    in(i,j,k,n) = -in(i,hi.y-glevel[1],k,n);
                else if(BCUtil::IsPeriodic(bc_type_yhi) || (nodal && j == hi.y)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
#if AMREX_SPACEDIM>2
            else if ((glevel[2]<0 || (nodal && k == lo.z)) && (face == Orientation::zlo || face == Orientation::All))
            {
                if (BCUtil::IsDirichlet(bc_type_zlo))
                {
                    const Set::Scalar val = bc_func_zlo(x,y,prob_lo[2],t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[2]<0) in(i,j,k,n) = 2.0 * val - in(i,j,2 * lo.z - k - 1,n);
                }
                else if (glevel[2]<0 && BCUtil::IsNeumann(bc_type_zlo))
                    in(i,j,k,n) = in(i,j,k-glevel[2],n) - bc_func_zlo(x,y,z,t)*DX[2];
                else if (glevel[2]<0 && BCUtil::IsReflectEven(bc_type_zlo))
                    in(i,j,k,n) = in(i,j,1-glevel[2],n);
                else if (glevel[2]<0 && BCUtil::IsReflectOdd(bc_type_zlo))
                    in(i,j,k,n) = -in(i,j,1-glevel[2],n);
                else if(BCUtil::IsPeriodic(bc_type_zlo) || (nodal && k == lo.z)) {}
                else Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[2]>0 || (nodal && k == hi.z)) && (face == Orientation::zhi || face == Orientation::All))
            {
                if (BCUtil::IsDirichlet(bc_type_zhi))
                {
                    const Set::Scalar val = bc_func_zhi(x,y,prob_hi[2],t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[2]>0) in(i,j,k,n) = 2.0 * val - in(i,j,2 * hi.z - k + 1,n);
                }
                else if(glevel[2]>0 && BCUtil::IsNeumann(bc_type_zhi))
                    in(i,j,k,n) = in(i,j,k-glevel[2],n) - bc_func_zhi(x,y,z,t)*DX[2];
                else if(glevel[2]>0 && BCUtil::IsReflectEven(bc_type_zhi))
                    in(i,j,k,n) = in(i,j,hi.z-glevel[2],n);
                else if(glevel[2]>0 && BCUtil::IsReflectOdd(bc_type_zhi))
                    in(i,j,k,n) = -in(i,j,hi.z-glevel[2],n);
                else if(BCUtil::IsPeriodic(bc_type_zhi) || (nodal && k == hi.z)) {}
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
Expression::GetBCRec(int component)
{
    int bc_lo[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XLO][component],m_bc_type[Face::YLO][component],m_bc_type[Face::ZLO][component])};
    int bc_hi[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XHI][component],m_bc_type[Face::YHI][component],m_bc_type[Face::ZHI][component])};

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
