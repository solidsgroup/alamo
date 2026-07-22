#include "Expression.H"

namespace BC
{

void
Expression::FillBoundary (amrex::BaseFab<Set::Scalar> &a_in,
                        const amrex::Box &a_box,
                        int ngrow, int /*dcomp*/, int /*ncomp*/, Set::Scalar time,
                        Orientation face, const amrex::Mask * /*mask*/)
{
    const amrex::Real* DX = m_geom.CellSize();

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
        amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Vector pos = Set::Position(i, j, k, m_geom, type);
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
                if (BCUtil::IsDirichlet(m_bc_type[Face::XLO][n]))
                {
                    const Set::Scalar val = m_bc_func[Face::XLO][n](m_geom.ProbLo()[0],y,z,t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[0]<0) in(i,j,k,n) = 2.0 * val - in(2 * lo.x - i - 1,j,k,n);
                }
                else if(glevel[0]<0 && BCUtil::IsNeumann(m_bc_type[Face::XLO][n]))
                    in(i,j,k,n) = in(i-glevel[0],j,k,n) - (m_bc_func[Face::XLO].size() > 0 ? m_bc_func[Face::XLO][n](x,y,z,t)*DX[0] : 0);
                else if(glevel[0]<0 && BCUtil::IsReflectEven(m_bc_type[Face::XLO][n]))
                    in(i,j,k,n) = in(1-glevel[0],j,k,n);
                else if(glevel[0]<0 && BCUtil::IsReflectOdd(m_bc_type[Face::XLO][n]))
                    in(i,j,k,n) = -in(1-glevel[0],j,k,n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::XLO][n]) || (nodal && i == lo.x)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[0]>0 || (nodal && i == hi.x)) && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::XHI][n]))
                {
                    const Set::Scalar val = m_bc_func[Face::XHI][n](m_geom.ProbHi()[0],y,z,t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[0]>0) in(i,j,k,n) = 2.0 * val - in(2 * hi.x - i + 1,j,k,n);
                }
                else if(glevel[0]>0 && BCUtil::IsNeumann(m_bc_type[Face::XHI][n]))
                    in(i,j,k,n) = in(i-glevel[0],j,k,n) - (m_bc_func[Face::XHI].size() > 0 ? m_bc_func[Face::XHI][n](x,y,z,t)*DX[0] : 0);
                else if(glevel[0]>0 && BCUtil::IsReflectEven(m_bc_type[Face::XHI][n]))
                    in(i,j,k,n) = in(hi.x-glevel[0],j,k,n);
                else if(glevel[0]>0 && BCUtil::IsReflectOdd(m_bc_type[Face::XHI][n]))
                    in(i,j,k,n) = -in(hi.x-glevel[0],j,k,n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::XHI][n]) || (nodal && i == hi.x)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[1]<0 || (nodal && j == lo.y)) && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::YLO][n]))
                {
                    const Set::Scalar val = m_bc_func[Face::YLO][n](x,m_geom.ProbLo()[1],z,t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[1]<0) in(i,j,k,n) = 2.0 * val - in(i,2 * lo.y - j - 1,k,n);
                }
                else if (glevel[1]<0 && BCUtil::IsNeumann(m_bc_type[Face::YLO][n]))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n) - (m_bc_func[Face::YLO].size() > 0 ? m_bc_func[Face::YLO][n](x,y,z,t)*DX[1] : 0);
                else if (glevel[1]<0 && BCUtil::IsReflectEven(m_bc_type[Face::YLO][n]))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n);
                else if (glevel[1]<0 && BCUtil::IsReflectOdd(m_bc_type[Face::YLO][n]))
                    in(i,j,k,n) = -in(i,j-glevel[1],k,n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::YLO][n]) || (nodal && j == lo.y)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[1]>0 || (nodal && j == hi.y)) && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::YHI][n]))
                {
                    const Set::Scalar val = m_bc_func[Face::YHI][n](x,m_geom.ProbHi()[1],z,t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[1]>0) in(i,j,k,n) = 2.0 * val - in(i,2 * hi.y - j + 1,k,n);
                }
                else if (glevel[1]>0 && BCUtil::IsNeumann(m_bc_type[Face::YHI][n]))
                    in(i,j,k,n) = in(i,j-glevel[1],k,n) - (m_bc_func[Face::YHI].size() > 0 ? m_bc_func[Face::YHI][n](x,y,z,t)*DX[1] : 0);
                else if (glevel[1]>0 && BCUtil::IsReflectEven(m_bc_type[Face::YHI][n]))
                    in(i,j,k,n) = in(i,hi.y-glevel[1],k,n);
                else if (glevel[1]>0 && BCUtil::IsReflectOdd(m_bc_type[Face::YHI][n]))
                    in(i,j,k,n) = -in(i,hi.y-glevel[1],k,n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::YHI][n]) || (nodal && j == hi.y)) {}
                else
                    Util::Abort(INFO, "Incorrect boundary conditions");
            }
#if AMREX_SPACEDIM>2
            else if ((glevel[2]<0 || (nodal && k == lo.z)) && (face == Orientation::zlo || face == Orientation::All))
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::ZLO][n]))
                {
                    const Set::Scalar val = m_bc_func[Face::ZLO][n](x,y,m_geom.ProbLo()[2],t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[2]<0) in(i,j,k,n) = 2.0 * val - in(i,j,2 * lo.z - k - 1,n);
                }
                else if (glevel[2]<0 && BCUtil::IsNeumann(m_bc_type[Face::ZLO][n]))
                    in(i,j,k,n) = in(i,j,k-glevel[2],n) - (m_bc_func[Face::ZLO].size() > 0 ? m_bc_func[Face::ZLO][n](x,y,z,t)*DX[2] : 0);
                else if (glevel[2]<0 && BCUtil::IsReflectEven(m_bc_type[Face::ZLO][n]))
                    in(i,j,k,n) = in(i,j,1-glevel[2],n);
                else if (glevel[2]<0 && BCUtil::IsReflectOdd(m_bc_type[Face::ZLO][n]))
                    in(i,j,k,n) = -in(i,j,1-glevel[2],n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::ZLO][n]) || (nodal && k == lo.z)) {}
                else Util::Abort(INFO, "Incorrect boundary conditions");
            }
            else if ((glevel[2]>0 || (nodal && k == hi.z)) && (face == Orientation::zhi || face == Orientation::All))
            {
                if (BCUtil::IsDirichlet(m_bc_type[Face::ZHI][n]))
                {
                    const Set::Scalar val = m_bc_func[Face::ZHI][n](x,y,m_geom.ProbHi()[2],t);
                    if (nodal) in(i,j,k,n) = val;
                    else if (glevel[2]>0) in(i,j,k,n) = 2.0 * val - in(i,j,2 * hi.z - k + 1,n);
                }
                else if(glevel[2]>0 && BCUtil::IsNeumann(m_bc_type[Face::ZHI][n]))
                    in(i,j,k,n) = in(i,j,k-glevel[2],n) - (m_bc_func[Face::ZHI].size() > 0 ? m_bc_func[Face::ZHI][n](x,y,z,t)*DX[2] : 0);
                else if(glevel[2]>0 && BCUtil::IsReflectEven(m_bc_type[Face::ZHI][n]))
                    in(i,j,k,n) = in(i,j,hi.z-glevel[2],n);
                else if(glevel[2]>0 && BCUtil::IsReflectOdd(m_bc_type[Face::ZHI][n]))
                    in(i,j,k,n) = -in(i,j,hi.z-glevel[2],n);
                else if(BCUtil::IsPeriodic(m_bc_type[Face::ZHI][n]) || (nodal && k == hi.z)) {}
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
    int bc_lo[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XLO][component],m_bc_type[Face::YLO][component],m_bc_type[Face::XLO][component])};
    int bc_hi[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XHI][component],m_bc_type[Face::YHI][component],m_bc_type[Face::XHI][component])};

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
